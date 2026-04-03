# Lotka-Volterra: Metropolis–Hastings MCMC targeting π(θ|y) ∝ L(y|θ) π(θ)
# with Gaussian observation model (log-likelihood via R/ODE_Gaussian_loglikelihood.R)
# and the same uniform priors as vignettes/LV_random_search.R.
#
# Parallel chains: uses (detectCores() - 1) workers (minimum 1 chain). Each chain has its
# own initial draw from the prior and its own RNG stream (seed offset per chain).
#
# Run from repository root:
#   Rscript vignettes/LV_MCMC.R
#
# Output: vignettes/MCMC/ (combined chain CSV, per-chain acceptance, trace/histogram plots).

resolve_repo_root <- function() {
    ca <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", ca[grep("^--file=", ca)])
    if (length(f)) {
        return(normalizePath(file.path(dirname(normalizePath(f[1])), "..")))
    }
    if (file.exists(file.path(getwd(), "R", "ODE_Gaussian_loglikelihood.R"))) {
        return(normalizePath(getwd()))
    }
    up <- normalizePath(file.path(getwd(), ".."))
    if (file.exists(file.path(up, "R", "ODE_Gaussian_loglikelihood.R"))) {
        return(up)
    }
    normalizePath(getwd())
}

root <- resolve_repo_root()
source(file.path(root, "R", "ODE_Gaussian_loglikelihood.R"))
suppressPackageStartupMessages(library(parallel))

# ----- 1. Priors (match LV_random_search.R) -----
# Independent uniforms:  theta_j ~ Uniform(prior_min[j], prior_max[j]), j = alpha…delta.
prior_min <- c(alpha = 0, beta = 0, gamma = 0, delta = 0)
prior_max <- c(alpha = 2, beta = 2, gamma = 2, delta = 2)
stopifnot(all(prior_max > prior_min))
par_names <- names(prior_min)

log_prior_uniform_box <- function(theta, lo = prior_min, hi = prior_max) {
    theta <- as.numeric(theta)
    if (length(theta) != length(lo)) stop("theta length must match prior dimension.")
    if (any(theta < lo | theta > hi, na.rm = TRUE) || anyNA(theta)) {
        return(-Inf)
    }
    sum(-log(hi - lo))
}

# ----- Modular MH core (log-space, symmetric Gaussian random walk) -----
mcmc_metropolis_hastings <- function(log_lik, log_pi, theta_init, n_iter,
                                     proposal_sd = 1) {
    nm_init <- base::names(theta_init)
    theta_init <- as.numeric(theta_init)
    if (!is.null(nm_init) && length(nm_init) == length(theta_init)) {
        base::names(theta_init) <- nm_init
    }
    p <- length(theta_init)
    stopifnot(n_iter >= 1L, length(proposal_sd) == 1L, proposal_sd > 0, p >= 1L)

    chain <- matrix(NA_real_, nrow = n_iter, ncol = p)
    colnames(chain) <- base::names(theta_init)
    if (is.null(colnames(chain))) colnames(chain) <- paste0("p", seq_len(p))

    logpost <- function(th) log_lik(th) + log_pi(th)
    lp_cur <- logpost(theta_init)
    if (!is.finite(lp_cur)) {
        stop("Initial theta has non-finite log posterior; choose theta_init inside support.")
    }

    cur <- theta_init
    n_accept <- 0L
    for (t in seq_len(n_iter)) {
        prop <- cur + stats::rnorm(p, mean = 0, sd = proposal_sd)
        lp_prop <- logpost(prop)
        log_alpha <- min(0, lp_prop - lp_cur)
        if (log(stats::runif(1L)) < log_alpha) {
            cur <- prop
            lp_cur <- lp_prop
            n_accept <- n_accept + 1L
        }
        chain[t, ] <- cur
    }

    list(
        chain            = chain,
        acceptance_rate  = n_accept / n_iter,
        n_iter           = n_iter,
        proposal_sd      = proposal_sd,
        n_accept         = n_accept
    )
}

drop_burnin <- function(chain, n_burn) {
    n_burn <- as.integer(n_burn)
    stopifnot(n_burn >= 0L, n_burn < nrow(chain))
    if (n_burn == 0L) return(chain)
    chain[(n_burn + 1L):nrow(chain), , drop = FALSE]
}

# ----- Data -----
data_path <- file.path(root, "data", "lv_ssa_fake_counts.csv")
meta_path <- file.path(root, "data", "lv_ssa_metadata.csv")
stopifnot(file.exists(data_path), file.exists(meta_path))
obs <- read.csv(data_path, check.names = FALSE)
meta <- read.csv(meta_path, check.names = FALSE)

variance <- 200^2

log_lik_lv <- function(theta) {
    p <- stats::setNames(as.numeric(theta), par_names)
    ll <- tryCatch(
        ODE_Gaussian_loglikelihood(as.list(p), obs, variance = variance, metadata = meta),
        error = function(e) -Inf
    )
    if (!is.finite(ll)) -Inf else ll
}

log_pi_lv <- function(theta) log_prior_uniform_box(theta)

# ----- Parallelism: number of chains = max(1, n_cores - 1) -----
n_cores_detected <- parallel::detectCores()
if (is.na(n_cores_detected) || n_cores_detected < 1L) n_cores_detected <- 1L
n_chains <- max(1L, as.integer(n_cores_detected - 1L))
message(
    "Detected cores: ", n_cores_detected,
    " | parallel MCMC chains: ", n_chains,
    " (workers = max(1, cores - 1))"
)

# ----- MCMC settings -----
proposal_sd <- 1
n_iter <- 10000L
base_seed <- 2L

# One chain: initial draw from prior; parallel: each worker draws independently
run_single_chain <- function(chain_id, seed) {
    set.seed(seed)
    theta_init <- stats::setNames(
        mapply(function(lo, hi) stats::runif(1L, min = lo, max = hi), prior_min, prior_max),
        par_names
    )
    fit <- mcmc_metropolis_hastings(
        log_lik     = log_lik_lv,
        log_pi      = log_pi_lv,
        theta_init  = theta_init,
        n_iter      = n_iter,
        proposal_sd = proposal_sd
    )
    fit$chain_id <- chain_id
    fit$theta_init <- theta_init
    fit
}

# ----- Run chains in parallel (PSOCK: works on macOS / Windows / Linux) -----
chain_seeds <- base_seed + seq_len(n_chains) * 100000L
t_run <- proc.time()

if (n_chains == 1L) {
    results <- list(run_single_chain(1L, chain_seeds[1L]))
} else {
    cl <- parallel::makeCluster(n_chains)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    ge <- .GlobalEnv
    parallel::clusterExport(cl, "root", envir = ge)
    parallel::clusterEvalQ(cl, {
        source(file.path(root, "R", "ODE_Gaussian_loglikelihood.R"))
    })
    varlist <- c(
        "run_single_chain", "mcmc_metropolis_hastings",
        "log_lik_lv", "log_pi_lv", "log_prior_uniform_box",
        "prior_min", "prior_max", "par_names",
        "obs", "meta", "variance", "n_iter", "proposal_sd", "chain_seeds"
    )
    parallel::clusterExport(cl, varlist, envir = ge)
    results <- parallel::parLapply(cl, seq_len(n_chains), function(k) {
        run_single_chain(k, chain_seeds[k])
    })
}

runtime_sec <- as.numeric((proc.time() - t_run)["elapsed"])

# ----- Combine chains for storage / plots -----
chain_rows <- do.call(rbind, lapply(results, function(fit) {
    ch <- fit$chain
    data.frame(
        chain     = fit$chain_id,
        iteration = seq_len(nrow(ch)),
        as.data.frame(ch),
        stringsAsFactors = FALSE
    )
}))

out_dir <- file.path(root, "vignettes", "MCMC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

chain_path <- file.path(out_dir, "lv_mcmc_chains.csv")
utils::write.csv(chain_rows, file = chain_path, row.names = FALSE)

acc_df <- data.frame(
    chain           = vapply(results, function(z) z$chain_id, integer(1)),
    acceptance_rate = vapply(results, function(z) z$acceptance_rate, numeric(1)),
    n_accept        = vapply(results, function(z) z$n_accept, integer(1)),
    n_iter          = vapply(results, function(z) z$n_iter, integer(1)),
    proposal_sd     = vapply(results, function(z) z$proposal_sd, numeric(1))
)
acc_df$runtime_sec <- NA_real_
acc_df$runtime_sec[1] <- runtime_sec
utils::write.csv(acc_df, file = file.path(out_dir, "lv_mcmc_acceptance.csv"), row.names = FALSE)

message(
    "MCMC done. Combined chain: ", chain_path,
    " | mean acceptance = ", signif(mean(acc_df$acceptance_rate), 4),
    " | runtime_sec = ", signif(runtime_sec, 5)
)

# ----- Plots: colour by chain; facets by parameter -----
n_burn <- min(200L, as.integer(n_iter %/% 5L))

if (requireNamespace("ggplot2", quietly = TRUE)) {
    library(ggplot2)
    n_ch <- length(results)
    long_tr <- do.call(rbind, lapply(results, function(fit) {
        ch <- fit$chain
        nr <- nrow(ch)
        nc <- ncol(ch)
        do.call(rbind, lapply(seq_len(nc), function(j) {
            data.frame(
                chain     = fit$chain_id,
                iteration = seq_len(nr),
                parameter = colnames(ch)[j],
                value     = ch[, j],
                stringsAsFactors = FALSE
            )
        }))
    }))

    p_trace <- ggplot(long_tr, aes(x = iteration, y = value, colour = factor(chain))) +
        geom_line(linewidth = 0.25, alpha = 0.85) +
        facet_wrap(~parameter, scales = "free_y", ncol = 1) +
        labs(
            title = paste0("MCMC traces (", n_ch, " chains, coloured by chain)"),
            x = "Iteration", y = NULL, colour = "chain"
        ) +
        theme_bw(base_size = 11) +
        theme(legend.position = "bottom")

    long_hi <- do.call(rbind, lapply(results, function(fit) {
        ch <- drop_burnin(fit$chain, n_burn)
        nr <- nrow(ch)
        nc <- ncol(ch)
        do.call(rbind, lapply(seq_len(nc), function(j) {
            data.frame(
                chain     = fit$chain_id,
                parameter = colnames(ch)[j],
                value     = ch[, j],
                stringsAsFactors = FALSE
            )
        }))
    }))

    p_hist <- ggplot(long_hi, aes(x = value, fill = factor(chain))) +
        geom_histogram(position = "identity", alpha = 0.45, bins = 40, colour = NA) +
        facet_wrap(~parameter, scales = "free", ncol = 2) +
        labs(
            title = paste0("Posterior marginals (burn-in = ", n_burn, " per chain)"),
            x = NULL, y = "Count", fill = "chain"
        ) +
        theme_bw(base_size = 11) +
        theme(legend.position = "bottom")

    ggsave(file.path(out_dir, "lv_mcmc_trace.png"), p_trace, width = 9, height = 10, dpi = 150)
    ggsave(file.path(out_dir, "lv_mcmc_hist.png"), p_hist, width = 9, height = 6, dpi = 150)

    # Optional: overlaid density (post burn-in), easier to compare chains
    p_den <- ggplot(long_hi, aes(x = value, colour = factor(chain))) +
        geom_density(linewidth = 0.6, adjust = 1.1) +
        facet_wrap(~parameter, scales = "free", ncol = 2) +
        labs(
            title = paste0("Posterior marginals — density (burn-in = ", n_burn, ")"),
            x = NULL, y = "Density", colour = "chain"
        ) +
        theme_bw(base_size = 11) +
        theme(legend.position = "bottom")
    ggsave(file.path(out_dir, "lv_mcmc_density.png"), p_den, width = 9, height = 6, dpi = 150)

    message("Saved plots: lv_mcmc_trace.png, lv_mcmc_hist.png, lv_mcmc_density.png in ", out_dir)
} else {
    message("Package ggplot2 not installed; skip plots.")
}

# =============================================================================
# Example stub (toy Gaussian): not run
# =============================================================================
if (FALSE) {
    set.seed(1)
    y_obs <- stats::rnorm(20, mean = 2.5, sd = 1)
    log_lik_toy <- function(theta) sum(stats::dnorm(y_obs, mean = theta, sd = 1, log = TRUE))
    log_pi_toy <- function(theta) stats::dnorm(theta, mean = 0, sd = 10, log = TRUE)
    fit_toy <- mcmc_metropolis_hastings(
        log_lik_toy,
        log_pi_toy,
        theta_init = 0,
        n_iter = 3000L,
        proposal_sd = 1
    )
    message("Toy example acceptance rate = ", fit_toy$acceptance_rate)
}
