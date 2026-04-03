# Gaussian observation model; returns log-likelihood log L (sum of Normal log-densities).
# Uses lv_populations_at_time() from R/ODE_solver_at_time.R; initial state at t = 0
# matches ODE_SSE_loss.R data handling; numerically stable in log-space.

.local_r_dir <- function() {
    w <- getwd()
    if (file.exists(file.path(w, "R", "ODE_solver_at_time.R"))) {
        return(normalizePath(file.path(w, "R")))
    }
    if (file.exists(file.path(w, "ODE_solver_at_time.R"))) {
        return(normalizePath(w))
    }
    cmd <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", cmd[grep("^--file=", cmd)])
    if (length(f)) {
        d <- dirname(normalizePath(f[1]))
        if (file.exists(file.path(d, "ODE_solver_at_time.R"))) {
            return(d)
        }
    }
    stop(
        "Cannot find ODE_solver_at_time.R. ",
        "Use the repository root or R/ as working directory, or run via Rscript on this file."
    )
}

source(file.path(.local_r_dir(), "ODE_solver_at_time.R"))

.as_param_list <- function(parameters_estimate) {
    if (is.data.frame(parameters_estimate)) {
        stopifnot(nrow(parameters_estimate) == 1L)
        parameters_estimate <- as.list(parameters_estimate[1, , drop = TRUE])
    }
    stopifnot(is.list(parameters_estimate))
    nm <- c("alpha", "beta", "gamma", "delta")
    stopifnot(all(nm %in% names(parameters_estimate)))
    out <- list(
        alpha = as.numeric(parameters_estimate$alpha),
        beta  = as.numeric(parameters_estimate$beta),
        gamma = as.numeric(parameters_estimate$gamma),
        delta = as.numeric(parameters_estimate$delta)
    )
    if (!is.null(parameters_estimate$N)) {
        out$N <- as.numeric(parameters_estimate$N)
    }
    out
}

# parameters_estimate : list or one-row data.frame with alpha, beta, gamma, delta (and optionally N)
# data               : data.frame with columns time, rabbit, fox
# variance            : scalar variance sigma^2 for each component (default 200^2);
#                       rabbit and fox counts at each time assumed conditionally independent
# metadata           : one-row data.frame with N for ODE scaling when N missing in parameters
#
# Model: for each time t, rabbit_obs(t), fox_obs(t) | mu(t) ~ Normal(ODE mean, variance) i.i.d. components
# Returns: scalar log-likelihood log L (same Gaussian model, all constant terms included).
ODE_Gaussian_loglikelihood <- function(parameters_estimate, data, variance = 200^2, metadata = NULL) {
    stopifnot(length(variance) == 1L, is.finite(variance), variance > 0)
    parms <- .as_param_list(parameters_estimate)
    if (is.null(parms$N) && !is.null(metadata)) {
        parms$N <- as.numeric(metadata$N[1])
    }
    if (is.null(parms$N)) parms$N <- 1

    stopifnot(all(c("time", "rabbit", "fox") %in% names(data)))
    dat <- data[, c("time", "rabbit", "fox")]
    dat <- dat[order(dat$time), , drop = FALSE]
    rownames(dat) <- NULL

    t0 <- dat$time
    i0 <- which.min(abs(t0 - 0))
    stopifnot(abs(t0[i0]) < 1e-8)
    initial_state <- list(
        rabbit = as.numeric(dat$rabbit[i0]),
        fox    = as.numeric(dat$fox[i0])
    )

    pred <- lv_populations_at_time(parms, times = t0, initial_state = initial_state)
    res_r <- dat$rabbit - pred$rabbit
    res_f <- dat$fox - pred$fox
    sse_part <- sum(res_r^2 + res_f^2)
    n <- 2L * nrow(dat)
    log_lik <- -(n / 2) * log(2 * pi * variance) - sse_part / (2 * variance)
    as.numeric(log_lik)
}
