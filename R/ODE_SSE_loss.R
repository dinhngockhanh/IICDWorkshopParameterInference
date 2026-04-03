# Sum of squared errors (SSE) between observed counts and Lotka-Volterra ODE fit.
# Uses lv_populations_at_time() from R/ODE_solver_at_time.R; initial state is taken
# from the data row with time == 0 (within numerical tolerance).

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
# data               : data.frame with columns time, rabbit, fox (e.g. lv_ssa_fake_counts.csv)
# metadata           : one-row data.frame from lv_ssa_metadata.csv (must contain N);
#                      if NULL the function falls back to N inside parameters_estimate, then 1
# Returns: scalar SSE = sum_t (rabbit_obs - rabbit_pred)^2 + (fox_obs - fox_pred)^2
ODE_SSE_loss <- function(parameters_estimate, data, metadata = NULL) {
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
    sse <- sum((dat$rabbit - pred$rabbit)^2 + (dat$fox - pred$fox)^2)
    as.numeric(sse)
}
