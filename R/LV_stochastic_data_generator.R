# =============================================================================
# Lotka-Volterra Stochastic Data Generator
# -----------------------------------------------------------------------------
# Implements a Gillespie-type Stochastic Simulation Algorithm (SSA) for the
# Lotka-Volterra predator-prey system.  This script is intended to generate
# "observed" (fake) data for use as the inference target in ABC / MCMC
# workshops.
#
# Model:
#   dx/dt = alpha * x - beta  * x * y      (prey:     birth - predation loss)
#   dy/dt = delta * x * y    - gamma * y   (predator: predation gain - death)
#
# Default parameterisation matches the ODE reference script:
#   alpha = a  (to infer),  beta  = 1 (fixed)
#   delta = b  (to infer),  gamma = 1 (fixed)
#
# Volume / system-size scaling (N):
#   Integer copy numbers  X = N * x,  Y = N * y  are propagated by SSA.
#   At the end, X/N and Y/N are returned so that output is on the same
#   density scale as the ODE model.
#
# Stochastic reactions (4-reaction formulation):
#   R1:  X  ->  X+1          propensity = alpha * X         (prey  birth)
#   R2:  X  ->  X-1          propensity = (beta/N)  * X * Y (prey  death by predation)
#   R3:  Y  ->  Y+1          propensity = (delta/N) * X * Y (predator birth from predation)
#   R4:  Y  ->  Y-1          propensity = gamma * Y         (predator death)
# =============================================================================

# -------------------------------------------------------------------------
# SSA — Gillespie Direct Method
# -------------------------------------------------------------------------
#   initial_state          : named list of integer copy numbers
#   parameters             : named list of rate constants
#   reaction_propensities  : function(state, parameters) -> numeric vector
#   reaction_stoichiometries: list of named lists (stoichiometric changes)
#   time_points            : strictly increasing numeric vector of output times
#
#   Returns a matrix [length(time_points) x length(initial_state)] with the
#   state recorded at each requested time point (last-value carried forward).
# -------------------------------------------------------------------------
SSA_LV <- function(initial_state, parameters,
                   reaction_propensities, reaction_stoichiometries,
                   time_points) {
    stopifnot(all(diff(time_points) > 0))

    state               <- initial_state
    time                <- 0
    n_species           <- length(initial_state)
    state_out           <- matrix(0, nrow = length(time_points), ncol = n_species)
    colnames(state_out) <- names(initial_state)
    next_tp_idx         <- 1

    while (time < time_points[length(time_points)]) {
        props      <- reaction_propensities(state, parameters)
        total_prop <- sum(props)

        if (total_prop == 0) {
            # Absorbing state: fill remaining time points with current state
            while (next_tp_idx <= length(time_points)) {
                state_out[next_tp_idx, ] <- unlist(state)
                next_tp_idx <- next_tp_idx + 1
            }
            break
        }

        delta_time <- rexp(1, rate = total_prop)
        reaction   <- sample.int(length(props), size = 1, prob = props / total_prop)
        state      <- mapply(
            function(s, delta) s + delta,
            state, reaction_stoichiometries[[reaction]],
            SIMPLIFY = FALSE
        )

        time_next <- time + delta_time
        # Record state for every time point the trajectory passes through
        while (next_tp_idx <= length(time_points) &&
               time  <  time_points[next_tp_idx] &&
               time_next >= time_points[next_tp_idx]) {
            state_out[next_tp_idx, ] <- unlist(state)
            next_tp_idx <- next_tp_idx + 1
        }
        time <- time_next
    }

    return(state_out)
}

# -------------------------------------------------------------------------
# LV_stochastic — single-run stochastic simulation
# -------------------------------------------------------------------------
#   a      : prey growth rate alpha
#   b      : predator efficiency delta
#   N      : system size (integer copy numbers = N * density)
#   t_obs  : observation time points (must be < t_max)
#   t_max  : end time of simulation
#
#   Returns a one-row data.frame with columns:
#     a, b,  X_1 … X_k,  Y_1 … Y_k
#   where X_i / Y_i are density values (copy number / N) at t_obs[i].
# -------------------------------------------------------------------------
LV_stochastic <- function(a, b, N = 100,
                           t_obs = c(1.1, 2.4, 3.9, 5.6, 7.5, 9.6, 11.9, 14.4),
                           t_max = 15) {
    stopifnot(N > 0, all(t_obs < t_max))

    initial_state <- list(
        X = round(1.0 * N),   # prey  density 1.0
        Y = round(0.5 * N)    # predator density 0.5
    )

    parameters <- list(
        alpha = a,
        beta  = 1,
        delta = b,
        gamma = 1,
        N     = N
    )

    # 4-reaction Lotka-Volterra
    reaction_propensities <- function(state, parms) {
        with(c(state, parms), c(
            alpha * X,              # R1: prey birth
            (beta  / N) * X * Y,   # R2: prey death by predation
            (delta / N) * X * Y,   # R3: predator birth from predation
            gamma * Y               # R4: predator death
        ))
    }

    reaction_stoichiometries <- list(
        list(X =  1L, Y =  0L),    # R1
        list(X = -1L, Y =  0L),    # R2
        list(X =  0L, Y =  1L),    # R3
        list(X =  0L, Y = -1L)     # R4
    )

    raw <- SSA_LV(
        initial_state           = initial_state,
        parameters              = parameters,
        reaction_propensities   = reaction_propensities,
        reaction_stoichiometries = reaction_stoichiometries,
        time_points             = t_obs
    )

    # Convert copy numbers back to density scale
    X_obs <- raw[, "X"] / N
    Y_obs <- raw[, "Y"] / N

    n_obs  <- length(t_obs)
    result <- data.frame(matrix(c(a, b, X_obs, Y_obs), nrow = 1))
    colnames(result) <- c("a", "b",
                          paste0("X_", seq_len(n_obs)),
                          paste0("Y_", seq_len(n_obs)))
    return(result)
}

# -------------------------------------------------------------------------
# generate_fake_data — wrapper for one or more parameter sets
# -------------------------------------------------------------------------
#   parameters : data.frame with columns "a" and "b" (one row per run)
#   N          : system size passed to LV_stochastic
#   seed       : optional RNG seed for reproducibility
#
#   Returns a data.frame in the same column layout as LV_stochastic.
# -------------------------------------------------------------------------
generate_fake_data <- function(parameters, N = 100, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    stopifnot(all(c("a", "b") %in% colnames(parameters)))

    results <- lapply(seq_len(nrow(parameters)), function(i) {
        LV_stochastic(a = parameters$a[i], b = parameters$b[i], N = N)
    })
    do.call(rbind, results)
}

# =============================================================================
# Demo — generate one "observed" dataset at the true parameters
# =============================================================================
if (sys.nframe() == 0) {                    # only runs when sourced directly

    parameters_truth <- data.frame(a = 1, b = 1)

    set.seed(42)
    observed <- generate_fake_data(parameters_truth, N = 100)

    cat("=== Stochastic observation at true parameters (a=1, b=1) ===\n")
    Xt <- unlist(observed[, paste0("X_", 1:8)])
    Yt <- unlist(observed[, paste0("Y_", 1:8)])
    cat("Xt:", paste(round(Xt, 3), collapse = ", "), "\n")
    cat("Yt:", paste(round(Yt, 3), collapse = ", "), "\n\n")

    # Summary statistics target (drop parameter columns)
    statistics_target <- observed[, -(1:2)]
    cat("statistics_target ready for inference.\n")

    # Optional: plot trajectories
    if (requireNamespace("ggplot2", quietly = TRUE)) {
        library(ggplot2)
        t_obs <- c(1.1, 2.4, 3.9, 5.6, 7.5, 9.6, 11.9, 14.4)
        df_plot <- data.frame(
            time    = rep(t_obs, 2),
            density = c(Xt, Yt),
            species = rep(c("Prey (X)", "Predator (Y)"), each = 8)
        )
        p <- ggplot(df_plot, aes(x = time, y = density,
                                 colour = species, shape = species)) +
            geom_line(linewidth = 0.8) +
            geom_point(size = 3) +
            labs(title = "Lotka-Volterra: stochastic fake data  (N = 100)",
                 x = "Time", y = "Population density",
                 colour = NULL, shape = NULL) +
            theme_bw(base_size = 14) +
            theme(legend.position = "top")
        print(p)
    }
}
