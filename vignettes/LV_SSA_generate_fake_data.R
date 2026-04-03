# Lotka–Volterra fake observations via Gillespie SSA (source: R/SSA.R).
#
# Parameter choices (aligned with source/LV_stochastic_data_generator.R):
#   - alpha = 1, beta = 1, gamma = 1, delta = 1: symmetric rates so mean-field
#     matches the workshop ODE scripts (a = 1, b = 1 with beta, gamma fixed at 1).
#   - N = 100: system-size scaling; integer copy numbers are N times the densities
#     (x0 = 1, y0 = 0.5), giving reachable states without excessive demographic noise.
#   - Time grid: 10 equally spaced points from 0 to 14 (endpoint included), comparable
#     to t_max = 15 used elsewhere while keeping a round number of observation times.
#   - set.seed(42): reproducible single trajectory for teaching / debugging.
#
# Output: data/lv_ssa_fake_counts.csv with columns time, rabbit, fox (integer counts).

resolve_repo_root <- function() {
    ca <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", ca[grep("^--file=", ca)])
    if (length(f)) {
        return(normalizePath(file.path(dirname(normalizePath(f[1])), "..")))
    }
    if (file.exists(file.path(getwd(), "R", "SSA.R"))) {
        return(normalizePath(getwd()))
    }
    up <- normalizePath(file.path(getwd(), ".."))
    if (file.exists(file.path(up, "R", "SSA.R"))) {
        return(up)
    }
    normalizePath(getwd())
}

root <- resolve_repo_root()
source(file.path(root, "R", "SSA.R"))

N <- 1000 #
parameters <- list(
    alpha = 1,
    beta  = 1,
    gamma = 1,
    delta = 1,
    N     = N
)

initial_state <- list(
    rabbit = as.integer(max(1L, round(0.7 * N))),
    fox    = as.integer(max(0L, round(0.3 * N)))
)

time_points <- seq(0, 14, length.out = 15L)

reaction_propensities <- function(state, parms) {
    with(c(state, parms), c(
        alpha * rabbit,
        (beta  / N) * rabbit * fox,
        (delta / N) * rabbit * fox,
        gamma * fox
    ))
}

reaction_stoichiometries <- list(
    list(rabbit =  1L, fox =  0L),
    list(rabbit = -1L, fox =  0L),
    list(rabbit =  0L, fox =  1L),
    list(rabbit =  0L, fox = -1L)
)

set.seed(42)
raw <- SSA(
    initial_state            = initial_state,
    parameters               = parameters,
    reaction_propensities    = reaction_propensities,
    reaction_stoichiometries = reaction_stoichiometries,
    time_points              = time_points
)

out_df <- data.frame(
    time   = time_points,
    rabbit = as.integer(round(raw[, "rabbit"])),
    fox    = as.integer(round(raw[, "fox"]))
)

data_dir <- file.path(root, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(data_dir, "lv_ssa_fake_counts.csv")
write.csv(out_df, file = out_path, row.names = FALSE)

meta <- data.frame(
    alpha = parameters$alpha,
    beta  = parameters$beta,
    gamma = parameters$gamma,
    delta = parameters$delta,
    N     = N
)
meta_path <- file.path(data_dir, "lv_ssa_metadata.csv")
write.csv(meta, file = meta_path, row.names = FALSE)

message("Wrote ", out_path)
message("Wrote ", meta_path)
print(out_df)
