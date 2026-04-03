# Tau-leaping approximate stochastic simulation.
# Drop-in replacement for SSA() with an additional tau parameter.
# Faster but approximate: each step draws Poisson-distributed reaction
# counts instead of simulating individual events.

tau_leaping <- function(initial_state, parameters, reaction_propensities, reaction_stoichiometries, time_points, tau = 0.01, max_pop = 5000L) {
    stopifnot(all(diff(time_points) > 0))
    state <- initial_state
    time <- 0
    n_reactions <- length(reaction_stoichiometries)
    state_out <- matrix(0, nrow = length(time_points), ncol = length(initial_state))
    colnames(state_out) <- names(initial_state)
    next_time_point_index <- 1L

    while (next_time_point_index <= length(time_points) &&
        time_points[next_time_point_index] <= time) {
        state_out[next_time_point_index, ] <- unlist(state)
        next_time_point_index <- next_time_point_index + 1L
    }

    while (time < time_points[length(time_points)]) {
        props <- reaction_propensities(state, parameters)
        total_prop <- sum(props)

        if (total_prop == 0) {
            while (next_time_point_index <= length(time_points)) {
                state_out[next_time_point_index, ] <- unlist(state)
                next_time_point_index <- next_time_point_index + 1L
            }
            break
        }

        firings <- rpois(n_reactions, pmax(props, 0) * tau)

        for (r in seq_len(n_reactions)) {
            if (firings[r] > 0) {
                state <- mapply(
                    function(x, y) x + y * firings[r],
                    state, reaction_stoichiometries[[r]],
                    SIMPLIFY = FALSE
                )
            }
        }

        state <- lapply(state, function(x) min(max(x, 0), max_pop))

        time_next <- time + tau
        while (next_time_point_index <= length(time_points) &&
            time < time_points[next_time_point_index] &&
            time_next >= time_points[next_time_point_index]) {
            state_out[next_time_point_index, ] <- unlist(state)
            next_time_point_index <- next_time_point_index + 1L
        }
        time <- time_next
    }
    return(state_out)
}
