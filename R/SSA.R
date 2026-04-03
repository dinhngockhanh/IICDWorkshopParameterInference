# Gillespie Stochastic Simulation Algorithm (SSA).
# Extracted for reuse across workshop models.

SSA <- function(initial_state, parameters, reaction_propensities, reaction_stoichiometries, time_points) {
    stopifnot(all(diff(time_points) > 0))
    state <- initial_state
    time <- 0
    state_out <- matrix(0, nrow = length(time_points), ncol = length(initial_state))
    colnames(state_out) <- names(initial_state)
    next_time_point_index <- 1L
    # Record any output times <= initial time (typically t = 0)
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

        delta_time <- rexp(1, rate = total_prop)
        reaction <- sample.int(length(props), size = 1, prob = props / total_prop)
        state <- mapply(function(x, y) x + y, state, reaction_stoichiometries[[reaction]], SIMPLIFY = FALSE)

        time_next <- time + delta_time
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
