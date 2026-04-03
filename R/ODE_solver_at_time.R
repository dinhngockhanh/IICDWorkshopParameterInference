library(deSolve)

# Forward Lotka–Volterra ODE at requested times.
#   parameters     : list with alpha, beta, gamma, delta, N
#   times          : numeric vector of times >= 0 (initial state is at t = 0)
#   initial_state  : list with rabbit (prey x) and fox (predator y) at t = 0
#   Returns data.frame with columns time, rabbit, fox (rows follow order of times).
#
#   N is the system-size parameter. The interaction terms are scaled by 1/N so
#   that state variables can be integer counts consistent with the SSA.
#   When N = 1 the equations reduce to the standard density form.
lv_populations_at_time <- function(parameters, times, initial_state) {
    alpha <- parameters$alpha
    beta <- parameters$beta
    gamma <- parameters$gamma
    delta <- parameters$delta
    N <- if (!is.null(parameters$N)) parameters$N else 1
    x0 <- initial_state$rabbit
    y0 <- initial_state$fox

    stopifnot(is.numeric(times), all(times >= 0))

    LV_rhs <- function(Time, State, parms) {
        with(as.list(c(State, parms)), {
            dx <- x * (alpha - (beta / N) * y)
            dy <- -y * (gamma - (delta / N) * x)
            list(c(dx, dy))
        })
    }
    parms <- c(alpha = alpha, beta = beta, gamma = gamma, delta = delta, N = N)
    y0v <- c(x = x0, y = y0)

    t_ode <- sort(unique(c(0, as.numeric(times))))
    if (length(t_ode) < 2L) {
        eps <- 1e-12 * (1 + abs(t_ode[1]))
        t_ode <- c(t_ode[1], t_ode[1] + max(eps, 1e-15))
    }
    out <- as.data.frame(
    ode(
        func = LV_rhs,
        y = y0v,
        parms = parms,
        times = t_ode,
        rtol = 1e-3,   # default ~1e-6
        atol = 1e-6    # default ~1e-8
    )
    )

    idx <- vapply(times, function(t) which.min(abs(out$time - t)), integer(1))
    data.frame(
        time = times,
        rabbit = as.numeric(out$x[idx]),
        fox = as.numeric(out$y[idx])
    )
}
