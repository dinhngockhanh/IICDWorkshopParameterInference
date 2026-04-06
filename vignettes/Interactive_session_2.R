library(IICDWorkshopParameterInference)
library(ggplot2)
data <- read.csv(system.file("extdata", "data.csv", package = "IICDWorkshopParameterInference"))
set.seed(1)
N_trials <- 30000
parameter_ranges <- data.frame(
    alpha = c(0, 2),
    beta = c(0, 2),
    gamma = c(0, 2),
    delta = c(0, 2),
    row.names = c("min", "max")
)
initial_state <- list(
    rabbit = 700,
    fox = 300
)
residual_sum <- function(parameters, initial_state, data) {
    predicted_populations <- Lotka_Volterra_ODE_solver(
        parameters = parameters,
        times = data$time,
        initial_state = initial_state
    )
    residual_sum <- sum((data$rabbit - predicted_populations$rabbit)^2 + (data$fox - predicted_populations$fox)^2)
    return(residual_sum)
}
trials <- data.frame(
    alpha = runif(N_trials, min = parameter_ranges["min", "alpha"], max = parameter_ranges["max", "alpha"]),
    beta = runif(N_trials, min = parameter_ranges["min", "beta"], max = parameter_ranges["max", "beta"]),
    gamma = runif(N_trials, min = parameter_ranges["min", "gamma"], max = parameter_ranges["max", "gamma"]),
    delta = runif(N_trials, min = parameter_ranges["min", "delta"], max = parameter_ranges["max", "delta"]),
    residual_sum = numeric(N_trials)
)
pb <- txtProgressBar(
    min = 0, max = N_trials,
    style = 3, width = 50, char = "+"
)
for (i in seq_len(N_trials)) {
    setTxtProgressBar(pb, i)
    trials$residual_sum[i] <- residual_sum(as.list(trials[i, ]), initial_state, data)
}
best_parameters <- trials[which.min(trials$residual_sum), , drop = FALSE]
predicted_populations <- Lotka_Volterra_ODE_solver(
    parameters = best_parameters,
    times = seq(0, 14, by = 0.1),
    initial_state = initial_state
)
#-----------------------------------------------------------------------Plot fitted & observed population dynamics
color_rabbit <- "#EFC000"
color_fox <- "#BC3C29"
obs_df <- rbind(
    data.frame(time = data$time, species = "Rabbit", count = data$rabbit),
    data.frame(time = data$time, species = "Fox", count = data$fox)
)
obs_df$species <- factor(obs_df$species, levels = c("Rabbit", "Fox"))
fit_df <- rbind(
    data.frame(time = predicted_populations$time, species = "Rabbit", count = predicted_populations$rabbit),
    data.frame(time = predicted_populations$time, species = "Fox", count = predicted_populations$fox)
)
fit_df$species <- factor(fit_df$species, levels = c("Rabbit", "Fox"))
p_fit <- ggplot() +
    geom_line(
        data = fit_df,
        aes(x = time, y = count, colour = species, group = species),
        linewidth = 1.2
    ) +
    geom_point(
        data = obs_df,
        aes(x = time, y = count, colour = species),
        size = 7,
        stroke = 0.3
    ) +
    scale_colour_manual(
        values = c(Rabbit = color_rabbit, Fox = color_fox),
        name = NULL
    ) +
    labs(x = "Time", y = "Count") +
    theme_bw(base_size = 20) +
    theme(
        plot.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.background = element_rect(
            fill = "white",
            colour = "grey35",
            linewidth = 0.5
        ),
        legend.margin = margin(6, 8, 6, 8)
    )
ggsave(
    filename = "Interactive_session_2_population_dynamics.png",
    plot = p_fit,
    width = 10,
    height = 6,
    dpi = 300,
    bg = "white"
)
#-----------------------------------------------------------------------Plot inferred parameter values
param_names <- c("alpha", "beta", "gamma", "delta")
long_pe <- data.frame(
    parameter = factor(param_names, levels = param_names),
    value = as.numeric(best_parameters[1, param_names])
)
p_pe <- ggplot(long_pe, aes(x = value)) +
    geom_vline(aes(xintercept = value), colour = "#08519C", linewidth = 1.2) +
    geom_point(aes(y = 0), colour = "#08519C", size = 4) +
    geom_text(
        aes(y = 0, label = sprintf("%.3f", value)),
        colour = "#08519C", vjust = -1.5, size = 4.5, fontface = "bold"
    ) +
    facet_wrap(~parameter, scales = "free_y", ncol = 2) +
    coord_cartesian(xlim = c(0, 2)) +
    labs(x = NULL, y = "Point estimate") +
    theme_bw(base_size = 20) +
    theme(
        axis.text.y  = element_blank(),
        panel.grid   = element_blank()
    )
ggsave(
    filename = "Interactive_session_2_parameters.png",
    plot = p_pe,
    width = 10, height = 6, dpi = 300, bg = "white"
)
