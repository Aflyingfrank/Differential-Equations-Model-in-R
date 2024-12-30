#############################################################################################################################################
############################################## Differential equations model with common usages ##############################################
#############################################################################################################################################
rm(list=ls(all=TRUE)) # clear the environment

################################################## Library all the packages to be used ##################################################
library(deSolve) # for solving ODEs
library(ggplot2)
library(tidyverse)

########################################## Malthusian Population Model & Logistic Population Model ##########################################
########## Raw data input ##########
# Initial raw data
population_data <- data.frame(time = seq(0, 210, by = 10),
                              observed = c(3.9, 5.3, 7.2, 9.6, 12.9, 17.1, 23.2, 31.4, 38.6, 50.2, 62.9, 76.0, 92.0, 106.5, 123.2, 131.7, 150.7, 179.3, 204.0, 226.5, 251.4, 281.4))

# Initial population size
initial_population <- population_data$observed[1]

########## Malthusian Growth Model ODE ##########
# Define the Malthusian model as a function for the ODE solver
malthusian_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dP <- r * P  # dP/dt = r * P
    list(c(dP))
  })
}

# Objective function for Malthusian model to minimize SSR
malthusian_objective_ode <- function(params, data) {
  r <- params[1]
  initial_state <- c(P = data$observed[1])  # Initial population size
  
  # Solve the Malthusian ODE
  times <- data$time
  solution <- ode(y = initial_state, times = times, func = malthusian_model, parms = list(r = r))
  
  # Calculate the residual sum of squares (RSS)
  predicted <- solution[, "P"]
  observed <- data$observed
  ssr <- sum((observed - predicted)^2)
  
  return(ssr)
}

# Estimate the growth rate r for the Malthusian model using optim
malthusian_fit_ode <- optim(par = c(r = 0.02), fn = malthusian_objective_ode, data = population_data)
r_malthusian <- malthusian_fit_ode$par["r"]

# Solve the Malthusian ODE using the estimated parameter
malthusian_solution <- ode(y = c(P = initial_population), times = population_data$time, func = malthusian_model, parms = list(r = r_malthusian))

########## Logistic Growth Model ODE ##########
# Define the Logistic model as a function for the ODE solver
logistic_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dP <- r * P * (1 - P / K)  # dP/dt = r * P * (1 - P / K)
    list(c(dP))
  })
}

# Objective function for Logistic model to minimize SSR
logistic_objective_ode <- function(params, data) {
  r <- params[1]
  K <- params[2]
  initial_state <- c(P = data$observed[1])  # Initial population size
  
  # Solve the Logistic ODE
  times <- data$time
  solution <- ode(y = initial_state, times = times, func = logistic_model, parms = list(r = r, K = K))
  
  # Calculate the residual sum of squares (RSS)
  predicted <- solution[, "P"]
  observed <- data$observed
  ssr <- sum((observed - predicted)^2)
  
  return(ssr)
}

# Initial guess for logistic model parameters
initial_guess_logistic <- c(r = 0.02, K = max(population_data$observed) * 1.2)

# Estimate the growth rate r and carrying capacity K for the Logistic model using optim
logistic_fit_ode <- optim(par = initial_guess_logistic, fn = logistic_objective_ode, data = population_data)
r_logistic <- logistic_fit_ode$par["r"]
K_logistic <- logistic_fit_ode$par["K"]

# Solve the Logistic ODE using the estimated parameters
logistic_solution <- ode(y = c(P = initial_population), times = population_data$time, func = logistic_model, parms = list(r = r_logistic, K = K_logistic))

########## Visualize the Results ##########
# Create a data frame for plotting
malthusian_predicted <- data.frame(time = malthusian_solution[, "time"], predicted = malthusian_solution[, "P"])
logistic_predicted <- data.frame(time = logistic_solution[, "time"], predicted = logistic_solution[, "P"])

rm(list=ls()[!ls() %in% c("malthusian_predicted", "logistic_predicted", "r_malthusian", "r_logistic", "K_logistic", "population_data", "malthusian_model", "logistic_model")])

# Plot the observed and model predictions
population.model.plot <- ggplot() +
  geom_point(data = population_data, aes(x = time, y = observed), color = "blue", size = 2) +
  geom_line(data = malthusian_predicted, aes(x = time, y = predicted), color = "coral", linewidth = 1, linetype = "solid") +
  geom_line(data = logistic_predicted, aes(x = time, y = predicted), color = "darkgreen", linewidth = 1, linetype = "solid") +
  labs(title = "",
       x = "Time (Years since 1790)", y = "Population Size (Millions)") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, max(population_data$observed) * 1.2)) +
  theme(legend.position = "bottom") +
  annotate("text", x = 50, y = max(population_data$observed) * 0.8, 
           label = paste("Malthusian Growth Rate (r) =", round(r_malthusian, 4)), color = "coral") +
  annotate("text", x = 50, y = max(population_data$observed) * 0.7, 
           label = paste("Logistic Growth Rate (r) =", round(r_logistic, 4)), color = "darkgreen") +
  annotate("text", x = 50, y = max(population_data$observed) * 0.6, 
           label = paste("Logistic Carrying Capacity (K) =", round(K_logistic, 4)), color = "darkgreen") +
  theme(
    plot.title = element_text(size = 16),      # Title size
    axis.title.x = element_text(size = 14),    # X-axis title size
    axis.title.y = element_text(size = 14),    # Y-axis title size
    axis.text.x = element_text(size = 12),     # X-axis labels size
    axis.text.y = element_text(size = 12)      # Y-axis labels size
  )

ggsave("Malthusian & Logistic Population Model Plot.png", plot = population.model.plot, width = 10, height = 9, dpi = 300, bg='white')

########## Numerical predictions ##########
predict.time <- seq(0, 240, by = 10) # predict 30 more years
initial.population <- c(P = population_data$observed[1])

malthusian.prediction <- ode(y = initial.population, times = predict.time, func = malthusian_model, parms = r_malthusian)
malthusian.prediction.df <- as.data.frame(malthusian.prediction)

logistic.prediction <- ode(y = initial.population, times = predict.time, func = logistic_model, parms = c(r_logistic, K_logistic))
logistic.prediction.df <- as.data.frame(logistic.prediction)

rm(list=ls()[!ls() %in% c("malthusian.prediction.df", "logistic.prediction.df", "population_data")])

prediction <- population_data %>%
  full_join(malthusian.prediction.df, by = "time") %>%
  full_join(logistic.prediction.df, by = "time") %>%
  rename(observed = observed,   
         malthusian = P.x,  
         logistic = P.y)   

########## Visualize the predictions ##########
population.prediction.plot <- ggplot() +
  geom_point(data = population_data, aes(x = time, y = observed), color = "blue", size = 2) +
  geom_line(data = prediction, aes(x = time, y = malthusian), color = "coral", linewidth = 1, linetype = "solid") +
  geom_line(data = prediction, aes(x = time, y = logistic), color = "darkgreen", linewidth = 1, linetype = "solid") +
  labs(title = "",
       x = "Years since 1970", y = "Population Size (Millions)") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, max(population_data$observed) * 1.2)) +
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = seq(0, 240, by = 20)) +
  scale_y_continuous(breaks = seq(0, 650, by = 100)) +
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  )
population.prediction.plot <- population.prediction.plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
population.prediction.plot

# Real time scale (1790 until 2030)
prediction.real <- prediction %>% mutate(time = 1790 + prediction$time)

population.prediction.plot <- ggplot() +
  geom_point(data = prediction.real, aes(x = time, y = observed), color = "blue", size = 2) +
  geom_line(data = prediction.real, aes(x = time, y = malthusian), color = "coral", linewidth = 1.5, linetype = "solid") +
  geom_line(data = prediction.real, aes(x = time, y = logistic), color = "darkgreen", linewidth = 1.5, linetype = "solid") +
  labs(title = "",
       x = "Years", y = "Population Size (Millions)") +
  theme_classic() +
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = seq(1790, 2030, by = 10), labels = ifelse(seq(1790, 2030, by = 10) %% 20 == 10, seq(1790, 2030, by = 10), ""))+
  scale_y_continuous(breaks = seq(0, 650, by = 50), labels = ifelse(seq(0, 650, by = 50) %% 100 == 0, seq(0, 650, by = 50), "")) +
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  )
population.prediction.plot <- population.prediction.plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))

xs <- 1810
ys <- 400
population.prediction.plot <- population.prediction.plot + 
  annotate("text", x = xs+6.8, y = ys+110, label = "Observed data", color = "black", size = 6, hjust = 0) +
  annotate("point", x = xs+3, y = ys+110, color = "blue", size = 2) +
  annotate("text", x = xs+6.8, y = ys+70, label = "Malthusian model", color = "black", size = 6, hjust = 0) +
  annotate("segment", x = xs, xend = xs+6, y = ys+70, yend = ys+70, color = "coral", linewidth = 1.5, linetype = "solid") +
  annotate("text", x = xs+6.8, y = ys+30, label = "Logistic model", color = "black", size = 6, hjust = 0) +
  annotate("segment", x = xs, xend = xs+6, y = ys+30, yend = ys+30, color = "darkgreen", linewidth = 1.5, linetype = "solid")
population.prediction.plot

#################################################### Lotka-Volterra Predator Prey Model ####################################################
rm(list=ls(all=TRUE)) # clear the environment

########## Define the Lotka-Volterra Predator Prey Model (x1 = prey, x2 = predator) ##########
volterra_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dx1 <- x1 * (r1 - lamda1 * x2 - e)  # dx1/dt = x1 * (r1 - lamda * x1)
    dx2 <- x2 * (-r2 + lamda2 * x1 - e)  # dx1/dt = x1 * (r1 - lamda * x1)
    list(c(dx1, dx2))
  })
}

########## Set the initial conditions, time and parameters ##########
predict.time <- seq(0, 20, by = 1)
initial.condition <- c(x1 = 25, x2 = 2)
parameters.bw <- c(r1 = 1,
                   lamda1 = 0.1,
                   r2 = 0.5,
                   lamda2 = 0.02,
                   e = 0.3)
parameters.aw <- c(r1 = 1,
                   lamda1 = 0.1,
                   r2 = 0.5,
                   lamda2 = 0.02,
                   e = 0.1)

########## Solve for the ODEs ##########
volterra.bw <- ode(y = initial.condition, times = predict.time, func = volterra_model, parms = parameters.bw)
volterra.bw.df <- as.data.frame(volterra.bw)

volterra.aw <- ode(y = initial.condition, times = predict.time, func = volterra_model, parms = parameters.aw)
volterra.aw.df <- as.data.frame(volterra.aw)

rm(list=ls()[!ls() %in% c("volterra.aw.df", "volterra.bw.df")])

########## Visualize the results ##########
volterra.plot <- ggplot() +
  geom_line(data = volterra.bw.df, aes(x = time, y = x1), color = "blue", linewidth = 0.8, linetype = "solid") +
  geom_line(data = volterra.bw.df, aes(x = time, y = x2), color = "red", linewidth = 0.8, linetype = "solid") +
  geom_line(data = volterra.aw.df, aes(x = time, y = x1), color = "blue", linewidth = 0.8, linetype = "dashed") +
  geom_line(data = volterra.aw.df, aes(x = time, y = x2), color = "red", linewidth = 0.8, linetype = "dashed") +
  labs(title = "",
       x = "Time", y = "Population Size") +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, 20, by = 1), labels = ifelse(seq(0, 20, by = 1) %% 2 == 0, seq(0, 20, by = 1), "")) +
  scale_y_continuous(breaks = seq(0, 100, by = 10)) +
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )
xs <- 15.5
ys <- 83
volterra.plot <- volterra.plot + 
  annotate("text", x = xs+1, y = ys+9, label = "Prey before war", color = "black", size = 4, hjust = 0) +
  annotate("segment", x = xs, xend = xs+0.9, y = ys+9, yend = ys+9, color = "blue", linewidth = 1, linetype = "solid") +
  annotate("text", x = xs+1, y = ys+6, label = "Predator before war", color = "black", size = 4, hjust = 0) +
  annotate("segment", x = xs, xend = xs+0.9, y = ys+6, yend = ys+6, color = "red", linewidth = 1, linetype = "solid") +
  annotate("text", x = xs+1, y = ys+3, label = "Prey after war", color = "black", size = 4, hjust = 0) +
  annotate("segment", x = xs, xend = xs+0.9, y = ys+3, yend = ys+3, color = "blue", linewidth = 1, linetype = "dashed") +
  annotate("text", x = xs+1, y = ys, label = "Predator after war", color = "black", size = 4, hjust = 0) +
  annotate("segment", x = xs, xend = xs+0.9, y = ys, yend = ys, color = "red", linewidth = 1, linetype = "dashed")
volterra.plot

################################################ Lotka-Volterra + Species Competition Model ################################################
rm(list=ls(all=TRUE)) # clear the environment

########## Define the Lotka-Volterra + Species Competition Model ##########
model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dx1 <- x1 * (r1 - lamda1 * x2 - e) - r1 * x1 * (x1 / N1 + sigma1 * x3 / N3)
    dx2 <- x2 * (-r2 + lamda2 * x1 - e)
    dx3 <- r3 * x3 * (1 - x3 / N3 - sigma3 * x1 / N3)
    list(c(dx1, dx2, dx3))
  })
}

########## Set the initial conditions, time and parameters ##########
predict.time <- seq(0, 20, by = 1)
initial.condition <- c(x1 = 20, x2 = 20, x3 = 20)
parameters <- c(r1 = 1,
                lamda1 = 0.1,
                sigma1 = 0.5,
                N1 = 500,
                r2 = 0.5,
                lamda2 = 0.02,
                r3 = 0.2,
                sigma3 = 2,
                N3 = 300,
                e = 0.3)

########## Solve for the ODEs ##########
volterra.sc <- ode(y = initial.condition, times = predict.time, func = model, parms = parameters)
volterra.sc.df <- as.data.frame(volterra.sc)

rm(list=ls()[!ls() %in% c("volterra.sc.df")])

########## Visualize the results ##########
prediction.plot <- ggplot() +
  geom_line(data = volterra.sc.df, aes(x = time, y = x1), color = "blue", linewidth = 0.8, linetype = "solid") +
  geom_line(data = volterra.sc.df, aes(x = time, y = x2), color = "red", linewidth = 0.8, linetype = "solid") +
  geom_line(data = volterra.sc.df, aes(x = time, y = x3), color = "green", linewidth = 0.8, linetype = "solid") +
  labs(title = "",
       x = "Time", y = "Population Size") +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, 20, by = 1), labels = ifelse(seq(0, 20, by = 1) %% 2 == 0, seq(0, 20, by = 1), "")) +
  scale_y_continuous(breaks = seq(0, 160, by = 10), labels = ifelse(seq(0, 160, by = 10) %% 20 == 0, seq(0, 160, by = 10), "")) +
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )
xs <- 1
ys <- 120
prediction.plot <- prediction.plot + 
  annotate("text", x = xs+1, y = ys+16, label = "Species A", color = "black", size = 4, hjust = 0) +
  annotate("segment", x = xs, xend = xs+0.9, y = ys+16, yend = ys+16, color = "blue", linewidth = 1, linetype = "solid") +
  annotate("text", x = xs+1, y = ys+8, label = "Species B", color = "black", size = 4, hjust = 0) +
  annotate("segment", x = xs, xend = xs+0.9, y = ys+8, yend = ys+8, color = "red", linewidth = 1, linetype = "solid") +
  annotate("text", x = xs+1, y = ys, label = "Species C", color = "black", size = 4, hjust = 0) +
  annotate("segment", x = xs, xend = xs+0.9, y = ys, yend = ys, color = "green", linewidth = 1, linetype = "solid")
prediction.plot














