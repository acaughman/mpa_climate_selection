library(ggplot2)
library(tidyverse)


# Model Parameter set up --------------------------------------------------

SST <- seq(0, 50) # sea surface temperatures to model

opt.temp <- 25 # optimal temperature ranges

temp.range <- 4 # thermal performance range, w


# Mortality Function ------------------------------------------------------

calc_mortality <- function(SST, opt.temp, temp.range) {
  nat.m <- 1 - exp((-(SST - opt.temp)^2) / (temp.range^2)) # temperature based mortality function from Walsworth et al.
  return(1 - nat.m)
}


# Simulate Mortality Rates Function -----------------------------------------------

simulate_mortality <- function(SST, opt.temp, temp.range) {
  results <- matrix(NaN, nrow = length(SST), ncol = 1) # empty matrix to hold results


  for (i in 1:length(SST)) {
    for (j in 1:length(opt.temp)) {
      mort <- calc_mortality(SST[i], opt.temp, temp.range)
      results[i][j] <- mort
    }
  }

  results_plot <- as.data.frame(results) %>%
    rename(survival = V1) %>%
    mutate(survival.prop = case_when(
      survival > .7 ~ .7,
      survival < .7 ~ survival
    ))

  results_plot$SST <- SST
  results_plot$opt.temp <- opt.temp
  results_plot$range <- temp.range

  return(results_plot)
}



# Simulate multiple optimal temperatures ----------------------------------

simulation_results <- data.frame()

for (k in 1:length(opt.temp)) {
  sub_result <- simulate_mortality(SST, opt.temp[k], temp.range)
  simulation_results <- bind_rows(sub_result, simulation_results)
}

simulation_results <- simulation_results %>%
  mutate(opt.temp = as.factor(opt.temp))

ggplot(simulation_results, aes(SST, survival.prop, color = opt.temp)) +
  geom_line() +
  theme_bw()


# Simulate multiple range widths ------------------------------------------

opt.temp <- 22 # optimal temperature ranges

temp.range <- c(2, 4, 6, 8, 10) # thermal performance range, w

simulation_results <- data.frame()

for (k in 1:length(temp.range)) {
  sub_result <- simulate_mortaility(SST, opt.temp, temp.range[k])
  simulation_results <- bind_rows(sub_result, simulation_results)
}

simulation_results <- simulation_results %>%
  mutate(range = as.factor(range))

ggplot(simulation_results, aes(SST, survival.prop, color = range)) +
  geom_path() +
  theme_bw()
