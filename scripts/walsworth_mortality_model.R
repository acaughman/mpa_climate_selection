library(ggplot2)
library(tidyverse)


# Model Parameter set up --------------------------------------------------

SST = seq(0, 35) #sea surface temperatures to model

opt.temp = seq(14,30, by = 2) # optimal temperature ranges

temp.range = 10 # thermal performance range, w


# Mortality Function ------------------------------------------------------

calc_mortality <- function(SST, opt.temp, temp.range) {
  nat.m = 1 - exp((-(SST - opt.temp)^2)/(temp.range^2)) # temperature based mortality function from Walsworth et al.
  return(nat.m)
}


# Simulate Mortality Rates Function -----------------------------------------------

simulate_mortaility <- function(SST, opt.temp, temp.range) {
  results = matrix(NaN, nrow = length(SST), ncol = 1) #empty matrix to hold results
  
  
  for (i in 1:length(SST)) {
    for (j in 1:length(opt.temp)) {
      
      mort = calc_mortality(SST[i], opt.temp, temp.range)
      results[i][j] = 1 - exp(-mort)
    }
  }
  
  results_plot = as.data.frame(results) %>% 
    rename(mortality = V1) %>% 
    mutate(survival.prop = case_when(
      exp(-mortality) > exp(-.37) ~ exp(-.37),
      exp(-mortality) < exp(-.37) ~ exp(-mortality)
    ))
  
  results_plot$SST = SST
  results_plot$opt.temp = opt.temp
  results_plot$range = temp.range
  
  return(results_plot)
}



# Simulate multiple optimal temperatures ----------------------------------

simulation_results = data.frame()

for (k in 1:length(opt.temp)) {
  sub_result = simulate_mortaility(SST, opt.temp[k], temp.range)
  simulation_results = bind_rows(sub_result, simulation_results)
}

simulation_results = simulation_results %>% 
  mutate(opt.temp = as.factor(opt.temp))

ggplot(simulation_results, aes(SST, survival.prop, color = opt.temp)) +
  geom_line() +
  theme_bw()


# Simulate multiple range widths ------------------------------------------

opt.temp = 22 # optimal temperature ranges

temp.range = c(2,4,6,8,10) # thermal performance range, w

simulation_results = data.frame()

for (k in 1:length(temp.range)) {
  sub_result = simulate_mortaility(SST, opt.temp, temp.range[k])
  simulation_results = bind_rows(sub_result, simulation_results)
}

simulation_results = simulation_results %>% 
  mutate(range = as.factor(range))

ggplot(simulation_results, aes(SST, survival.prop, color = range)) +
  geom_path() +
  theme_bw() 

