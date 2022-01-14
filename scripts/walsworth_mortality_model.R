library(ggplot2)
library(tidyverse)


# Model Parameter set up --------------------------------------------------

SST = seq(0, 35) #sea surface temperatures to model

opt.temp = 25 # optimal temperature ranges

temp.range = 8 # thermal performance range, w

# Morality Rate Function --------------------------------------------------

calc_mortality <- function(SST, opt.temp, temp.range) {
  
  return(nat.m)
}

# Simulate Mortality Rates -----------------------------------------------

results = matrix(NaN, nrow = length(SST), ncol = length(opt.temp))

for (i in 1:length(SST)) {
  for (j in 1:length(opt.temp)) {
    temp.range = 8
    mort = calc_mortality(SST[i], opt.temp[j], temp.range)
    results[i][j] = mort
  }
}

results = as.data.frame(results) %>% 
  rename(mortality = V1) %>% 
  mutate(survival.prop = exp(-mortality))

results$SST = SST

ggplot(results, aes(SST, survival.prop)) +
  geom_path() +
  theme_bw()

