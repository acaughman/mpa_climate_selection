library(here)
library(tidyverse)
library(ggplot2)
library(patchwork)


# Read in Data ------------------------------------------------------------

MPA2x2F8 = read_csv(here("intermediate_data", "2x2MPA.8F_freq.csv")) %>% 
  mutate(MPA = "2x2") %>% 
  mutate(Fish = .8)
MPA2x2F5 = read_csv(here("intermediate_data", "2x2MPA.5F_freq.csv")) %>% 
  mutate(MPA = "2x2") %>% 
  mutate(Fish = .5)
MPA2x2F2 = read_csv(here("intermediate_data", "2x2MPA.2F_freq.csv")) %>% 
  mutate(MPA = "2x2") %>% 
  mutate(Fish = .5)

MPA1x1F8 = read_csv(here("intermediate_data", "1x1MPA.8F_freq.csv")) %>% 
  mutate(MPA = "1x1") %>% 
  mutate(Fish = .8)
MPA1x1F5 = read_csv(here("intermediate_data", "1x1MPA.5F_freq.csv")) %>% 
  mutate(MPA = "1x1") %>% 
  mutate(Fish = .8)
MPA1x1F2 = read_csv(here("intermediate_data", "1x1MPA.2F_freq.csv")) %>% 
  mutate(MPA = "1x1") %>% 
  mutate(Fish = .8)

output_sum = do.call("rbind", list(MPA2x2F8,MPA2x2F5,MPA2x2F2,MPA1x1F8,MPA1x1F2,MPA1x1F5))

# Plots -------------------------------------------------------------------

years = c(50, 100, 125, 150, 175, 200)

plot_sum = output_sum %>% 
  filter(generation %in% c(50, 100, 125, 150, 175, 200)) %>% 
  mutate(generation = as.numeric(generation))

p1 = ggplot(plot_sum, aes(lon, lat, color = freq, fill = freq)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Genotype Frequency", color = "Genotype Frequency") 

p2 = ggplot(plot_sum, aes(lon, lat, color = geno_pop_sum, fill = geno_pop_sum)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Population Size", color = "Population Size")

p2 / p1

plot = p2/p1

#ggsave(plot, file = "test_fig.png",path = here("outputs"))
