# Set Up ------------------------------------------------------------------

library(here)
library(tidyverse)
library(ggplot2)
library(patchwork)

# Read in Data ------------------------------------------------------------

MPA2x2F8 = read_csv(here("intermediate_data", "2x2MPA.8F_freq.csv")) %>% 
  mutate(MPA = "2x2") %>% 
  mutate(Fish = .8) %>% 
  filter(generation %in% c(50, 100, 125, 150, 175, 200)) %>% 
  mutate(generation = as.numeric(generation))
MPA2x2F5 = read_csv(here("intermediate_data", "2x2MPA.5F_freq.csv")) %>% 
  mutate(MPA = "2x2") %>% 
  mutate(Fish = .5) %>% 
  filter(generation %in% c(50, 100, 125, 150, 175, 200)) %>% 
  mutate(generation = as.numeric(generation))
MPA2x2F2 = read_csv(here("intermediate_data", "2x2MPA.2F_freq.csv")) %>% 
  mutate(MPA = "2x2") %>% 
  mutate(Fish = .5) %>% 
  filter(generation %in% c(50, 100, 125, 150, 175, 200)) %>% 
  mutate(generation = as.numeric(generation))

MPA1x1F8 = read_csv(here("intermediate_data", "1x1MPA.8F_freq.csv")) %>% 
  mutate(MPA = "1x1") %>% 
  mutate(Fish = .8) %>% 
  filter(generation %in% c(50, 100, 125, 150, 175, 200)) %>% 
  mutate(generation = as.numeric(generation))
MPA1x1F5 = read_csv(here("intermediate_data", "1x1MPA.5F_freq.csv")) %>% 
  mutate(MPA = "1x1") %>% 
  mutate(Fish = .8) %>% 
  filter(generation %in% c(50, 100, 125, 150, 175, 200)) %>% 
  mutate(generation = as.numeric(generation))
MPA1x1F2 = read_csv(here("intermediate_data", "1x1MPA.2F_freq.csv")) %>% 
  mutate(MPA = "1x1") %>% 
  mutate(Fish = .8) %>% 
  filter(generation %in% c(50, 100, 125, 150, 175, 200)) %>% 
  mutate(generation = as.numeric(generation))

MPA4x4F8 = read_csv(here("intermediate_data", "MPA4x4F8_freq.csv")) %>% 
  mutate(MPA = "4x4") %>% 
  mutate(Fish = .8) %>% 
  filter(generation %in% c(50, 100, 125, 150, 175, 200)) %>% 
  mutate(generation = as.numeric(generation))

output_sum = do.call("rbind", list(MPA2x2F8,MPA2x2F5,MPA2x2F2,MPA1x1F8,MPA1x1F2,MPA1x1F5))

no_MPA = read_csv(here("intermediate_data", "no_MPA_freq.csv")) %>%
  mutate(generation = as.numeric(generation)) %>% 
  filter(generation %in% c(50, 100, 125, 150, 175, 200)) %>% 
  mutate(generation = as.numeric(generation))

no_F = read_csv(here("intermediate_data", "no_F_freq.csv")) %>%
  mutate(generation = as.numeric(generation)) %>% 
  filter(generation %in% c(50, 100, 125, 150, 175, 200)) %>% 
  mutate(generation = as.numeric(generation))


# Plots 2x2 MPA different Fishing Pressure -------------------------------------------------------------------

pop1 = ggplot(MPA2x2F8, aes(lon, lat, color = freq, fill = freq)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Genotype Frequency", color = "Genotype Frequency") 

freq1 = ggplot(MPA2x2F8, aes(lon, lat, color = geno_pop_sum, fill = geno_pop_sum)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Population Size", color = "Population Size")

pop2 = ggplot(MPA2x2F5, aes(lon, lat, color = freq, fill = freq)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Genotype Frequency", color = "Genotype Frequency") 

freq2 = ggplot(MPA2x2F5, aes(lon, lat, color = geno_pop_sum, fill = geno_pop_sum)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Population Size", color = "Population Size")

pop3 = ggplot(MPA2x2F2, aes(lon, lat, color = freq, fill = freq)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Genotype Frequency", color = "Genotype Frequency") 

freq3 = ggplot(MPA2x2F2, aes(lon, lat, color = geno_pop_sum, fill = geno_pop_sum)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Population Size", color = "Population Size")


p1 = (pop1 / freq1) 
p2 = (pop2 / freq2)
p3 = (pop3 / freq3)

ggsave(p1, file = "MPA2x2F8.png",path = here("outputs"))
ggsave(p3, file = "MPA2x2F2.png",path = here("outputs"))


# Plot 1x1 MPA different Fishing Pressure ---------------------------------

pop1 = ggplot(MPA1x1F8, aes(lon, lat, color = freq, fill = freq)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Genotype Frequency", color = "Genotype Frequency") 

freq1 = ggplot(MPA1x1F8, aes(lon, lat, color = geno_pop_sum, fill = geno_pop_sum)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Population Size", color = "Population Size")

pop2 = ggplot(MPA1x1F5, aes(lon, lat, color = freq, fill = freq)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Genotype Frequency", color = "Genotype Frequency") 

freq2 = ggplot(MPA1x1F5, aes(lon, lat, color = geno_pop_sum, fill = geno_pop_sum)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Population Size", color = "Population Size")

pop3 = ggplot(MPA1x1F2, aes(lon, lat, color = freq, fill = freq)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Genotype Frequency", color = "Genotype Frequency") 

freq3 = ggplot(MPA1x1F2, aes(lon, lat, color = geno_pop_sum, fill = geno_pop_sum)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Population Size", color = "Population Size")


p1 = (pop1 / freq1) 
p2 = (pop2 / freq2)
p3 = (pop3 / freq3)

plot2 = p1 + p2 + p3

#ggsave(plot, file = "test_fig.png",path = here("outputs"))


# Fishing Pressure .8 -----------------------------------------------------

pop1 = ggplot(MPA1x1F8, aes(lon, lat, color = freq, fill = freq)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Genotype Frequency", color = "Genotype Frequency") 

freq1 = ggplot(MPA1x1F8, aes(lon, lat, color = geno_pop_sum, fill = geno_pop_sum)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Population Size", color = "Population Size")

pop2 = ggplot(MPA2x2F8, aes(lon, lat, color = freq, fill = freq)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Genotype Frequency", color = "Genotype Frequency") 

freq2 = ggplot(MPA2x2F8, aes(lon, lat, color = geno_pop_sum, fill = geno_pop_sum)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Population Size", color = "Population Size")

pop3 = ggplot(MPA2x2F8, aes(lon, lat, color = freq, fill = freq)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Genotype Frequency", color = "Genotype Frequency") 

freq3 = ggplot(MPA4x4F8, aes(lon, lat, color = geno_pop_sum, fill = geno_pop_sum)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Population Size", color = "Population Size")

p1 = pop1 / freq1 
p2 = pop2 / freq2
p3 = pop3 / freq3

plot3 = p1 + p2

ggsave(p1, file = "MPA1x1F8.png",path = here("outputs"))


# Fishing Pressure .5 -----------------------------------------------------
pop1 = ggplot(MPA1x1F5, aes(lon, lat, color = freq, fill = freq)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Genotype Frequency", color = "Genotype Frequency") 

freq1 = ggplot(MPA1x1F5, aes(lon, lat, color = geno_pop_sum, fill = geno_pop_sum)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Population Size", color = "Population Size")

pop2 = ggplot(MPA2x2F5, aes(lon, lat, color = freq, fill = freq)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Genotype Frequency", color = "Genotype Frequency") 

freq2 = ggplot(MPA2x2F5, aes(lon, lat, color = geno_pop_sum, fill = geno_pop_sum)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Population Size", color = "Population Size")

p1 = (pop1 / freq1) 
p2 = (pop2 / freq2)

plot4 = p1 + p2

#ggsave(plot, file = "test_fig.png",path = here("outputs"))

# Fishing Pressure .2 -----------------------------------------------------

pop1 = ggplot(MPA1x1F2, aes(lon, lat, color = freq, fill = freq)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Genotype Frequency", color = "Genotype Frequency") 

freq1 = ggplot(MPA1x1F2, aes(lon, lat, color = geno_pop_sum, fill = geno_pop_sum)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Population Size", color = "Population Size")

pop2 = ggplot(MPA2x2F2, aes(lon, lat, color = freq, fill = freq)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Genotype Frequency", color = "Genotype Frequency") 

freq2 = ggplot(MPA2x2F2, aes(lon, lat, color = geno_pop_sum, fill = geno_pop_sum)) +
  geom_tile() + facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Population Size", color = "Population Size")

p1 = (pop1 / freq1) 
p2 = (pop2 / freq2)

plot5 = p1 + p2

#ggsave(plot, file = "test_fig.png",path = here("outputs"))


# No MPA ------------------------------------------------------------------

# pop1 = ggplot(no_MPA, aes(lon, lat, color = freq, fill = freq)) +
#   geom_tile() + facet_grid(genotype~generation) + 
#   labs(x = "Longitude", y = "Latitude", fill = "Genotype Frequency", color = "Genotype Frequency") 
# 
# freq1 = ggplot(no_MPA, aes(lon, lat, color = geno_pop_sum, fill = geno_pop_sum)) +
#   geom_tile() + facet_grid(genotype~generation) + 
#   labs(x = "Longitude", y = "Latitude", fill = "Population Size", color = "Population Size")
# 
# p1 = (pop1 / freq1) 

# ggplot(no_MPA, aes(freq)) + 
#   geom_histogram(bins = 15) +
#   facet_grid(~genotype)
# 
# ggplot(no_MPA, aes(geno_pop_sum)) + 
#   geom_histogram(bins = 15) +
#   facet_grid(~genotype)


# No Fishing --------------------------------------------------------------

pop1 = ggplot(no_F, aes(lon, lat, color = freq, fill = freq)) +
  geom_tile() + facet_grid(genotype~generation) +
  labs(x = "Longitude", y = "Latitude", fill = "Genotype Frequency", color = "Genotype Frequency")

freq1 = ggplot(no_F, aes(lon, lat, color = geno_pop_sum, fill = geno_pop_sum)) +
  geom_tile() + facet_grid(genotype~generation) +
  labs(x = "Longitude", y = "Latitude", fill = "Population Size", color = "Population Size")

p1 = pop1 / freq1

