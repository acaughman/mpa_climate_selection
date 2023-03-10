library(rfishbase)
library(tidyverse)

data <- read_csv(here::here("01_raw_data", "species.csv")) %>%
  select(species = Species)
validate_names(data$species)

fish <- load_taxa()

species <- species(fish$Species)
stocks <- stocks(fish$Species)
mat <- maturity(fish$Species)
fec <- fecundity(fish$Species)
pop <- popgrowth(fish$Species)

dat <- full_join(species, stocks)
dat <- full_join(dat, mat)
dat <- full_join(dat, fec)
dat <- full_join(dat, pop) %>%
  filter(Species %in% data$species) %>%
  select(
    Species, M,
    AgeMatMin,
    TempMin, TempMax
  ) %>%
  rename(species = Species) %>% 
  mutate(ThermalBreadth = TempMax - TempMin)

dat_full <- full_join(dat, data) %>% 
  group_by(species) %>% 
  summarize(M = mean(M, na.rm = TRUE),
            AgeMatMin = mean(AgeMatMin, na.rm =TRUE),
            ThermalBreadth = mean(ThermalBreadth, na.rm = TRUE)) %>% 
  mutate(SurvivalRate = exp(-M)) %>% 
  select(-M)

write_csv(dat_full, here::here("03_generated_data", "species_sum.csv"))
