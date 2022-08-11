library(rfishbase)
library(tidyverse)

data = read_csv(here::here("data","species.csv")) %>% 
  select(species, homerange, measure)
validate_names(data$species)

fish = load_taxa()

species = species(fish$Species)
stocks = stocks(fish$Species)
mat = maturity(fish$Species)
fec = fecundity(fish$Species)
pop = popgrowth(fish$Species)

dat = full_join(species, stocks)
dat = full_join(dat, mat)
dat = full_join(dat, fec)
dat = full_join(dat, pop) %>% 
  filter(Species %in% data$species) %>% 
  select(Species, M,
         FecundityMax,tm,AgeMatMin,
         TempPreferred,TempPref50,
         TempMin,TempMax,FBname) %>% 
  rename(species = Species)

dat_full = full_join(dat, data)

write_csv(dat_full, here::here("species_sum.csv"))
