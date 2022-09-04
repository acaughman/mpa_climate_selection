library(rfishbase)
library(tidyverse)

data = read_csv(here::here("01_raw_data","species.csv")) %>% 
  select(species)
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

write_csv(dat_full, here::here("03_generated_data","species_sum.csv"))
