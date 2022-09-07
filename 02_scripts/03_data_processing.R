library(tidyverse)
library(patchwork)

# Allie Explore -----------------------------------------------------------

reps = 1
NS.patches <- 100 # the number of patches on the north-south axis
EW.patches <- 20
NUM.age.classes <- 3 #babies, juvenile, adult
NUM.sexes <- 2 #female male
NUM.genotypes <- 3 #AA,Aa,aa
NUM.gens.pre.fishing <- 15 # The number of generations before any fishery
NUM.gens.pre.reserve <- 10 # The number of generations of fishing before reserves are installed
NUM.gens.post.reserve <- 150 # The number of generations with the reserve installed
gens = NUM.gens.pre.fishing+NUM.gens.pre.reserve+NUM.gens.post.reserve

load(file = here::here("sensitivity_analysis", "time", "enso8.rda"))
load(file = here::here("03_generated_data", "climate_layer" ,"enso.rda"))

# Output results into a dataframe
output_df = data.frame() #create dataframe to hold results

for(a in 1:reps) {
  world_sub <- array(0, c(NS.patches, EW.patches))
  for(b in 1:gens) {
    for(c in 1:NUM.genotypes) {
      for(d in 1:NUM.sexes) {
        for(e in 1:NUM.age.classes) {
          world_sub = output.array[,,e,d,c,b,a] %>% 
            as.data.frame()
          world_sub$rep = paste0(a)
          world_sub$generation = paste0(b)
          world_sub$genotype = paste0(c)
          world_sub$sex = paste0(d)
          world_sub$age = paste0(e)
          world_sub$lat = c(1:NS.patches)
          world_sub$max_temp = SST.patches[1,1,b]
          world_sub$min_temp = SST.patches[100,1,b]
          output_df = bind_rows(output_df, world_sub)
        }
      }
    }
  }
}

# Wrangle dataframe into plottable format
output_df = output_df %>% 
  pivot_longer(V1:V20,
               names_to = "lon",
               values_to = "pop") %>% 
  mutate(lon = case_when(
    lon == "V1" ~ 1,
    lon == "V2" ~ 2,
    lon == "V3" ~ 3,
    lon == "V4" ~ 4,
    lon == "V5" ~ 5,
    lon == "V6" ~ 6,
    lon == "V7" ~ 7,
    lon == "V8" ~ 8,
    lon == "V9" ~ 9,
    lon == "V10" ~ 10,
    lon == "V11" ~ 11,
    lon == "V12" ~ 12,
    lon == "V13" ~ 13,
    lon == "V14" ~ 14,
    lon == "V15" ~ 15,
    lon == "V16" ~ 16,
    lon == "V17" ~ 17,
    lon == "V18" ~ 18,
    lon == "V19" ~ 19,
    lon == "V20" ~ 20
  )) %>% 
  mutate(genotype = case_when(
    genotype == 1 ~ "AA",
    genotype == 2 ~ "Aa",
    genotype == 3 ~ "aa"
  )) %>% 
  mutate(genotype = as.factor(genotype)) %>% 
  mutate(sex = case_when(
    sex == 1 ~ "female",
    sex == 2 ~ "male"
  )) %>% 
  mutate(sex = as.factor(sex)) %>%
  mutate(age = case_when(
    age == 1 ~ "baby",
    age == 2 ~ "juvenile",
    age == 3 ~ "adult"
  )) %>% 
  mutate(age = as.factor(age)) %>%
  mutate(lat = as.numeric(lat)) %>% 
  mutate(lon = as.numeric(lon)) %>% 
  mutate(rep = as.numeric(rep))

write_csv(output_df, here::here("sensitivity_analysis","time", "enso8.csv"))