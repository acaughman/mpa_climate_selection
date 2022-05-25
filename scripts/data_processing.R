library(tidyverse)
library(patchwork)

# Allie Explore -----------------------------------------------------------

reps = 10
NS.patches <- 100 # the number of patches on the north-south axis
EW.patches <- 20
NUM.age.classes <- 3 #babies, juvenile, adult
NUM.sexes <- 2 #female male
NUM.genotypes <- 3 #AA,Aa,aa
NUM.gens.pre.fishing <- 25 # The number of generations before any fishery
NUM.gens.pre.reserve <- 25 # The number of generations of fishing before reserves are installed
NUM.gens.post.reserve <- 100 # The number of generations with the reserve installed
gens = NUM.gens.pre.fishing+NUM.gens.pre.reserve+NUM.gens.post.reserve

load(file = here::here("data", "3x3null8F.rda"))

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
  mutate(lon = as.numeric(lon))

write_csv(output_df, here::here("output", "3x3null8F.csv"))

#Summarize pop size and frequency by genotype
geno_sum = output_df %>% 
  group_by(lat, lon, generation,genotype,age) %>% 
  summarise(geno_pop_sum = sum(pop)) 

pop_sum = output_df %>% 
  group_by(lat, lon, generation,age) %>% 
  summarise(pop_sum = sum(pop))


output_sum = full_join(geno_sum, pop_sum) %>%
  mutate(freq = geno_pop_sum/pop_sum)

write_csv(output_sum, here("output", "3x3null8F_sum.csv"))

# output_sum = read_csv(here("output", "3x3null8F.csv"))

plot_sum = output_sum %>% 
  filter(generation %in% c(50,70,90,110,130,150)) %>% 
  mutate(generation = as.numeric(generation)) %>% 
  filter(age == "adult")

p1 = ggplot(plot_sum, aes(lon, lat, fill = freq)) +
  geom_tile() + 
  facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Genotype Frequency", color = "Genotype Frequency") +
  theme_bw() +
  scale_fill_gradient2(low = "gainsboro", high = "midnightblue", mid = "skyblue3", midpoint = 0.5)

p2 = ggplot(plot_sum, aes(lon, lat, fill = geno_pop_sum)) +
  geom_tile() + 
  facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Population Size", color = "Population Size") +
  theme_bw() +
  scale_fill_gradient2(low = "gainsboro", high = "midnightblue", mid = "skyblue3", midpoint = 25)

p2 / p1

plot = p2 / p1

ggsave(plot, file=paste0("3x3null8Ftest.pdf"), path = here::here("figs", "test"), height = 11, width = 8)

