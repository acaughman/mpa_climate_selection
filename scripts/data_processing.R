library(tidyverse)
library(patchwork)

# Allie Explore -----------------------------------------------------------

reps = 1
NS.patches <- 100 # the number of patches on the north-south axis
EW.patches <- 20
NUM.age.classes <- 3 #babies, juvenile, adult
NUM.sexes <- 2 #female male
NUM.genotypes <- 3 #AA,Aa,aa
NUM.gens.pre.fishing <- 25 # The number of generations before any fishery
NUM.gens.pre.reserve <- 25 # The number of generations of fishing before reserves are installed
NUM.gens.post.reserve <- 100 # The number of generations with the reserve installed
gens = NUM.gens.pre.fishing+NUM.gens.pre.reserve+NUM.gens.post.reserve

load(file = here::here("data", "tests", "test_mean_noevo.rda"))
load(file = here::here("data","mean.rda"))

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
  mutate(lon = as.numeric(lon)) 

write_csv(output_df, here::here("output", "test_mean_noevo.csv"))

#Summarize pop size and frequency by genotype
geno_sum = output_df %>% 
  filter(rep == 1) %>% #FOR NOW
  group_by(lat, lon, generation,genotype,age) %>% 
  summarise(geno_pop_sum = sum(pop))

pop_sum = output_df %>% 
  group_by(lat, lon, generation,age) %>% 
  summarise(pop_sum = sum(pop), max_temp = max_temp, min_temp = min_temp)

output_sum = full_join(geno_sum, pop_sum) %>%
  mutate(freq = geno_pop_sum/pop_sum)

#write_csv(output_sum, here("output", "3x3null8F_sum.csv"))

# output_sum = read_csv(here("output", "3x3null8F.csv"))

plot_sum = output_sum %>% 
  filter(generation %in% c(25,50,75,100,125,150)) %>% 
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
  scale_fill_gradient2(low = "gainsboro", high = "midnightblue", mid = "skyblue3", midpoint = 20)

p2 / p1

plot = p2 / p1

ggsave(plot, file=paste0("test_mean_noevo_h.pdf"), path = here::here("figs", "test"), height = 11, width = 8)



# Line plots --------------------------------------------------------------



line_df = output_sum %>% 
  group_by(generation, age, genotype) %>% 
  summarise(location_sum = sum(geno_pop_sum),
            max_temp = max(max_temp),
            min_temp = min(min_temp)) %>% 
  filter(age == "adult") %>% 
  mutate(generation = as.numeric(generation)) %>% 
  mutate(min_m = 1 - (1 - exp((-(max_temp - 25)^2)/(4^2)))) %>% 
  mutate(min_survival = ifelse(min_m > 0.59, 0.59, min_m)) %>% 
  mutate(max_m = 1 - (1 - exp((-(min_temp - 25)^2)/(4^2)))) %>% 
  mutate(max_survival = ifelse(max_m > 0.59, 0.59, max_m)) %>% 
  select(-max_m, -min_m) %>% 
  filter(generation > 25)
  
mpa_df = output_sum %>% 
  filter(lon %in% c(9, 10, 11)) %>% 
  filter(lat %in% c(10, 11, 12)) %>% 
  group_by(generation, age, genotype) %>% 
  summarise(location_sum = sum(geno_pop_sum),
            max_temp = max(max_temp),
            min_temp = min(min_temp)) %>% 
  filter(age == "adult") %>% 
  mutate(generation = as.numeric(generation)) %>% 
  mutate(min_m = 1 - (1 - exp((-(max_temp - 25)^2)/(4^2)))) %>% 
  mutate(min_survival = ifelse(min_m > 0.59, 0.59, min_m)) %>% 
  mutate(max_m = 1 - (1 - exp((-(min_temp - 25)^2)/(4^2)))) %>% 
  mutate(max_survival = ifelse(max_m > 0.59, 0.59, max_m)) %>%
  select(-max_m, -min_m)

p3 = ggplot(mpa_df, aes(generation, location_sum)) +
  geom_line() +
  facet_wrap(~genotype, nrow = 3, scales = "free_y") +
  theme_bw() +
  labs(x = "Year", 
       y = "Population Density", 
       color = "Age") +
  geom_vline(xintercept = 11, alpha = 0.3) +
  geom_vline(xintercept = 26, alpha = 0.3) +
  #geom_vline(xintercept = c(136,66,124,79,35,143,130,71,82,90,25,73,80,103,51,34,142,40,86,141,97,137,113,144,122,132,149,129), alpha = 0.3, color= "red") +
  scale_x_continuous(breaks=c(11, 26), labels=c("fishing starts","MPA establishment")) +
  theme(axis.text.x = element_text(angle = 30, hjust=1))
p3

#mean 66, enso 117, shock c(136,66,124,79,35,143,130,71,82,90,25,73,80,103,51,34,142,40,86,141,97,137,113,144,122,132,149,129)

ggsave(p3, file=paste0("test_mean_noevo.pdf"), path = here::here("figs", "test"), height = 11, width = 8)

p4 = ggplot(line_df, aes(generation, location_sum)) +
  geom_line() +
  facet_wrap(~genotype, nrow = 3, scales = "free_y") +
  theme_bw() +
  labs(x = "Year", 
       y = "Population Density", 
       color = "Age") +
  geom_vline(xintercept = 73, alpha = 0.2) +
  scale_x_continuous(breaks=c(73), labels=c("Mortality < 0.2")) +
  theme(axis.text.x = element_text(angle = 60, hjust=1))
p4
