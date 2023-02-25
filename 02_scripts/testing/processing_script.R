library(tidyverse)

output_df = read_csv(here::here("03_generated_data", "model_outputs", "model31.csv"))

# Summarize pop size and frequency by genotype
geno_sum2 <- output_df %>%
  filter(age == "adult") %>% 
  group_by(lat, lon, generation, rep, genotype) %>%
  summarise(geno_pop_sum = sum(pop, na.rm = TRUE))

geno_mean2 <- geno_sum2 %>% 
  group_by(lat, lon, generation, genotype) %>%
  summarise(geno_pop_mean = mean(geno_pop_sum, na.rm=TRUE),
            geno_pop_sd = sd(geno_pop_sum, na.rm=TRUE)) %>% 
  ungroup() %>% 
  mutate(geno_pop_se = geno_pop_sd / sqrt(10))

pop_sum2 <- output_df %>%
  filter(age == "adult") %>% 
  group_by(lat, lon, generation, rep) %>%
  summarise(pop_sum = sum(pop, na.rm = TRUE), 
            max_temp = max_temp, 
            min_temp = min_temp,
            fished_sum = sum(fished))

pop_mean2 <- pop_sum2 %>%
  group_by(lat, lon, generation) %>%
  summarise(pop_mean = mean(pop_sum, na.rm = TRUE), 
            pop_sd = sd(pop_sum, na.rm = TRUE),
            max_temp = max_temp, 
            min_temp = min_temp,
            fished_mean = mean(fished_sum, na.rm=TRUE)) %>% 
  ungroup() %>% 
  mutate(pop_se = pop_sd / sqrt(10))

output_sum2 <- full_join(geno_mean2, pop_mean2)

output_sum2 = output_sum2 %>%
  distinct() %>% 
  mutate(freq = geno_pop_mean / pop_mean)

write_csv(output_sum2, here::here("summary_m31.csv"))
