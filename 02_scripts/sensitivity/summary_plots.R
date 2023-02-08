library(tidyverse)
library(patchwork)
library(lemon)


# load data ---------------------------------------------------------------

df1 <- read_csv(here::here("sensitivity_analysis", "density_dependence", "mean_shock_large_sum.csv")) %>%
  mutate(climate = "Mean Shock") %>%
  mutate(dd = as.factor(dd)) %>%
  mutate(max_temp = as.character(max_temp)) %>%
  filter(dd == 0.005) %>%
  mutate(size = "largeium")

df2 <- read_csv(here::here("sensitivity_analysis", "density_dependence", "mean_large_sum.csv")) %>%
  mutate(climate = "Mean") %>%
  mutate(dd = as.factor(dd)) %>%
  mutate(max_temp = as.character(max_temp)) %>%
  filter(dd == 0.005) %>%
  mutate(size = "largeium")

df3 <- read_csv(here::here("sensitivity_analysis", "density_dependence", "shock_large_sum.csv")) %>%
  mutate(climate = "Shock") %>%
  mutate(dd = as.factor(dd)) %>%
  mutate(max_temp = as.character(max_temp)) %>%
  filter(dd == 0.005) %>%
  mutate(size = "largeium")

df4 <- read_csv(here::here("sensitivity_analysis", "density_dependence", "enso_large_sum.csv")) %>%
  mutate(climate = "Enso") %>%
  mutate(dd = as.factor(dd)) %>%
  mutate(max_temp = as.character(max_temp)) %>%
  filter(dd == 0.005) %>%
  mutate(size = "largeium")

df5 <- read_csv(here::here("sensitivity_analysis", "density_dependence", "null_large_sum.csv")) %>%
  mutate(climate = "Null") %>%
  mutate(dd = as.factor(dd)) %>%
  mutate(max_temp = as.character(max_temp)) %>%
  filter(dd == 0.005) %>%
  mutate(size = "largeium")

df6 <- read_csv(here::here("sensitivity_analysis", "density_dependence", "mean_shock_small_sum.csv")) %>%
  mutate(climate = "Mean Shock") %>%
  mutate(dd = as.factor(dd)) %>%
  mutate(max_temp = as.character(max_temp)) %>%
  filter(dd == 0.005) %>%
  mutate(size = "Small")

df7 <- read_csv(here::here("sensitivity_analysis", "density_dependence", "mean_small_sum.csv")) %>%
  mutate(climate = "Mean") %>%
  mutate(dd = as.factor(dd)) %>%
  mutate(max_temp = as.character(max_temp)) %>%
  filter(dd == 0.005) %>%
  mutate(size = "Small")

df8 <- read_csv(here::here("sensitivity_analysis", "density_dependence", "shock_small_sum.csv")) %>%
  mutate(climate = "Shock") %>%
  mutate(dd = as.factor(dd)) %>%
  mutate(max_temp = as.character(max_temp)) %>%
  filter(dd == 0.005) %>%
  mutate(size = "Small")

df9 <- read_csv(here::here("sensitivity_analysis", "density_dependence", "enso_small_sum.csv")) %>%
  mutate(climate = "Enso") %>%
  mutate(dd = as.factor(dd)) %>%
  mutate(max_temp = as.character(max_temp)) %>%
  filter(dd == 0.005) %>%
  mutate(size = "Small")

df10 <- read_csv(here::here("sensitivity_analysis", "density_dependence", "null_small_sum.csv")) %>%
  mutate(climate = "Null") %>%
  mutate(dd = as.factor(dd)) %>%
  mutate(max_temp = as.character(max_temp)) %>%
  filter(dd == 0.005) %>%
  mutate(size = "Small")

df11 <- read_csv(here::here("sensitivity_analysis", "density_dependence", "mean_shock_large_sum.csv")) %>%
  mutate(climate = "Mean Shock") %>%
  mutate(dd = as.factor(dd)) %>%
  mutate(max_temp = as.character(max_temp)) %>%
  filter(dd == 0.005) %>%
  mutate(size = "Large")

df12 <- read_csv(here::here("sensitivity_analysis", "density_dependence", "mean_large_sum.csv")) %>%
  mutate(climate = "Mean") %>%
  mutate(dd = as.factor(dd)) %>%
  mutate(max_temp = as.character(max_temp)) %>%
  filter(dd == 0.005) %>%
  mutate(size = "Large")

df13 <- read_csv(here::here("sensitivity_analysis", "density_dependence", "shock_large_sum.csv")) %>%
  mutate(climate = "Shock") %>%
  mutate(dd = as.factor(dd)) %>%
  mutate(max_temp = as.character(max_temp)) %>%
  filter(dd == 0.005) %>%
  mutate(size = "Large")

df14 <- read_csv(here::here("sensitivity_analysis", "density_dependence", "enso_large_sum.csv")) %>%
  mutate(climate = "Enso") %>%
  mutate(dd = as.factor(dd)) %>%
  mutate(max_temp = as.character(max_temp)) %>%
  filter(dd == 0.005) %>%
  mutate(size = "Large")

df15 <- read_csv(here::here("sensitivity_analysis", "density_dependence", "null_large_sum.csv")) %>%
  mutate(climate = "Null") %>%
  mutate(dd = as.factor(dd)) %>%
  mutate(max_temp = as.character(max_temp)) %>%
  filter(dd == 0.005) %>%
  mutate(size = "Large")

# merge data --------------------------------------------------------------


large <- list(df1, df2, df3, df4, df5) %>% 
  reduce(full_join) %>% 
  select(-dd) %>% 
  filter(generation > 20)
                
small = list(df6, df7, df8, df9, df10)%>% 
  reduce(full_join) %>% 
  select(-dd) %>% 
  filter(generation > 20) 

large = list(df11, df12, df13, df14, df15)%>% 
  reduce(full_join) %>% 
  select(-dd) %>% 
  filter(generation > 20)


# calculate death points --------------------------------------------------

small_zero_ms = small %>% 
  filter(location_sum == 0) %>% 
  filter(climate == "Mean Shock") %>% 
  pull(generation)

position = which((diff(small_zero_ms) == 1) == FALSE)

mean_shock_death_s = small_zero_ms[max(position) + 1]

small_zero_e = small %>% 
  filter(location_sum == 0) %>% 
  filter(climate == "Enso") %>% 
  pull(generation)

position = which((diff(small_zero_e) == 1) == FALSE)

enso_death_s = small_zero_e[1]

small_zero_m = small %>% 
  filter(location_sum == 0) %>% 
  filter(climate == "Mean") %>% 
  pull(generation)

position = which((diff(small_zero_m) == 1) == FALSE)

mean_death_s = small_zero_m[max(position) + 1]


large_zero_ms = large %>% 
  filter(location_sum == 0) %>% 
  filter(climate == "Mean Shock") %>% 
  pull(generation)

position = which((diff(large_zero_ms) == 1) == FALSE)

mean_shock_death_m = large_zero_ms[max(position) + 1]

large_zero_e = large %>% 
  filter(location_sum == 0) %>% 
  filter(climate == "Enso") %>% 
  pull(generation)

position = which((diff(large_zero_e) == 1) == FALSE)

enso_death_m = large_zero_e[max(position) + 1]

large_zero_m = large %>% 
  filter(location_sum == 0) %>% 
  filter(climate == "Mean") %>% 
  pull(generation)

position = which((diff(large_zero_m) == 1) == FALSE)

mean_death_m = large_zero_m[max(position) + 1]


large_zero_ms = large %>% 
  filter(location_sum == 0) %>% 
  filter(climate == "Mean Shock") %>% 
  pull(generation)

position = which((diff(large_zero_ms) == 1) == FALSE)

mean_shock_death_l = large_zero_ms[max(position) + 1]

large_zero_e = large %>% 
  filter(location_sum == 0) %>% 
  filter(climate == "Enso") %>% 
  pull(generation)

position = which((diff(large_zero_e) == 1) == FALSE)

enso_death_l = large_zero_e[max(position) + 1]

large_zero_m = large %>% 
  filter(location_sum == 0) %>% 
  filter(climate == "Mean") %>% 
  pull(generation)

position = which((diff(large_zero_m) == 1) == FALSE)

mean_death_l = large_zero_m[max(position) + 1]

# plot data ---------------------------------------------------------------

p1 = ggplot(small, aes(generation, location_sum)) +
  geom_line(aes(color = climate)) +
  theme_bw() +
  labs(
    title = "3x3",
    x = "Year",
    y = "Population Density",
    color = "Climate Scenario"
) + 
  scale_color_viridis_d() +
  geom_vline(xintercept = mean_shock_death_s, alpha = 0.5, color = "#21918c") +
  geom_vline(xintercept = enso_death_s, alpha = 0.5, color = "#440154") +
  geom_vline(xintercept = mean_death_s, alpha = 0.5, color = "#3b528b")

p2 = ggplot(large, aes(generation, location_sum)) +
  geom_line(aes(color = climate)) +
  theme_bw() +
  labs(
    title = "6x6",
    x = "Year",
    y = "Population Density",
    color = "Climate Scenario"
  ) + 
  scale_color_viridis_d() +
  geom_vline(xintercept = mean_shock_death_m, alpha = 0.5, color = "#21918c") +
  geom_vline(xintercept = enso_death_m, alpha = 0.5, color = "#440154") +
  geom_vline(xintercept = mean_death_m, alpha = 0.5, color = "#3b528b")

p3 = ggplot(large, aes(generation, location_sum)) +
  geom_line(aes(color = climate)) +
  theme_bw() +
  labs(
    title = "12x12",
    x = "Year",
    y = "Population Density",
    color = "Climate Scenario"
  ) + 
  scale_color_viridis_d() +
  geom_vline(xintercept = mean_shock_death_l, alpha = 0.5, color = "#21918c") +
  geom_vline(xintercept = enso_death_l, alpha = 0.5, color = "#440154") +
  geom_vline(xintercept = mean_death_l, alpha = 0.5, color = "#3b528b")


plot = p1 / p2 / p3 + plot_layout(guides = "collect")

ggsave(plot, file = paste0("sums.pdf"), path = here::here("sensitivity_analysis"), height = 8, width = 15)
