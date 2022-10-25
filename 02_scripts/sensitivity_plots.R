library(tidyverse)
library(patchwork)

# Read in Data ------------------------------------------------------------

df1 <- read_csv(here::here("sensitivity_analysis", "fishing_pressure", "mean_shock_small_sum.csv")) %>% 
  mutate(climate = "mean shock") %>% 
  mutate(dd = as.factor(dd)) %>% 
  mutate(max_temp = as.character(max_temp))

df2 <- read_csv(here::here("sensitivity_analysis", "fishing_pressure", "mean_small_sum.csv")) %>% 
  mutate(climate = "mean") %>% 
  mutate(dd = as.factor(dd)) %>% 
  mutate(max_temp = as.character(max_temp))

df3 <- read_csv(here::here("sensitivity_analysis", "fishing_pressure", "shock_small_sum.csv")) %>% 
  mutate(climate = "shock") %>% 
  mutate(dd = as.factor(dd)) %>% 
  mutate(max_temp = as.character(max_temp))

df4 <- read_csv(here::here("sensitivity_analysis", "fishing_pressure", "enso_small_sum.csv")) %>% 
  mutate(climate = "enso") %>% 
  mutate(dd = as.factor(dd)) %>% 
  mutate(max_temp = as.character(max_temp))

df5 <- read_csv(here::here("sensitivity_analysis", "fishing_pressure", "null_small_sum.csv")) %>% 
  mutate(climate = "null") %>% 
  mutate(dd = as.factor(dd)) %>% 
  mutate(max_temp = as.character(max_temp))

data = full_join(df1, df2)
data = full_join(data, df3) 
data = full_join(data, df4) 
data = full_join(data, df5) 


# Biomass Plot ------------------------------------------------------------


p <- ggplot(data, aes(generation, location_sum)) +
  geom_line(aes(color = dd)) +
  facet_wrap(~climate) +
  theme_bw() +
  labs(
    x = "Year",
    y = "Population Density",
    color = "Density Dependence"
  ) +
  geom_vline(xintercept = 16, alpha = 0.3) +
  geom_vline(xintercept = 26, alpha = 0.3) +
  # geom_vline(xintercept = c(136,66,124,79,35,143,130,71,82,90,25,73,80,103,51,34,142,40,86,141,97,137,113,144,122,132,149,129), alpha = 0.3, color= "red") +
  scale_x_continuous(breaks = c(16, 26), labels = c("fishing starts", "MPA establishment")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_color_viridis_d() +
  theme(panel.grid.minor = element_blank()) +ylim(c(0, 10000))

ggsave(p, file = paste0("small.pdf"), path = here::here("sensitivity_analysis", "fishing_pressure"), height = 8, width = 15)


# Frequency Plot ----------------------------------------------------------

p <- ggplot(data, aes(generation, freq_avg)) +
  geom_line(aes(color=dd)) +
  theme_bw() +
  facet_grid(climate~genotype) +
  labs(
    x = "Year",
    y = "Population Density",
    color = "Fishing Pressure"
  ) +
  geom_vline(xintercept = 16, alpha = 0.3) +
  geom_vline(xintercept = 26, alpha = 0.3) +
  # geom_vline(xintercept = c(136,66,124,79,35,143,130,71,82,90,25,73,80,103,51,34,142,40,86,141,97,137,113,144,122,132,149,129), alpha = 0.3, color= "red") +
  scale_x_continuous(breaks = c(16, 26), labels = c("fishing starts", "MPA establishment")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_color_viridis_d() +
  theme(panel.grid.minor = element_blank())

ggsave(p, file = paste0("small.pdf"), path = here::here("sensitivity_analysis", "fishing_pressure"), height = 8, width = 15)

