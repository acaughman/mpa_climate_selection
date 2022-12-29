library(tidyverse)
library(patchwork)
library(lemon)

# Read in Data ------------------------------------------------------------

df1 <- read_csv(here::here("sensitivity_analysis", "initial_SST", "mean_shock_large_sum.csv")) %>% 
  mutate(climate = "Mean Shock") %>% 
  mutate(dd = as.factor(dd)) %>% 
  mutate(max_temp = as.character(max_temp))

df2 <- read_csv(here::here("sensitivity_analysis", "initial_SST", "mean_large_sum.csv")) %>% 
  mutate(climate = "Mean") %>% 
  mutate(dd = as.factor(dd)) %>% 
  mutate(max_temp = as.character(max_temp))

df3 <- read_csv(here::here("sensitivity_analysis", "initial_SST", "shock_large_sum.csv")) %>% 
  mutate(climate = "Shock") %>% 
  mutate(dd = as.factor(dd)) %>% 
  mutate(max_temp = as.character(max_temp))

df4 <- read_csv(here::here("sensitivity_analysis", "initial_SST", "enso_large_sum.csv")) %>% 
  mutate(climate = "Enso") %>% 
  mutate(dd = as.factor(dd)) %>% 
  mutate(max_temp = as.character(max_temp))

df5 <- read_csv(here::here("sensitivity_analysis", "initial_SST", "null_large_sum.csv")) %>%
  mutate(climate = "Null") %>%
  mutate(dd = as.factor(dd)) %>%
  mutate(max_temp = as.character(max_temp))

data = full_join(df1, df2)
data = full_join(data, df3) 
# data = full_join(data, df4)
data = full_join(data, df4) %>% 
  mutate(climate=fct_relevel(climate,c("Mean", "Mean Shock", "Shock", "Enso"))) %>%
  mutate(dd = as.factor(dd))


# Biomass Plot ------------------------------------------------------------

p <- ggplot(data, aes(generation, location_sum)) +
  geom_line(aes(color = dd)) +
  facet_wrap(~climate) +
  theme_bw() +
  labs(
    x = "Year",
    y = "Population Density",
    color = "Initial SST"
  ) +
  geom_vline(xintercept = 16, alpha = 0.3) +
  geom_vline(xintercept = 26, alpha = 0.3) +
  # geom_vline(xintercept = c(136,66,124,79,35,143,130,71,82,90,25,73,80,103,51,34,142,40,86,141,97,137,113,144,122,132,149,129), alpha = 0.3, color= "red") +
  scale_x_continuous(breaks = c(16, 26), labels = c("fishing starts", "MPA establishment")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_viridis_d() +
  theme(panel.grid.minor = element_blank()) + ylim(c(0, 300000))
  

ggsave(p, file = paste0("large.pdf"), path = here::here("sensitivity_analysis", "initial_SST"), height = 8, width = 15)
#reposition_legend(p, 'center', panel='panel-3-2')

# Frequency Plot ----------------------------------------------------------

p <- ggplot(data, aes(generation, freq_avg)) +
  geom_line(aes(color=dd)) +
  theme_bw() +
  facet_grid(climate~genotype) +
  labs(
    x = "Year",
    y = "Allele Frequency",
    color = "Fishing Pressure"
  ) +
  geom_vline(xintercept = 16, alpha = 0.3) +
  geom_vline(xintercept = 26, alpha = 0.3) +
  # geom_vline(xintercept = c(136,66,124,79,35,143,130,71,82,90,25,73,80,103,51,34,142,40,86,141,97,137,113,144,122,132,149,129), alpha = 0.3, color= "red") +
  scale_x_continuous(breaks = c(16, 26), labels = c("fishing starts", "MPA establishment")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_color_viridis_d() +
  theme(panel.grid.minor = element_blank())

ggsave(p, file = paste0("large.pdf"), path = here::here("sensitivity_analysis", "initial_SST"), height = 8, width = 15)

