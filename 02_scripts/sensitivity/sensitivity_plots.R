library(tidyverse)
library(patchwork)
library(lemon)

# Read in Data ------------------------------------------------------------

df1 <- read_csv(here::here("sensitivity_analysis", "climate_rate", "mean_shock_med_sum.csv")) %>%
  mutate(climate = "Shocks with Mean Shift") %>%
  mutate(dd = as.factor(dd)) %>%
  mutate(max_temp = as.character(max_temp))

df2 <- read_csv(here::here("sensitivity_analysis", "climate_rate", "mean_med_sum.csv")) %>%
  mutate(climate = "Mean Shift") %>%
  mutate(dd = as.factor(dd)) %>%
  mutate(max_temp = as.character(max_temp))

# df3 <- read_csv(here::here("sensitivity_analysis", "climate_rate", "shock_med_sum.csv")) %>%
#   mutate(climate = "Shocks") %>%
#   mutate(dd = as.factor(dd)) %>%
#   mutate(max_temp = as.character(max_temp))

df4 <- read_csv(here::here("sensitivity_analysis", "climate_rate", "enso_med_sum.csv")) %>%
  mutate(climate = "El Nino/La Nina") %>%
  mutate(dd = as.factor(dd)) %>%
  mutate(max_temp = as.character(max_temp))

# df5 <- read_csv(here::here("sensitivity_analysis", "climate_rate", "null_med_sum.csv")) %>%
#   mutate(climate = "Null") %>%
#   mutate(dd = as.factor(dd)) %>%
#   mutate(max_temp = as.character(max_temp))

data <- full_join(df1, df2)
# data <- full_join(data, df3)
# data = full_join(data, df4)
data <- full_join(data, df4) %>%
  mutate(climate = fct_relevel(climate, c("El Nino/La Nina", "Mean Shift", "Shocks with Mean Shift", "Shocks"))) %>% 
  mutate(dd = as.factor(dd)) 

# Biomass Plot ------------------------------------------------------------

p <- ggplot(data, aes(generation, location_sum)) +
  geom_line(aes(color = dd)) +
  facet_wrap(~climate, ncol = 2) +
  theme_bw() +
  labs(
    x = "Year",
    y = "Population Density",
    color = "Rate of Climate Increase Per Year"
  ) +
  geom_vline(xintercept = 16, alpha = 0.3) +
  geom_vline(xintercept = 26, alpha = 0.3) +
  # geom_vline(xintercept = c(136,66,124,79,35,143,130,71,82,90,25,73,80,103,51,34,142,40,86,141,97,137,113,144,122,132,149,129), alpha = 0.3, color= "red") +
  scale_x_continuous(breaks = c(16, 26), labels = c("fishing starts", "MPA establishment")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_viridis_d() +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(fill = "white")) +
  ylim(c(0, 21000))

ggsave(reposition_legend(p, 'center', panel='panel-2-2'), file = paste0("med.pdf"), path = here::here("sensitivity_analysis", "climate_rate"), height = 8, width = 15)
# reposition_legend(p, 'center', panel='panel-3-2')

# Frequency Plot ----------------------------------------------------------

p <- ggplot(data, aes(generation, freq_avg)) +
  geom_line(aes(color = dd)) +
  theme_bw() +
  facet_grid(genotype ~ climate) +
  labs(
    x = "Year",
    y = "Allele Frequency",
    color = "Fishing Pressure"
  ) +
  geom_vline(xintercept = 16, alpha = 0.3) +
  geom_vline(xintercept = 26, alpha = 0.3) +
  # geom_vline(xintercept = c(136,66,124,79,35,143,130,71,82,90,25,73,80,103,51,34,142,40,86,141,97,137,113,144,122,132,149,129), alpha = 0.3, color= "red") +
  scale_x_continuous(breaks = c(16, 26), labels = c("fishing starts", "MPA establishment")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_viridis_d() +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(fill = "white"))

ggsave(p, file = paste0("med.pdf"), path = here::here("sensitivity_analysis", "climate_rate"), height = 8, width = 15)
