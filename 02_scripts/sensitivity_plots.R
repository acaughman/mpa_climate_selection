
# Read in Data ------------------------------------------------------------



# Biomass Plot ------------------------------------------------------------


p <- ggplot(data, aes(generation, location_sum)) +
  geom_line(aes(color = dd)) +
  theme_bw() +
  labs(
    x = "Year",
    y = "Population Density",
    color = "Movement Pattern"
  ) +
  geom_vline(xintercept = 16, alpha = 0.3) +
  geom_vline(xintercept = 26, alpha = 0.3) +
  # geom_vline(xintercept = c(136,66,124,79,35,143,130,71,82,90,25,73,80,103,51,34,142,40,86,141,97,137,113,144,122,132,149,129), alpha = 0.3, color= "red") +
  scale_x_continuous(breaks = c(16, 26), labels = c("fishing starts", "MPA establishment")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_color_viridis_d() +
  theme(panel.grid.minor = element_blank())

ggsave(p, file = paste0("null_med.pdf"), path = here::here("sensitivity_analysis", "movement_pattern"), height = 8, width = 15)


# Frequency Plot ----------------------------------------------------------

p <- ggplot(data, aes(generation, freq_avg)) +
  geom_line(aes(color=dd)) +
  theme_bw() +
  facet_wrap(~genotype) +
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

ggsave(p, file = paste0("null_large.pdf"), path = here::here("sensitivity_analysis", "fishing_pressure"), height = 8, width = 15)


# Nested Plot -------------------------------------------------------------

p <- ggplot(data, aes(generation, location_sum)) +
  geom_line(aes(color = dd)) +
  theme_bw() +
  labs(
    x = "Year",
    y = "Population Density",
    color = "Movement Pattern"
  ) +
  geom_vline(xintercept = 16, alpha = 0.3) +
  geom_vline(xintercept = 26, alpha = 0.3) +
  # geom_vline(xintercept = c(136,66,124,79,35,143,130,71,82,90,25,73,80,103,51,34,142,40,86,141,97,137,113,144,122,132,149,129), alpha = 0.3, color= "red") +
  scale_x_continuous(breaks = c(16, 26), labels = c("fishing starts", "MPA establishment")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_color_viridis_d() +
  theme(panel.grid.minor = element_blank())

sub <- ggplot(data, aes(generation, location_sum)) +
  geom_line(aes(color = dd)) +
  theme_bw() +
  labs(
    x = "",
    y = "",
    color = ""
  ) +
  geom_vline(xintercept = 16, alpha = 0.3) +
  geom_vline(xintercept = 26, alpha = 0.3) +
  # geom_vline(xintercept = c(136,66,124,79,35,143,130,71,82,90,25,73,80,103,51,34,142,40,86,141,97,137,113,144,122,132,149,129), alpha = 0.3, color= "red") +
  theme(axis.text = element_blank()) +
  theme(axis.title = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(axis.line = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_viridis_d() +
  ylim(c(0, 100000)) +
  theme(legend.position = "none")

plot <- p + inset_element(sub,
                          right = .95,
                          bottom = 0.5,
                          left = 0.3,
                          top = 0.96
) + plot_annotation(tag_levels = "a")

ggsave(plot, file = paste0("null_med.pdf"), path = here::here("sensitivity_analysis", "movement_pattern"), height = 8, width = 15)
