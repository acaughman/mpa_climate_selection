library(tidyverse)
library(here)

addTaskCallback(function(...) {set.seed(42);TRUE})
options(warn=-1)

years = 150
NS.patches = 100
EW.patches = 20
opt.temp = 25

calc_temp_mortality <- function(SST, opt.temp, temp.range, s) {
  nat.m = array(0, c(nrow(SST), ncol(SST)))
  m = 1 - exp((-(SST[,] - opt.temp)^2)/(temp.range^2)) # temperature based mortality function from Walsworth et al.
  m = 1 - m
  nat.m = ifelse(m > s, s, m)
  return(nat.m)
}

# 0.3678794 29
# 0.2096114 30
# 0.1053992 31
# 0.0467706 32
# 0.0183156 33


### UNCOMMENT FOR CONSTANT MEAN SHIFT SST
SST.patches.mean <- array(0, c(NS.patches, EW.patches, years))
start_SST = (opt.temp + 2) + NS.patches*0.01

for (i in 1:years) {
  SST = start_SST
  for (lat in 1:NS.patches) {
    SST.patches.mean[lat,,i] = SST
    SST = SST - 0.01
  }
  start_SST = start_SST + 0.018
}

### UNCOMMENT FOR ENSO VARIABLE MEAN SST
SST.patches.enso <- array(0, c(NS.patches, EW.patches, years))
start_SST = (opt.temp + 2) + NS.patches*0.01

t=seq(1,years,1)
enso.value = 0.5*sin(t) + 0.018

for (i in 1:years) {
  SST = start_SST
  for (lat in 1:NS.patches) {
    SST.patches.enso[lat,,i] = SST
    SST = SST - 0.01
  }
  start_SST = start_SST + enso.value[i]
}

### UNCOMMENT FOR SHOCK SST CHANGES
SST.patches.shock <- array(0, c(NS.patches, EW.patches, years))
start_SST = (opt.temp + 2) + NS.patches*0.01

for (i in 1:years) {
  heat_prob = runif(1, 0, 1)
  if ((i < 75 & heat_prob < 0.1) | (i >= 75 & heat_prob < 0.35)) {
    intensity <- runif(1, 1, ifelse(i < 75, 3, 5))
    SST = start_SST + intensity
  } else {
    SST = start_SST
  }
  for (lat in 1:NS.patches) {
    SST.patches.shock[lat,,i] = SST
    SST = SST - 0.01
  }
}

output_df_mean = data.frame() #create dataframe to hold results
output_df_enso = data.frame()
output_df_shock = data.frame()


# Mean --------------------------------------------------------------------

for(b in 1:years) {
  SSTdf_mean = SST.patches.mean[,,b] %>% 
    as.data.frame() 
  
  SSTdf_mean$lat = c(1:NS.patches)
  SSTdf_mean$year = paste0(b)
  
  output_df_mean = bind_rows(output_df_mean, SSTdf_mean)
}

SSTdf_mean = output_df_mean %>% 
  pivot_longer(V1:V20,
               names_to = "lon",
               values_to = "sst") %>% 
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
  mutate(year = as.numeric(year))

mt = ggplot(SSTdf_mean, aes(lon, lat, fill = sst)) +
  geom_tile() +
  labs(x = "Longitude", y = "Latitude", fill = "SST") +
  theme_bw() +
  scale_fill_gradient2(low = "white", high = "midnightblue", mid = "lightskyblue", midpoint = 29) +
  facet_wrap(~year) + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank())

SSTdf_mean = SSTdf_mean %>% 
  mutate(mortality = 1 - (1 - exp((-(sst - 25)^2)/(4^2)))) %>% 
  mutate(survival = ifelse(mortality > 0.59, 0.59, mortality))

mm = ggplot(SSTdf_mean, aes(lon, lat, fill = survival)) +
  geom_tile() +
  labs(x = "Longitude", y = "Latitude", fill = "s") +
  theme_bw() +
  scale_fill_gradient2(low = "black", high = "red", mid = "pink", midpoint = .3) +
  facet_wrap(~year) + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank())

# ENSO --------------------------------------------------------------------

for(b in 1:years) {
  SSTdf_enso = SST.patches.enso[,,b] %>% 
    as.data.frame() 
  
  SSTdf_enso$lat = c(1:NS.patches)
  SSTdf_enso$year = paste0(b)
  
  output_df_enso = bind_rows(output_df_enso, SSTdf_enso)
}

SSTdf_enso = output_df_enso %>% 
  pivot_longer(V1:V20,
               names_to = "lon",
               values_to = "sst") %>% 
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
  mutate(year = as.numeric(year))

et = ggplot(SSTdf_enso, aes(lon, lat, fill = sst)) +
  geom_tile() +
  labs(x = "Longitude", y = "Latitude", fill = "SST") +
  theme_bw() +
  scale_fill_gradient2(low = "white", high = "midnightblue", mid = "lightskyblue", midpoint = 29) +
  facet_wrap(~year) + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank())

SSTdf_enso = SSTdf_enso %>% 
  mutate(mortality = 1 - (1 - exp((-(sst - 25)^2)/(4^2)))) %>% 
  mutate(survival = ifelse(mortality > 0.59, 0.59, mortality))

em = ggplot(SSTdf_enso, aes(lon, lat, fill = survival)) +
  geom_tile() +
  labs(x = "Longitude", y = "Latitude", fill = "Natural Survival Rate") +
  theme_bw() +
  scale_fill_gradient2(low = "black", high = "red", mid = "pink", midpoint = .3) +
  facet_wrap(~year) + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank())

# shock -------------------------------------------------------------------

for(b in 1:years) {
  SSTdf_shock = SST.patches.shock[,,b] %>% 
    as.data.frame() 
  
  SSTdf_shock$lat = c(1:NS.patches)
  SSTdf_shock$year = paste0(b)
  
  output_df_shock = bind_rows(output_df_shock, SSTdf_shock)
}

SSTdf_shock = output_df_shock %>% 
  pivot_longer(V1:V20,
               names_to = "lon",
               values_to = "sst") %>% 
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
  mutate(year = as.numeric(year))

st = ggplot(SSTdf_shock, aes(lon, lat, fill = sst)) +
  geom_tile() +
  labs(x = "Longitude", y = "Latitude", fill = "SST") +
  theme_bw() +
  scale_fill_gradient2(low = "white", high = "midnightblue", mid = "lightskyblue", midpoint = 29) +
  facet_wrap(~year) + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank())

SSTdf_shock = SSTdf_shock %>% 
  mutate(mortality = 1 - (1 - exp((-(sst - 25)^2)/(4^2)))) %>% 
  mutate(survival = ifelse(mortality > 0.59, 0.59, mortality))

sm = ggplot(SSTdf_shock, aes(lon, lat, fill = survival)) +
  geom_tile() +
  labs(x = "Longitude", y = "Latitude", fill = "Natural Survival Rate") +
  theme_bw() +
  scale_fill_gradient2(low = "black", high = "red", mid = "pink", midpoint = .3) +
  facet_wrap(~year) + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank())

ggsave(mt, file=paste0("mean_temp.pdf"), path = here("figs", "climate"), height = 10, width = 15)
ggsave(mm, file=paste0("mean_mortality.pdf"), path = here("figs", "climate"), height = 10, width = 15)
ggsave(et, file=paste0("enso_temp.pdf"), path = here("figs", "climate"), height = 10, width = 15)
ggsave(em, file=paste0("enso_mortality.pdf"), path = here("figs", "climate"), height = 10, width = 15)
ggsave(st, file=paste0("shock_temp.pdf"), path = here("figs", "climate"), height = 10, width = 15)
ggsave(sm, file=paste0("shock_mortality.pdf"), path = here("figs", "climate"), height = 10, width = 15)
