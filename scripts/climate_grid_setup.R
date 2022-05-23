library(tidyverse)

addTaskCallback(function(...) {set.seed(42);TRUE})
options(warn=-1)

years = 150
NS.patches = 100
EW.patches = 20
opt.temp = 25

# 0.37      29
# 0.2909605 30
# 0.1690133 31
# 0.0889436 32
# 0.0424048 33
# 0.0183156 34
# 0.007167  35

### UNCOMMENT FOR CONSTANT MEAN SHIFT SST
# SST.patches <- array(0, c(NS.patches, EW.patches, years))
# start_SST = (opt.temp + 4) + NS.patches*0.01
# 
# for (i in 1:years) {
#   SST = start_SST
#   for (lat in 1:NS.patches) {
#     SST.patches[lat,,i] = SST
#     SST = SST - 0.01
#   }
#   start_SST = start_SST + 0.018
# }

### UNCOMMENT FOR ENSO VARIABLE MEAN SST
# SST.patches <- array(0, c(NS.patches, EW.patches, years))
# start_SST = (opt.temp + 4) + NS.patches*0.01
# 
# for (i in 1:years) {
#   SST = start_SST
#   for (lat in 1:NS.patches) {
#     SST.patches[lat,,i] = SST
#     SST = SST - 0.01
#   }
#   start_SST = start_SST + rnorm(1, mean = 0.018, sd = .5)
# }

### UNCOMMENT FOR SHOCK SST CHANGES
# SST.patches <- array(0, c(NS.patches, EW.patches, years))
start_SST = (opt.temp + 4) + NS.patches*0.01

for (i in 1:years) {
  heat_prob = runif(1, 0, 1)
  print(i)
  print(heat_prob)
  if ((i < 75 & heat_prob < 0.1) | (i >= 75 & heat_prob < 0.35)) {
    intensity <- runif(1, .5, ifelse(i < 75, 2, 4))
    SST = start_SST + intensity
  } else {
    SST = start_SST
  }
  for (lat in 1:NS.patches) {
    SST.patches[lat,,i] = SST
    SST = SST - 0.01
  }
}

output_df = data.frame() #create dataframe to hold results

for(b in 1:years) {
  SSTdf = SST.patches[,,b] %>% 
    as.data.frame() 
  
  SSTdf$lat = c(1:NS.patches)
  SSTdf$year = paste0(b)
  
  output_df = bind_rows(output_df, SSTdf)
}

SSTdf = output_df %>% 
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
  # filter(year %in% c(20,40,60,80,100,120,140,150)) %>% 
  # 

#SSTdf$year = factor(SSTdf$year, levels = c(20,40,60,80,100,120,140,150)) #reorder factor year

ggplot(SSTdf, aes(lon, lat, fill = sst)) +
  geom_tile() +
  labs(x = "Longitude", y = "Latitude", fill = "SST") +
  theme_bw() +
  scale_fill_gradient2(low = "white", high = "midnightblue", mid = "lightskyblue", midpoint = 31) +
  facet_wrap(~year) + 
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )