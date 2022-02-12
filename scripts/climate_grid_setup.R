library(tidyverse)

years = 15
NS.patches = 32
EW.patches = 8

SST.patches <- array(0, c(NS.patches, EW.patches, years))
start_SST = 28 + 32*.1

for (i in 1:years) {
  SST = start_SST
  for (lat in 1:NS.patches) {
    SST.patches[lat,,i] = SST
    SST = SST - .1
  }
  start_SST = start_SST + rnorm(1, mean = 0.018, sd = 1)
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
  pivot_longer(V1:V8,
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
    lon == "V8" ~ 8
  )) %>% 
  mutate(year = as.factor(year))

ggplot(SSTdf, aes(lon, lat, fill = sst)) +
  geom_tile() +
  labs(x = "Longitude", y = "Latitude", fill = "SST") +
  theme_bw() +
  scale_fill_gradient2(low = "white", high = "midnightblue", mid = "lightskyblue", midpoint = 30) +
  facet_grid(~year)

