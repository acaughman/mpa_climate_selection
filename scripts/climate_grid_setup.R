library(OpenImageR)

SST.patches <- array(0, c(32, 8, 1))
mean_SST = 28

for (lat in 1:32) {
  for(lon in 1:8) {
    SST = rnorm(1, mean = mean_SST, sd = 0.5)
    SST.patches[lat,lon,1] = SST
  }
}

SSTdf = SST.patches %>% 
  as.data.frame() 

SSTdf$lat = c(1:32)

SSTdf = SSTdf %>% 
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
  ))

ggplot(SSTdf, aes(lon, lat, fill = sst)) +
  geom_tile() +
  labs(x = "Longitude", y = "Latitude", fill = "SST") +
  theme_bw() +
  scale_fill_gradient2(low = "white", high = "midnightblue", mid = "lightskyblue", midpoint = 28)

