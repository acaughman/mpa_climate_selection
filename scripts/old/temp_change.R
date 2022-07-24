library(raster)
library(ncdf4)
library(tidyverse)

# 1 degree latitude = 111.32 km

# Raster Data (cds.climate.copernicus.eu) -------------------------------------------------------------

r1 = raster(here::here("climate_data", "20161201120000-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_CDR2.1-v02.0-fv01.0.nc"))
r2 = raster(here::here("climate_data", "20160701120000-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_CDR2.1-v02.0-fv01.0.nc"))

df1 = rasterToPoints(r1) %>% 
  as.data.frame() %>% 
  mutate(temp = -273.15 + analysed.sea.surface.temperature) 

df2 = rasterToPoints(r2) %>% 
  as.data.frame() %>% 
  mutate(temp = -273.15 + analysed.sea.surface.temperature) 

df1.1 = sample_n(df1, size=1000) %>% 
  filter(y > -50 & y < -30)
df1.2 = sample_n(df1, size=1000) %>% 
  filter(y < 50 & y > 30)

df2.1 = sample_n(df2, size=1000) %>% 
  filter(y > -50 & y < -30)
df2.2 = sample_n(df2, size=1000) %>% 
  filter(y < 50 & y > 30)

mod1.1 = lm(temp ~ y, data = df1.1)
mod2.1 = lm(temp ~ y, data = df2.1)
mod1.2 = lm(temp ~ y, data = df1.2)
mod2.2 = lm(temp ~ y, data = df2.2)

coef1.1 = mod1.1$coefficients %>% 
  as.data.frame()
coef2.1 = mod2.1$coefficients %>% 
  as.data.frame()
coef1.2 = mod1.2$coefficients %>% 
  as.data.frame()
coef2.2 = mod2.2$coefficients %>% 
  as.data.frame()

slope1.1 = coef1.1$.[2]
slope2.1 = coef2.1$.[2]
slope1.2 = abs(coef1.2$.[2])
slope2.2 = abs(coef2.2$.[2])

dkm = 111.12
dt = median(c(slope1.1,slope1.2,slope2.1,slope2.2))

temp_per_km =  dt/dkm #0.007611381
temp_per_100m = temp_per_km/1000*100 #0.0007611381

ggplot(sample_n(df1, size=1000), aes(y, temp)) +
  geom_point()

ggplot(sample_n(df2, size=1000), aes(y, temp)) +
  geom_point()

