library(here)
library(tidyverse)
library(patchwork)
library(beepr)

set.seed(42)

# Mee Simulation ----------------------------------------------------------

## Parameters:

NUM.reps <- 1 # The number of replicate simulations to run
## 150 years total
NUM.gens.pre.fishing <- 150 # The number of generations before any fishery
NUM.gens.pre.reserve <- 0 # The number of generations of fishing before reserves are installed
NUM.gens.post.reserve <- 0 # The number of generations with the reserve installed
years = NUM.gens.pre.fishing+NUM.gens.pre.reserve+NUM.gens.post.reserve

NS.patches <- 120 # the number of patches on the north-south axis
EW.patches <- 8 # the number of patches on the east-west axis
patch.size <- 100 # the width and height of each grid cell in nautical miles (COULD BE METERS?)
## View the "world" coordinates:
view.world <- array(seq(1,NS.patches*EW.patches),c(NS.patches,EW.patches))
view.world

init.a <- 0.5 # The initial frequency of the low movement allele

sb <- 0.37 # survival proportion for babies
s <- 0.37 # survival proportion
dd <- 0.0005 # density dependence of baby survival 
fecundity <- 1500 # The number of babies produced, on average, by each adult female each year.
maturity.age <- 1.5 # The average age at which individuals mature (i.e., the age at which 50% of individuals are mature)
fished.factor <- 0.5
#fished <- fished.factor*(1-s) # Fishing mortalty: the proportion of adults that get fished per year
fished <- fished.factor
buffer.fished <- 0 #buffer fishing pressure (lower than total = buffer zone, higher than total = fishing the line)
reserves.at <- c(102,134,103,135) # This determines which patches are marine reserves. Should be a list: e.g., for one reserve, c(369,370,371,372,389,390,391,392,409,410,411,412,429,430,431,432)
buffer.at <- c()
bold.mover.distance <- 200 # Individuals with AA genotype move this distance on average every year, in nautical miles
lazy.mover.distance <- 100 # Individuals with aa genotype move this distance on average every year, in nautical miles
Dominance.coefficient <- 0.5 # Dominance coefficient
Heritability.index <- 2 # Influences stochastic variation in movement distance. High numbers decrease variation by reducing the variance around the phenotypic mean in a negative binomial distribution. The phenotypic mean is determined by the genotype.
opt.temp = 25 #optimal temperature of species
temp.range = 5 #thermal breath of species

############################################################################
## Create the world

NUM.age.classes <- 3
NUM.sexes <- 2
NUM.genotypes <- 3

world <- array(0, c(NS.patches, EW.patches, NUM.age.classes, NUM.sexes, NUM.genotypes))

############################################################################
## This populates the world.

init <- function() {
  init.AA <- round(200*(1-init.a)^2)
  init.Aa <- round(200*2*(init.a)*(1-init.a))
  init.aa <- round(200*(init.a)^2)
  pop <- world
  for(lat in 1:NS.patches) {
    for(lon in 1:EW.patches) {
      for(i in 1:NUM.age.classes) {
        for(j in 1:NUM.sexes) {
          pop[lat,lon,i,j,1] <- init.AA
          pop[lat,lon,i,j,2] <- init.Aa
          pop[lat,lon,i,j,3] <- init.aa
        }
      }
    }
  }
  return(pop)
}

############################################################################
## This function creates an array to tell the simulation the reserve locations

where.reserves <- function(reserves.at) {
  reserve.patches <- array(0, c(NS.patches, EW.patches))
  for(i in 1:length(reserves.at)) {
    x <- ((reserves.at[i]-1) %/% NS.patches) + 1
    y <- ((reserves.at[i]-1) %% NS.patches) + 1
    reserve.patches[y,x] <- 1
  }
  return(reserve.patches)
}
reserve.patches <- where.reserves(reserves.at)

############################################################################
## This function creates an array to tell the simulation the buffer/fishing the line locations

where.buffer <- function(buffer.at) {
  buffer.patches <- array(0, c(NS.patches, EW.patches))
  for(i in 1:length(buffer.at)) {
    x <- ((buffer.at[i]-1) %/% NS.patches) + 1
    y <- ((buffer.at[i]-1) %% NS.patches) + 1
    buffer.patches[y,x] <- 1
  }
  return(buffer.patches)
}
buffer.patches <- where.buffer(buffer.at)

############################################################################
## This function sets up the Sea surface temperature grid

init_SST <- function(years) {
  
  ### UNCOMMENT FOR CONSISTENT SST
  # SST.patches <- array(25, c(NS.patches, EW.patches, years))
  
  ### UNCOMMENT FOR CONSTANT MEAN SHIFT SST
  SST.patches <- array(0, c(NS.patches, EW.patches, years))
  start_SST = opt.temp + NS.patches*.05
  
  for (i in 1:years) {
    SST = start_SST
    if (SST > 35) {
      SST = 35
    }
    for (lat in 1:NS.patches) {
      SST.patches[lat,,i] = SST
      SST = SST - .05
    }
    start_SST = start_SST + 0.018
  }

  
  ### UNCOMMENT FOR LOW VARIABLE MEAN SST
  # SST.patches <- array(0, c(NS.patches, EW.patches, years))
  # start_SST = opt.temp + NS.patches*.05
  # 
  # for (i in 1:years) {
  #   SST = start_SST
  #   if (SST > 35) {
  #     SST = 35
  #   }
  #   for (lat in 1:NS.patches) {
  #     SST.patches[lat,,i] = SST
  #     SST = SST - .05
  #   }
  #   start_SST = start_SST + rnorm(1, mean = 0.018, sd = 0.01)
  # }
  
  
  ### UNCOMMENT FOR ENSO VARIABLE MEAN SST
  # SST.patches <- array(0, c(NS.patches, EW.patches, years))
  # start_SST = opt.temp + NS.patches*.05
  # 
  # for (i in 1:years) {
  #   SST = start_SST
  #   if (SST > 35) {
  #     SST = 35
  #   }
  #   for (lat in 1:NS.patches) {
  #     SST.patches[lat,,i] = SST
  #     SST = SST - .05
  #   }
  #   start_SST = start_SST + rnorm(1, mean = 0.018, sd = 0.1)
  # }
  # 
  
  return(SST.patches)
  
}

############################################################################
## This function causes adults to reproduce in spawning areas

spawn <- function(pop) {
  
  fec <- fecundity
  
  for(lat in 1:NS.patches) {
    for(lon in 1:EW.patches) {
      num.females <- sum(pop[lat,lon,3,1,])
      num.males <- sum(pop[lat,lon,3,2,])
      # Spawning only occurs if there is at least one males and one females in the patch
      if(num.females > 0 && num.males > 0) {
        # All females produce the same mean number of eggs
        NUM.A.eggs <- rpois(1,fec*pop[lat,lon,3,1,1] + fec*pop[lat,lon,3,1,2]/2)
        NUM.a.eggs <- rpois(1,fec*pop[lat,lon,3,1,3] + fec*pop[lat,lon,3,1,2]/2)
        # Males produce sperm in proportion to their genotypes 
        freq.A.sperm <- pop[lat,lon,3,2,1]/num.males + (pop[lat,lon,3,2,2]/num.males)/2
        freq.a.sperm <- pop[lat,lon,3,2,3]/num.males + (pop[lat,lon,3,2,2]/num.males)/2
        # Sperm fertilize eggs in proportion to sperm genotype frequencies
        AA <- rbinom(1,NUM.A.eggs,freq.A.sperm)
        aa <- rbinom(1,NUM.a.eggs,freq.a.sperm)
        Aa <- NUM.A.eggs+NUM.a.eggs-AA-aa
        # Divide zygotes 50:50 among the sexes
        AA.f <- rbinom(1,AA,0.5)
        AA.m <- AA-AA.f
        Aa.f <- rbinom(1,Aa,0.5)
        Aa.m <- Aa-Aa.f
        aa.f <- rbinom(1,aa,0.5)
        aa.m <- aa-aa.f
        # Female babies
        pop[lat,lon,1,1,1] <- pop[lat,lon,1,1,1] + AA.f
        pop[lat,lon,1,1,2] <- pop[lat,lon,1,1,2] + Aa.f
        pop[lat,lon,1,1,3] <- pop[lat,lon,1,1,3] + aa.f
        # Male babies
        pop[lat,lon,1,2,1] <- pop[lat,lon,1,2,1] + AA.m
        pop[lat,lon,1,2,2] <- pop[lat,lon,1,2,2] + Aa.m
        pop[lat,lon,1,2,3] <- pop[lat,lon,1,2,3] + aa.m
      }
    }
  }
  return(pop)
}

############################################################################
## This function calculates temperature based mortality based on sea surface temperature, temperature range and optimal temperature of species

calc_temp_mortality <- function(SST, opt.temp, temp.range, s) {
  m = 1 - exp((-(SST - opt.temp)^2)/(temp.range^2)) # temperature based mortality function from Walsworth et al.
  m = 1 - m
  if(m > s) {
    nat.m = s
  } else if (m < s) {
    nat.m = m
  }
  return(nat.m)
}

############################################################################
## This function determines density dependent survival proportion for babies

survival_b <- function(num, SST) {
  s = calc_temp_mortality(SST, opt.temp, temp.range, sb)
  dd <- dd # density dependence of survival
  result <- s/(1 + dd * num)
  return(result)
}

############################################################################
## This function determines density dependent survival proportion for juveniles and adults

survival <- function(SST) {
  result = calc_temp_mortality(SST, opt.temp, temp.range, s)
  return(result)
}

############################################################################
## This function determines natural survival and recruitment within each grid cell.

p <- 1/(maturity.age)

recruit <- function(pop) {
  recruit.array <- world
  for(lat in 1:NS.patches) {
    for(lon in 1:EW.patches) {
      SST = SST.patches[lat, lon, t]
      for(i in 1:NUM.age.classes) {
        if(i == 1) {
          # Some babies survive and recruit to juvenile age class
          s1 <- survival_b(sum(pop[lat,lon,i,,]), SST)
          s <- survival(SST)
          for(j in 1:NUM.sexes) {
            for(k in 1:NUM.genotypes) {
              if(pop[lat,lon,i,j,k] > 0) {
                recruit.array[lat,lon,i+1,j,k] <- recruit.array[lat,lon,i+1,j,k] + rbinom(1,pop[lat,lon,i,j,k],s1)
              }
            }
          }
        }
        if(i == 2) {
          # Some juveniles survive
          for(j in 1:NUM.sexes) {
            for(k in 1:NUM.genotypes) {
              if(pop[lat,lon,i,j,k] > 0) {
                juvies.surviving <- rbinom(1,pop[lat,lon,i,j,k],s)
                # Some juveniles recruit to adult age class
                juvies.recruiting <- rbinom(1,juvies.surviving,p)
                recruit.array[lat,lon,i+1,j,k] <- recruit.array[lat,lon,i+1,j,k] + juvies.recruiting
                # The rest of the juveniles remain in the juvenile age class
                recruit.array[lat,lon,i,j,k] <- recruit.array[lat,lon,i,j,k] + juvies.surviving-juvies.recruiting
              }
            }
          }
        }
        if(i == 3) {
          # Some adults survive
          for(j in 1:NUM.sexes) {
            for(k in 1:NUM.genotypes) {
              if(pop[lat,lon,i,j,k] > 0) {
                recruit.array[lat,lon,i,j,k] <- recruit.array[lat,lon,i,j,k] + rbinom(1,pop[lat,lon,i,j,k],s)
              }
            }
          }
        }
      }
    }
  }
  return(recruit.array)
}

############################################################################
## This function determines fishing mortality within each grid cell, depending whether the cell is a reserve.

fishing <- function(pop,gen) {
  if(gen <= pre.reserve.gens+pre.fishing.gens) {
    each.patch.pop <- array(0,c(NS.patches,EW.patches))
    for(i in 2:NUM.age.classes) {
      for(j in 1:NUM.sexes) {
        for(k in 1:NUM.genotypes) {
          each.patch.pop[,] <- each.patch.pop[,] + pop[,,i,j,k]
        }
      }
    }
    mean.per.patch.pop <- mean(each.patch.pop)
    ff <- mean.per.patch.pop*(1/fished-1)
    if(mean.per.patch.pop > 0) {
      for(lat in 1:NS.patches) {
        for(lon in 1:EW.patches) {
          patch.pop <- sum(pop[lat,lon,c(2,3),,])
          if(patch.pop > 0) {
            f <- patch.pop/(ff+patch.pop)
            for(i in 2:NUM.age.classes) {
              for(j in 1:NUM.sexes) {
                for(k in 1:NUM.genotypes) {
                  pop[lat,lon,i,j,k] <- rbinom(1,pop[lat,lon,i,j,k],(1-f))
                }
              }
            }
          }
        }
      }
    }
  }
  if(gen > pre.reserve.gens+pre.fishing.gens) {
    reserve.area <- sum(reserve.patches)/(NS.patches*EW.patches)
    buffer.area <- sum(buffer.patches)/(NS.patches*EW.patches)
    fished.adj <- (fished - (buffer.area*buffer.fished)) * 1/(1-(reserve.area + buffer.area))
    each.patch.pop <- array(0,c(NS.patches,EW.patches))
    for(i in 2:NUM.age.classes) {
      for(j in 1:NUM.sexes) {
        for(k in 1:NUM.genotypes) {
          each.patch.pop[,] <- each.patch.pop[,] + pop[,,i,j,k]
        }
      }
    }
    for(lat in 1:NS.patches) {
      for(lon in 1:EW.patches) {
        if(reserve.patches[lat,lon] == 1) {
          each.patch.pop[lat,lon] <- NaN
        }
      }
    }
    mean.per.patch.pop <- mean(each.patch.pop,na.rm=TRUE)
    ff <- mean.per.patch.pop*(1/fished.adj-1)
    if(mean.per.patch.pop > 0) {
      for(lat in 1:NS.patches) {
        for(lon in 1:EW.patches) {
          if(reserve.patches[lat,lon] == 0 && buffer.patches == 0) {
            patch.pop <- sum(pop[lat,lon,c(2,3),,])
            if(patch.pop > 0) {
              f <- patch.pop/(ff+patch.pop)
              for(i in 2:NUM.age.classes) {
                for(j in 1:NUM.sexes) {
                  for(k in 1:NUM.genotypes) {
                    pop[lat,lon,i,j,k] <- rbinom(1,pop[lat,lon,i,j,k],(1-f))
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return(pop)
}

############################################################################
## This function determines how far each individual moves. Movement distance for each genotype is drawn from a negative bimonial function. Babies do not move between grid cells. 

move <- function(pop) {
  
  h <- Dominance.coefficient # Dominance coefficient
  Hom.A.movers <- bold.mover.distance # Individuals with AA genotype move this distance on average, in nautical miles
  Hom.a.movers <- lazy.mover.distance # Individuals with aa genotype move this distance on average, in nautical miles
  Het.movers <- min(Hom.A.movers,Hom.a.movers) + h * abs(Hom.A.movers - Hom.a.movers)
  Herit <- Heritability.index # Influences heritability of movement. High numbers increase heritability by reducing the variance around the phenotypic mean. The phenotypic mean is determined by the genotype.
  
  move.array <- world
  
  for(lat in 1:NS.patches) {
    for(lon in 1:EW.patches) {
      for(i in 2:NUM.age.classes) {
        for(j in 1:NUM.sexes) {
          for(k in 1:NUM.genotypes) {
            if(k == 1) { mean.dist <- Hom.A.movers }
            if(k == 2) { mean.dist <- Het.movers }
            if(k == 3) { mean.dist <- Hom.a.movers }
            if(pop[lat,lon,i,j,k] > 0) {
              # movers are subtracted from the present grid cell
              move.array[lat,lon,i,j,k] <- move.array[lat,lon,i,j,k] - pop[lat,lon,i,j,k]
              # determine the distribution of movement distances in nautical miles:
              dist <- rnbinom(pop[lat,lon,i,j,k], mu = mean.dist, size = Herit)
              # determine the direction of each move
              theta <- runif(pop[lat,lon,i,j,k],0,2*pi)
              # bias this movement in the north-south direction (along coasts) if this is a great white shark simulation (otherwise, comment out the next three lines):
              f.adj <- function(x, u) x-cos(x)*sin(x) - u
              my.uniroot <- function(x) uniroot(f.adj, c(0, 2*pi), tol = 0.0001, u = x)$root
              theta <- vapply(theta, my.uniroot, numeric(1))
              # convert direction and distance into a distance in the x-direction (longitude)
              x <- cos(theta)*dist
              # bounce off edges (assume fish start in centre of cell)
              for(m in 1:length(x)) {
                miss_x.edges <- FALSE
                while(miss_x.edges==FALSE) {
                  if(x[m] <= patch.size*(EW.patches-lon)+patch.size/2) {
                    if(x[m] >= patch.size/2-lon*patch.size) {
                      miss_x.edges <- TRUE }
                  }
                  if(x[m] > patch.size*(EW.patches-lon)+patch.size/2) {
                    x[m] <- -(x[m]-2*(patch.size*(EW.patches-lon)+patch.size/2))
                    # distance penalty for hitting an edge
                    #x[m] <- x[m] + 1
                  }
                  if(x[m] < patch.size/2-lon*patch.size) {
                    x[m] <- -(x[m]-2*(patch.size/2-lon*patch.size))
                    # distance penalty for hitting an edge
                    #x[m] <- x[m] - 1
                  }
                }
              }
              # convert direction and distance into a distance in the y-direction (latitude)
              y <- sin(theta)*dist
              # bounce off edges (assume fish start in centre of cell)
              for(m in 1:length(y)) {
                miss_y.edges <- FALSE
                while(miss_y.edges==FALSE) {
                  if(y[m] <= patch.size*(NS.patches-lat)+patch.size/2) {
                    if(y[m] >= patch.size/2-lat*patch.size) {
                      miss_y.edges <- TRUE }
                  }
                  if(y[m] > patch.size*(NS.patches-lat)+patch.size/2) {
                    y[m] <- -(y[m]-2*(patch.size*(NS.patches-lat)+patch.size/2))
                    # distance penalty for hitting an edge
                    #y[m] <- y[m] + 1
                  }
                  if(y[m] < patch.size/2-lat*patch.size) {
                    y[m] <- -(y[m]-2*(patch.size/2-lat*patch.size))
                    # distance penalty for hitting an edge
                    #y[m] <- y[m] - 1
                  }
                }
              }
              # convert movement distances into numbers of grid cells (assume fish start in centre of cell):
              xy <- as.data.frame(cbind(x,y))
              hx <- hist(xy$x, breaks = seq(round(min(x),-2)-patch.size/2,round(max(x),-2)+patch.size/2,patch.size), plot = FALSE)
              xbreaks <- hx$breaks
              xmids <- hx$mids
              hy <- hist(xy$y, breaks = seq(round(min(y),-2)-patch.size/2,round(max(y),-2)+patch.size/2,patch.size), plot = FALSE)
              ybreaks <- hy$breaks
              ymids <- hy$mids
              freq <-  as.data.frame(table(cut(xy[,2], ybreaks,labels=ymids),cut(xy[,1], xbreaks,labels=xmids)))
              freq2D <- as.data.frame(array(0,c(length(ymids),length(xmids))))
              names(freq2D) <- xmids/patch.size
              row.names(freq2D) <- ymids/patch.size
              freq2D[cbind(freq[,1], freq[,2])] <- freq[,3]
              # populate the move.array with movers (and stayers)
              for(xx in 1:length(xmids)) {
                for(yy in 1:length(ymids)) {
                  move.array[lat+as.numeric(row.names(freq2D)[yy]),lon+as.numeric(names(freq2D)[xx]),i,j,k] <- move.array[lat+as.numeric(row.names(freq2D)[yy]),lon+as.numeric(names(freq2D)[xx]),i,j,k] + freq2D[yy,xx]
                }
              }
            }
          }
        }
      }
    }
  }
  # add move.array to pop to finish movement
  return(pop+move.array)
}

############################################################################
## THE SIMULATION ##########################################################

reps <- NUM.reps

pre.fishing.gens <- NUM.gens.pre.fishing
pre.reserve.gens <- NUM.gens.pre.reserve
post.reserve.gens <- NUM.gens.post.reserve
gens <- pre.fishing.gens+pre.reserve.gens+post.reserve.gens

output.array <- array(0 ,c(NS.patches, EW.patches, NUM.age.classes, NUM.sexes, NUM.genotypes, gens, reps))

for(rep in 1:reps) {
  #print(rep)
  pop <- init()
  SST.patches <- init_SST(gens)
  for(t in 1:gens) {
    #print(SST.patches[,,t])
    output.array[,,,,,t,rep] <- pop
    pop <- spawn(pop)
    pop <- recruit(pop)
    if(t > pre.fishing.gens) {
      gen <- t
      pop <- fishing(pop,gen)
    }
    pop <- move(pop)
    print(t)
  }
}

beep(5)


# Allie Explore -----------------------------------------------------------

# Output results into a dataframe
output_df = data.frame() #create dataframe to hold results

for(a in 1:reps) {
  world_sub <- array(0, c(NS.patches, EW.patches))
  for(b in 1:gens) {
    for(c in 1:NUM.genotypes) {
      for(d in 1:NUM.sexes) {
        for(e in 1:NUM.age.classes) {
          world_sub = output.array[,,e,d,c,b,a] %>% 
            as.data.frame()
          world_sub$rep = paste0(a)
          world_sub$generation = paste0(b)
          world_sub$genotype = paste0(c)
          world_sub$sex = paste0(d)
          world_sub$age = paste0(e)
          world_sub$lat = c(1:NS.patches)
          output_df = bind_rows(output_df, world_sub)
        }
      }
    }
  }
}

# Wrangle dataframe into plottable format
output_df = output_df %>% 
  pivot_longer(V1:V8,
               names_to = "lon",
               values_to = "pop") %>% 
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
  mutate(genotype = case_when(
    genotype == 1 ~ "AA",
    genotype == 2 ~ "Aa",
    genotype == 3 ~ "aa"
  )) %>% 
  mutate(genotype = as.factor(genotype)) %>% 
  mutate(lat = as.numeric(lat)) %>% 
  mutate(lon = as.numeric(lon))


#Summarize pop size and frequency by genotype
geno_sum = output_df %>% 
  group_by(lat, lon, rep, generation,genotype) %>% 
  summarise(geno_pop_sum = sum(pop)) 

pop_sum = output_df %>% 
  group_by(lat, lon, rep, generation) %>% 
  summarise(pop_sum = sum(pop))


output_sum = full_join(geno_sum, pop_sum) %>%
  mutate(freq = geno_pop_sum/pop_sum) 

#write_csv(output_sum, here("test_output", "climate_test.csv"))


plot_sum = output_sum %>% 
  filter(generation %in% c(20,40,60,80,100,120,140,160,180,200)) %>% 
  mutate(generation = as.numeric(generation))

p1 = ggplot(plot_sum, aes(lon, lat, fill = freq)) +
  geom_tile() + 
  facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Genotype Frequency", color = "Genotype Frequency") +
  theme_bw() +
  scale_fill_gradient2(low = "white", high = "midnightblue", mid = "lightskyblue", midpoint = 0.4)

p2 = ggplot(plot_sum, aes(lon, lat, fill = geno_pop_sum)) +
  geom_tile() + 
  facet_grid(genotype~generation) + 
  labs(x = "Longitude", y = "Latitude", fill = "Population Size", color = "Population Size") +
  theme_bw() +
  scale_fill_gradient2(low = "white", high = "midnightblue", mid = "lightskyblue", midpoint = 300)

p2 / p1

plot = p2 / p1

#ggsave(plot, file=paste0("2x28FConstantTopt.pdf"), path = here("figs"))
