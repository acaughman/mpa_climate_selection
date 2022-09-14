library(tidyverse)
library(here)
library(pracma)
library(beepr)

addTaskCallback(function(...) {
  set.seed(42)
  TRUE
})
options(warn = -1)
options(dplyr.summarise.inform = FALSE)

# Mee Simulation ----------------------------------------------------------

## Parameters:

NUM.reps <- 1 # The number of replicate simulations to run
## 150 years total
NUM.gens.pre.fishing <- 15 # The number of generations before any fishery
NUM.gens.pre.reserve <- 10 # The number of generations of fishing before reserves are installed
NUM.gens.post.reserve <- 150 # The number of generations with the reserve installed
years <- NUM.gens.pre.fishing + NUM.gens.pre.reserve + NUM.gens.post.reserve

NS.patches <- 100 # the number of patches on the north-south axis
EW.patches <- 20 # the number of patches on the east-west axis
patch.size <- 1 # the width and height of each grid cell in nautical miles (COULD BE METERS?)
## View the "world" coordinates:
view.world <- array(seq(1, NS.patches * EW.patches), c(NS.patches, EW.patches))
view.world

init.a <- 0.3 # The initial frequency of the low movement allele
sb <- 0.70 # survival proportion for babies
s <- 0.70 # survival proportion
dd <- 0.002 # density dependence of baby survival
fecundity <- 20000 # The number of babies produced, on average, by each adult female each year.
maturity.age <- 3 # The average age at which individuals mature (i.e., the age at which 50% of individuals are mature)
fished <- 0.7
buffer.fished <- 0.2 # buffer fishing pressure (lower than total = buffer zone, higher than total = fishing the line)

reserves.at <- c(949, 1049, 1149, 950, 1050, 1150, 951, 1051, 1151)
#large_MPA  c(848, 948, 1048, 1148, 1248, 1348, 849, 949, 1049, 1149, 1249, 1349, 850, 950, 1050, 1150, 1250, 1350, 851, 951, 1051, 1151, 1251, 1351, 852, 952, 1052, 1152, 1252, 1352, 853, 953, 1053, 1153, 1253, 1353)
#small_MPA  c(949, 1049, 1149, 950, 1050, 1150, 951, 1051, 1151)
dynamic.reserve <- FALSE

buffer.at <- c()

bold.mover.distance <- 3 # Individuals with AA genotype move this distance on average every year
lazy.mover.distance <- 2 # Individuals with aa genotype move this distance on average every year
Dominance.coefficient <- 0.5 # Dominance coefficient
Heritability.index <- 2 # Influences stochastic variation in movement distance. High numbers decrease variation by reducing the variance around the phenotypic mean in a negative binomial distribution. The phenotypic mean is determined by the genotype.

opt.temp <- 25 # optimal temperature of species
temp.range <- 4 # thermal breath of species

############################################################################
## Create the world

NUM.age.classes <- 3 # babies, juvenile, adult
NUM.sexes <- 2 # female male
NUM.genotypes <- 3 # AA,Aa,aa

world <- array(0, c(NS.patches, EW.patches, NUM.age.classes, NUM.sexes, NUM.genotypes))

############################################################################
## This populates the world.

init <- function() {
  init.AA <- round(100 * (1 - init.a)^2) # initiate the AA
  init.Aa <- round(100 * 2 * (init.a) * (1 - init.a)) # initiate Aa
  init.aa <- round(100 * (init.a)^2) # initiate aa
  pop <- world # assign pop as empty world
  pop[, , , , 1] <- init.AA # add AA to world
  pop[, , , , 2] <- init.Aa # add Aa to word
  pop[, , , , 3] <- init.aa # add aa to world
  return(pop)
}


############################################################################
## This function sets up the Sea surface temperature grid

init_SST <- function(years, climate) {
  if (climate == "null") {
    SST.patches <- array(opt.temp, c(NS.patches, EW.patches, years))
    start_SST <- (opt.temp) + NS.patches * 0.02

    for (i in 1:years) {
      SST <- start_SST
      for (lat in 1:NS.patches) {
        SST.patches[lat, , i] <- SST
        SST <- SST - 0.02
      }
    }
  } else if (climate == "mean") {
    SST.patches <- array(opt.temp, c(NS.patches, EW.patches, years))
    start_SST <- (opt.temp) + NS.patches * 0.02

    for (i in 1:25) {
      SST <- start_SST
      for (lat in 1:NS.patches) {
        SST.patches[lat, , i] <- SST
        SST <- SST - 0.02
      }
    }

    for (i in 26:years) {
      SST <- start_SST
      for (lat in 1:NS.patches) {
        SST.patches[lat, , i] <- SST
        SST <- SST - 0.02
      }
      start_SST <- start_SST + 0.033
    }
  } else if (climate == "enso") {
    SST.patches <- array(opt.temp, c(NS.patches, EW.patches, years))
    start_SST <- (opt.temp) + NS.patches * 0.02

    for (i in 1:25) {
      SST <- start_SST
      for (lat in 1:NS.patches) {
        SST.patches[lat, , i] <- SST
        SST <- SST - 0.02
      }
    }

    t <- seq(1, years - 25, 1)
    b <- array(0, 25)
    change <- 0.5 * sin(t) + 0.033
    enso.value <- c(b, change)

    for (i in 26:years) {
      SST <- start_SST
      for (lat in 1:NS.patches) {
        SST.patches[lat, , i] <- SST
        SST <- SST - 0.02
      }
      start_SST <- start_SST + enso.value[i]
    }
  } else if (climate == "shock") {
    SST.patches <- array(opt.temp, c(NS.patches, EW.patches, years))
    start_SST <- (opt.temp) + NS.patches * 0.02

    for (i in 1:25) {
      SST <- start_SST
      for (lat in 1:NS.patches) {
        SST.patches[lat, , i] <- SST
        SST <- SST - 0.02
      }
    }

    for (i in 26:years) {
      heat_prob <- runif(1, 0, 1)
      if ((i < 75 & heat_prob < 0.1) | (i >= 75 & heat_prob < 0.35)) {
        intensity <- runif(1, 1, ifelse(i < 75, 3, 5))
        SST <- start_SST + intensity
      } else {
        SST <- start_SST
      }
      for (lat in 1:NS.patches) {
        SST.patches[lat, , i] <- SST
        SST <- SST - 0.02
      }
    }
  } else if (climate == "mean shock") {
    SST.patches <- array(opt.temp, c(NS.patches, EW.patches, years))
    start_SST <- (opt.temp) + NS.patches * 0.02

    for (i in 1:25) {
      SST <- start_SST
      for (lat in 1:NS.patches) {
        SST.patches[lat, , i] <- SST
        SST <- SST - 0.02
      }
    }

    for (i in 26:years) {
      SST <- start_SST
      heat_prob <- runif(1, 0, 1)
      if ((i < 75 & heat_prob < 0.1) | (i >= 75 & heat_prob < 0.35)) {
        intensity <- runif(1, 1, ifelse(i < 75, 3, 5))
        SST <- start_SST + intensity
      } else {
        SST <- start_SST
      }
      for (lat in 1:NS.patches) {
        SST.patches[lat, , i] <- SST
        SST <- SST - 0.02
      }
      start_SST <- start_SST + 0.033
    }
  }
  return(SST.patches)
}



############################################################################
## This function creates an array to tell the simulation the reserve locations

where.reserves <- function(reserves.at) {
  if (dynamic.reserve) {
    reserve.patches <- array(0, c(NS.patches, EW.patches, years))
    for (j in 1:years) {
      for (i in 1:length(reserves.at)) {
        x <- ((reserves.at[i] - 1) %/% NS.patches) + 1
        y <- ((reserves.at[i] - 1) %% NS.patches) + 1
        reserve.patches[y, x, j] <- 1
      }
      if (((j %% 3) == 0) & (j >= 40) & (j < 120)) {
        reserves.at <- reserves.at + 3
      }
    }
  } else {
    reserve.patches <- array(0, c(NS.patches, EW.patches, years))
    for (i in 1:length(reserves.at)) {
      x <- ((reserves.at[i] - 1) %/% NS.patches) + 1
      y <- ((reserves.at[i] - 1) %% NS.patches) + 1
      reserve.patches[y, x, ] <- 1
    }
  }
  return(reserve.patches)
}
reserve.patches <- where.reserves(reserves.at)

############################################################################
## This function creates an array to tell the simulation the buffer/fishing the line locations

where.buffer <- function(buffer.at) {
  buffer.patches <- array(0, c(NS.patches, EW.patches))
  for (i in 1:length(buffer.at)) {
    x <- ((buffer.at[i] - 1) %/% NS.patches) + 1
    y <- ((buffer.at[i] - 1) %% NS.patches) + 1
    buffer.patches[y, x] <- 1
  }
  return(buffer.patches)
}
buffer.patches <- where.buffer(buffer.at)

############################################################################
## This function causes adults to reproduce in spawning areas


spawn <- function(pop) {
  fec <- fecundity

  num.males <- rowSums(pop[, , 3, 2, ], dims = 2)

  # All females produce the same mean number of eggs
  NUM.A.eggs <- Reshape(rpois(NS.patches * EW.patches, fec * pop[, , 3, 1, 1] + fec * pop[, , 3, 1, 2] / 2), NS.patches, EW.patches)
  NUM.a.eggs <- Reshape(rpois(NS.patches * EW.patches, fec * pop[, , 3, 1, 3] + fec * pop[, , 3, 1, 2] / 2), NS.patches, EW.patches)
  # Males produce sperm in proportion to their genotypes
  freq.A.sperm <- ifelse(pop[, , 3, 2, 1] == 0, 0, pop[, , 3, 2, 1] / num.males) + ifelse(pop[, , 3, 2, 2] == 0, 0, (pop[, , 3, 2, 2] / num.males) / 2)
  freq.a.sperm <- ifelse(pop[, , 3, 2, 3] == 0, 0, pop[, , 3, 2, 3] / num.males) + ifelse(pop[, , 3, 2, 2] == 0, 0, (pop[, , 3, 2, 2] / num.males) / 2)
  # Sperm fertilize eggs in proportion to sperm genotype frequencies
  AA <- rbinom(NS.patches * EW.patches, NUM.A.eggs, freq.A.sperm)
  aa <- rbinom(NS.patches * EW.patches, NUM.a.eggs, freq.a.sperm)
  Aa <- NUM.A.eggs + NUM.a.eggs - AA - aa
  # Divide zygotes 50:50 among the sexes
  AA.f <- rbinom(NS.patches * EW.patches, AA, 0.5)
  AA.m <- AA - AA.f
  Aa.f <- rbinom(NS.patches * EW.patches, Aa, 0.5)
  Aa.m <- Aa - Aa.f
  aa.f <- rbinom(NS.patches * EW.patches, aa, 0.5)
  aa.m <- aa - aa.f
  # Female babies
  pop[, , 1, 1, 1] <- pop[, , 1, 1, 1] + Reshape(AA.f, NS.patches, EW.patches)
  pop[, , 1, 1, 2] <- pop[, , 1, 1, 2] + Reshape(Aa.f, NS.patches, EW.patches)
  pop[, , 1, 1, 3] <- pop[, , 1, 1, 3] + Reshape(aa.f, NS.patches, EW.patches)
  # Male babies
  pop[, , 1, 2, 1] <- pop[, , 1, 2, 1] + Reshape(AA.m, NS.patches, EW.patches)
  pop[, , 1, 2, 2] <- pop[, , 1, 2, 2] + Reshape(Aa.m, NS.patches, EW.patches)
  pop[, , 1, 2, 3] <- pop[, , 1, 2, 3] + Reshape(aa.m, NS.patches, EW.patches)
  return(pop)
}

############################################################################
## This function calculates temperature based mortality based on sea surface temperature, temperature range and optimal temperature of species

calc_temp_mortality <- function(SST, opt.temp, temp.range, s) {
  nat.m <- array(0, c(nrow(SST), ncol(SST)))
  m <- 1 - exp((-(SST[, ] - opt.temp)^2) / (temp.range^2)) # temperature based mortality function from Walsworth et al.
  m <- 1 - m
  nat.m <- ifelse(m > s, s, m)
  return(nat.m)
}

############################################################################
## This function determines density dependent survival proportion for babies

survival_b <- function(num, SST) {
  s <- calc_temp_mortality(SST, opt.temp, temp.range, sb)
  dd <- dd # density dependence of survival
  result <- s / (1 + dd * num)
  return(result)
}

############################################################################
## This function determines density dependent survival proportion for juveniles and adults

survival <- function(SST) {
  result <- calc_temp_mortality(SST, opt.temp, temp.range, s)
  return(result)
}

############################################################################
## This function determines natural survival and recruitment within each grid cell.

p <- 1 / (maturity.age)

recruit <- function(pop) {
  recruit.array <- world
  SST <- SST.patches[, , t]
  # Some babies survive and recruit to juvenile age class
  s1 <- survival_b(rowSums(pop[, , 1, , ], dim = 2), SST)
  s <- survival(SST)
  # baby dispersal here
  for (j in 1:NUM.sexes) {
    for (k in 1:NUM.genotypes) {
      recruit.array[, , 1 + 1, j, k] <- recruit.array[, , 1 + 1, j, k] + Reshape(rbinom(NS.patches * EW.patches, pop[, , 1, j, k], s1), NS.patches, EW.patches)
    }
  }
  # Some adults survive
  for (j in 1:NUM.sexes) {
    for (k in 1:NUM.genotypes) {
      recruit.array[, , 3, j, k] <- recruit.array[, , 3, j, k] + Reshape(rbinom(NS.patches * EW.patches, pop[, , 3, j, k], s), NS.patches, EW.patches)
    }
  }
  # Some juveniles survive
  for (j in 1:NUM.sexes) {
    for (k in 1:NUM.genotypes) {
      juvies.surviving <- Reshape(rbinom(NS.patches * EW.patches, pop[, , 2, j, k], s), NS.patches, EW.patches)
      # Some juveniles recruit to adult age class
      juvies.recruiting <- Reshape(rbinom(NS.patches * EW.patches, juvies.surviving, p), NS.patches, EW.patches)
      juvies.staying <- juvies.surviving - juvies.recruiting
      recruit.array[, , 2 + 1, j, k] <- recruit.array[, , 2 + 1, j, k] + juvies.recruiting
      # The rest of the juveniles remain in the juvenile age class
      recruit.array[, , 2, j, k] <- recruit.array[, , 2, j, k] + juvies.staying
    }
  }
  return(recruit.array)
}

############################################################################
## This function determines fishing mortality within each grid cell, depending whether the cell is a reserve.

fishing <- function(pop, gen) {
  if (gen <= pre.reserve.gens + pre.fishing.gens) {
    each.patch.pop <- array(0, c(NS.patches, EW.patches))
    for (i in 2:NUM.age.classes) {
      for (j in 1:NUM.sexes) {
        for (k in 1:NUM.genotypes) {
          each.patch.pop[, ] <- each.patch.pop[, ] + pop[, , i, j, k]
        }
      }
    }
    mean.per.patch.pop <- mean(each.patch.pop)
    ff <- mean.per.patch.pop * (1 / fished - 1)
    patch.pop <- rowSums(pop[, , c(2, 3), , ], dims = 2)
    f <- patch.pop / (ff + patch.pop)
    for (i in 2:NUM.age.classes) {
      for (j in 1:NUM.sexes) {
        for (k in 1:NUM.genotypes) {
          pop[, , i, j, k] <- Reshape(rbinom(NS.patches * EW.patches, pop[, , i, j, k], (1 - f)), NS.patches, EW.patches)
        }
      }
    }
  }
  if (gen > pre.reserve.gens + pre.fishing.gens) {
    reserve.area <- sum(reserve.patches[, , gen]) / (NS.patches * EW.patches)
    buffer.area <- sum(buffer.patches) / (NS.patches * EW.patches)
    fished.adj <- (fished - (buffer.area * buffer.fished)) * 1 / (1 - (reserve.area + buffer.area))
    each.patch.pop <- array(0, c(NS.patches, EW.patches))
    for (i in 2:NUM.age.classes) {
      for (j in 1:NUM.sexes) {
        for (k in 1:NUM.genotypes) {
          each.patch.pop[, ] <- each.patch.pop[, ] + pop[, , i, j, k]
        }
      }
    }
    each.patch.pop <- ifelse(reserve.patches[, , gen] == 1 | buffer.patches == 1, NaN, each.patch.pop)
    mean.per.patch.pop <- mean(each.patch.pop, na.rm = TRUE)
    ff <- mean.per.patch.pop * (1 / fished.adj - 1)
    patch.pop <- rowSums(pop[, , c(2, 3), , ], dims = 2)
    patch.pop <- ifelse(reserve.patches[, , gen] == 1 | buffer.patches == 1, NaN, patch.pop)
    f <- patch.pop / (ff + patch.pop)
    if (buffer.fished != 0) {
      f <- ifelse(buffer.patches == 1, buffer.fished, f)
    }
    for (i in 2:NUM.age.classes) {
      for (j in 1:NUM.sexes) {
        for (k in 1:NUM.genotypes) {
          fished.array <- Reshape(rbinom(NS.patches * EW.patches, pop[, , i, j, k], (1 - f)), NS.patches, EW.patches)
          pop[, , i, j, k] <- ifelse(is.na(fished.array[, ]), pop[, , i, j, k], fished.array[, ])
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
  Het.movers <- min(Hom.A.movers, Hom.a.movers) + h * abs(Hom.A.movers - Hom.a.movers)
  Herit <- Heritability.index # Influences heritability of movement. High numbers increase heritability by reducing the variance around the phenotypic mean. The phenotypic mean is determined by the genotype.

  move.array <- world

  for (lat in 1:NS.patches) {
    for (lon in 1:EW.patches) {
      for (i in 2:NUM.age.classes) {
        for (j in 1:NUM.sexes) {
          for (k in 1:NUM.genotypes) {
            if (k == 1) {
              mean.dist <- Hom.A.movers
            }
            if (k == 2) {
              mean.dist <- Het.movers
            }
            if (k == 3) {
              mean.dist <- Hom.a.movers
            }
            # pop[lat,lon,i,j,k] = ifelse(is.na(pop[lat,lon,i,j,k]), 0, pop[lat,lon,i,j,k])
            if (pop[lat, lon, i, j, k] > 0) {
              # movers are subtracted from the present grid cell
              move.array[lat, lon, i, j, k] <- move.array[lat, lon, i, j, k] - pop[lat, lon, i, j, k]
              # determine the distribution of movement distances in nautical miles:
              dist <- rnbinom(pop[lat, lon, i, j, k], mu = mean.dist, size = Herit)
              dist <- ifelse(dist > EW.patches, EW.patches, dist)
              # determine the direction of each move
              theta <- runif(pop[lat, lon, i, j, k], 0, 2 * pi)
              #### bias this movement in the north-south direction (along coasts) if this is a great white shark simulation (otherwise, comment out the next three lines):
              # f.adj <- function(x, u) x-cos(x)*sin(x) - u
              # my.uniroot <- function(x) uniroot(f.adj, c(0, 2*pi), tol = 0.0001, u = x)$root
              # theta <- vapply(theta, my.uniroot, numeric(1))
              # convert direction and distance into a distance in the x-direction (longitude)
              x <- cos(theta) * dist
              # bounce off edges
              x_bool <- ifelse((x <= patch.size * (EW.patches - lon) + patch.size / 2) & (x >= patch.size / 2 - lon * patch.size), TRUE, FALSE)
              x <- ifelse((x_bool == FALSE) & (x > patch.size * (EW.patches - lon) + patch.size / 2), (-(x - 2 * (patch.size * (EW.patches - lon) + patch.size / 2))), x)
              x <- ifelse((x_bool == FALSE) & (x < patch.size / 2 - lon * patch.size), (-(x - 2 * (patch.size / 2 - lon * patch.size))), x)
              # convert direction and distance into a distance in the y-direction (latitude)
              y <- sin(theta) * dist
              #### bounce off edges
              y_bool <- ifelse((y <= patch.size * (NS.patches - lat) + patch.size / 2) & (y >= patch.size / 2 - lat * patch.size), TRUE, FALSE)
              y <- ifelse((y_bool == FALSE) & (y > patch.size * (NS.patches - lat) + patch.size / 2), (-(y - 2 * (patch.size * (NS.patches - lat) + patch.size / 2))), y)
              y <- ifelse((y_bool == FALSE) & (y < patch.size / 2 - lat * patch.size), (-(y - 2 * (patch.size / 2 - lat * patch.size))), y)
              #### Loop NS
              # y_bool = ifelse((y <= patch.size*(NS.patches-lat)+patch.size/2) & (y >= patch.size/2-lat*patch.size), TRUE, FALSE)
              # y = ifelse((y_bool == FALSE) & (y > patch.size*(NS.patches-lat)+patch.size/2),(y - (patch.size * NS.patches)),y)
              # y = ifelse((y_bool == FALSE) & (y < patch.size/2-lat*patch.size),(y + (patch.size * NS.patches)),y)
              # convert movement distances into numbers of grid cells (assume fish start in centre of cell):
              x <- round(x) / patch.size
              y <- round(y) / patch.size
              xy <- as.data.frame(cbind(x, y))
              xy$count <- 1
              freq <- aggregate(count ~ x + y, data = xy, FUN = sum)
              freq2D <- as.data.frame(array(0, c(length(unique(xy$y)), length(unique(xy$x)))))
              names(freq2D) <- sort(unique(xy$x))
              row.names(freq2D) <- sort(unique(xy$y))
              for (row in 1:length(freq$x)) {
                freq2D[as.character(freq$y[row]), as.character(freq$x[row])] <- freq$count[row]
              }
              # populate the move.array with movers (and stayers)
              for (xx in 1:length(unique(xy$x))) {
                for (yy in 1:length(unique(xy$y))) {
                  move.array[lat + as.numeric(row.names(freq2D)[yy]), lon + as.numeric(names(freq2D)[xx]), i, j, k] <- move.array[lat + as.numeric(row.names(freq2D)[yy]), lon + as.numeric(names(freq2D)[xx]), i, j, k] + freq2D[yy, xx]
                }
              }
            }
          }
        }
      }
    }
  }
  # add move.array to pop to finish movement
  return(pop + move.array)
  gc()
}

############################################################################
## THE SIMULATION ##########################################################

reps <- NUM.reps

pre.fishing.gens <- NUM.gens.pre.fishing
pre.reserve.gens <- NUM.gens.pre.reserve
post.reserve.gens <- NUM.gens.post.reserve
gens <- pre.fishing.gens + pre.reserve.gens + post.reserve.gens

output.array <- array(0, c(NS.patches, EW.patches, NUM.age.classes, NUM.sexes, NUM.genotypes, gens, reps))

start_time <- Sys.time()

for (rep in 1:reps) {
  print(rep)
  SST.patches <- init_SST(years, "null") # null, mean, enso, shock, or mean shock
  # save(SST.patches, file = here::here("03_generated_data","climate_layer", "mean_shock.rda"))
  pop <- init()
  for (t in 1:gens) {
    output.array[, , , , , t, rep] <- pop
    pop <- spawn(pop)
    pop <- recruit(pop)
    if (t > pre.fishing.gens) {
      gen <- t
      pop <- fishing(pop, gen)
    }
    pop <- move(pop)
    print(t)
  }
  gc() # clear memory
}
gc()

end_time <- Sys.time()
end_time - start_time

beepr::beep(5)

save(output.array, file = here::here("sensitivity_analysis", "test.rda"))
