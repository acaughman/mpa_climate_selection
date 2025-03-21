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

NUM.reps <- 10 # The number of replicate simulations to run
## 175 years total
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
dd <- 0.005 # density dependence of baby survival
fecundity <- 20000 # The number of babies produced, on average, by each adult female each year.
maturity.age <- 3 # The average age at which individuals mature (i.e., the age at which 50% of individuals are mature)
fished <- 0.7 # base fishing pressure
buffer.fished <- 0 # buffer fishing pressure (lower than total = buffer zone, higher than total = fishing the line)

reserves.at <- c(849, 949, 1049, 850, 950, 1050, 851, 951, 1051, 829, 929, 1029, 928, 828, 1028, 827, 927, 1027, 871, 971, 1071, 872, 972, 1072, 873, 973, 1073)
# large MPA c(446, 546, 646, 746, 846, 946, 1046, 1146, 1246, 1346, 447, 547, 647, 747, 847, 947, 1047, 1147, 1247, 1347, 448, 548, 648, 748, 848, 948, 1048, 1148, 1248, 1348, 449, 549, 649, 749, 849, 949, 1049, 1149, 1249, 1349, 450, 550, 650, 750, 850, 950, 1050, 1150, 1250, 1350, 451, 551, 651, 751, 851, 951, 1051, 1151, 1251, 1351, 452, 552, 652, 752, 852, 952, 1052, 1152, 1252, 1352, 453, 553, 653, 753, 853, 953, 1053, 1153, 1253, 1353, 454, 554, 654, 754, 854, 954, 1054, 1154, 1254, 1354, 455, 555, 655, 755, 855, 955, 1055, 1155, 1255, 1355)
# medium MPA c(847, 947, 1047, 1147, 1247, 1347, 848, 948, 1048, 1148, 1248, 1348, 849, 949, 1049, 1149, 1249, 1349, 850, 950, 1050, 1150, 1250, 1350, 851, 951, 1051, 1151, 1251, 1351, 852, 952, 1052, 1152, 1252, 1352)
# small MPA c(849, 949, 1049, 850, 950, 1050, 851, 951, 1051)
# tiny MPA c(950)

buffer.at <- c()

bold.mover.distance <- 3 # Individuals with AA genotype move this distance on average every year
lazy.mover.distance <- 1 # Individuals with aa genotype move this distance on average every year
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
  # for Null climate
  if (climate == "null") {
    SST.patches <- array(opt.temp, c(NS.patches, EW.patches, years)) # initiate array for SST
    start_SST <- (opt.temp) + NS.patches * 0.02 # have the highest SST be 0.2 * NS.patches

    for (i in 1:years) { # loop through years
      SST <- start_SST
      for (lat in 1:NS.patches) {
        SST.patches[lat, , i] <- SST
        SST <- SST - 0.02 # subtract 0.02 for each latitude
      }
    }
    # for Mean Shifts
  } else if (climate == "mean") {
    SST.patches <- array(opt.temp, c(NS.patches, EW.patches, years)) # initiate array for SST
    start_SST <- (opt.temp) + NS.patches * 0.02 # have the highest SST be 0.2 * NS.patches

    for (i in 1:25) { # loop through pre MPA years
      SST <- start_SST
      for (lat in 1:NS.patches) {
        SST.patches[lat, , i] <- SST
        SST <- SST - 0.02 # subtract 0.02 for each latitude
      }
    }

    for (i in 26:years) { # loop through MPA years
      SST <- start_SST
      for (lat in 1:NS.patches) {
        SST.patches[lat, , i] <- SST
        SST <- SST - 0.02 # subtract 0.02 for each latitude
      }
      start_SST <- start_SST + 0.033 # add 0.033 degree C every year
    }

    # for El Nino La Nina
  } else if (climate == "enso") {
    SST.patches <- array(opt.temp, c(NS.patches, EW.patches, years)) # initiate array for SST
    start_SST <- (opt.temp) + NS.patches * 0.02 # have the highest SST be 0.2 * NS.patches

    for (i in 1:25) { # loop through pre MPA years
      SST <- start_SST
      for (lat in 1:NS.patches) {
        SST.patches[lat, , i] <- SST
        SST <- SST - 0.02 # subtract 0.02 for each latitude
      }
    }

    t <- seq(1, years - 25, 1) # set time array
    b <- array(0, 25) # set pre MPA times
    change <- 0.5 * sin(t) + 0.033 # calculate temp change per year based on sinusoidal equation
    enso.value <- c(b, change) # append no change pre MPA years to ENSO years

    for (i in 26:years) { # loop through MPA years
      SST <- start_SST
      for (lat in 1:NS.patches) {
        SST.patches[lat, , i] <- SST
        SST <- SST - 0.02 # subtract 0.02 for each latitude
      }
      start_SST <- start_SST + enso.value[i] # add change in temp for each year based on ENSO
    }
    # for Shocks
  } else if (climate == "shock") {
    SST.patches <- array(opt.temp, c(NS.patches, EW.patches, years)) # initiate array for SST
    start_SST <- (opt.temp) + NS.patches * 0.02 # have the highest SST be 0.2 * NS.patches

    for (i in 1:25) { # loop through pre MPA years
      SST <- start_SST
      for (lat in 1:NS.patches) {
        SST.patches[lat, , i] <- SST
        SST <- SST - 0.02 # subtract 0.02 for each latitude
      }
    }

    for (i in 26:years) { # loop through MPA years
      heat_prob <- runif(1, 0, 1) # draw heatwave probability
      if ((i < 75 & heat_prob < 0.1) | (i >= 75 & heat_prob < 0.35)) {
        intensity <- runif(1, 1, 4) # if heatwave, draw intensity
        SST <- start_SST + intensity # add intensity to SST
      } else {
        SST <- start_SST
      }
      for (lat in 1:NS.patches) {
        SST.patches[lat, , i] <- SST
        SST <- SST - 0.02 # subtract 0.02 for each latitude
      }
    }
    # for Shocks with Mean Shift
  } else if (climate == "mean shock") {
    SST.patches <- array(opt.temp, c(NS.patches, EW.patches, years)) # initiate array for SST
    start_SST <- (opt.temp) + NS.patches * 0.02 # have the highest SST be 0.2 * NS.patches

    for (i in 1:25) { # loop through pre MPA years
      SST <- start_SST
      for (lat in 1:NS.patches) {
        SST.patches[lat, , i] <- SST
        SST <- SST - 0.02 # subtract 0.02 for each latitude
      }
    }

    for (i in 26:years) { # loop through MPA years
      SST <- start_SST
      heat_prob <- runif(1, 0, 1) # draw heatwave probability
      if ((i < 75 & heat_prob < 0.1) | (i >= 75 & heat_prob < 0.35)) {
        intensity <- runif(1, 1, 4) # if heatwave, draw intensity
        SST <- start_SST + intensity # add intensity to SST
      } else {
        SST <- start_SST
      }
      for (lat in 1:NS.patches) {
        SST.patches[lat, , i] <- SST
        SST <- SST - 0.02 # subtract 0.02 for each latitude
      }
      start_SST <- start_SST + 0.033 # add 0.033 degree C every year
    }
  }
  return(SST.patches)
}

############################################################################
## This function creates an array to tell the simulation the reserve locations

where.reserves <- function(reserves.at) {
  reserve.patches <- array(0, c(NS.patches, EW.patches, years)) # get numbers for reserve
  for (i in 1:length(reserves.at)) {
    x <- ((reserves.at[i] - 1) %/% NS.patches) + 1 # get corresponding x values
    y <- ((reserves.at[i] - 1) %% NS.patches) + 1 # get corresponding y values
    reserve.patches[y, x, ] <- 1 # assign 1 to indicate reserve
  }
  return(reserve.patches)
}
reserve.patches <- where.reserves(reserves.at) # initiate reserves

############################################################################
## This function creates an array to tell the simulation the buffer/fishing the line locations

where.buffer <- function(buffer.at) {
  buffer.patches <- array(0, c(NS.patches, EW.patches)) # get numbers for buffer
  for (i in 1:length(buffer.at)) {
    x <- ((buffer.at[i] - 1) %/% NS.patches) + 1 # get corresponding x values
    y <- ((buffer.at[i] - 1) %% NS.patches) + 1 # get corresponding y values
    buffer.patches[y, x] <- 1 # assign 1 to indicate buffer
  }
  return(buffer.patches)
}
buffer.patches <- where.buffer(buffer.at) # initiate buffer

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
  nat.m <- array(0, c(nrow(SST), ncol(SST))) # create temp array
  m <- 1 - exp((-(SST[, ] - opt.temp)^2) / (temp.range^2)) # temperature based survival function from Norberg et al.
  m <- 1 - m # calculate mortality
  nat.m <- ifelse(m > s, s, m) # take natural mortality if it is less than temp based mortality
  return(nat.m)
}

############################################################################
## This function determines density dependent survival proportion for babies

survival_b <- function(num, SST) {
  s <- calc_temp_mortality(SST, opt.temp, temp.range, sb) # calculate temperature based mortality
  dd <- dd # density dependence of survival
  result <- s / (1 + dd * num)
  return(result)
}

############################################################################
## This function determines density dependent survival proportion for juveniles and adults

survival <- function(SST) {
  result <- calc_temp_mortality(SST, opt.temp, temp.range, s) # calculate temperature based mortality
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
  fished.num <- array(0, c(NUM.age.classes, NUM.sexes, NUM.genotypes)) # create fishing array
  if (gen <= pre.reserve.gens + pre.fishing.gens) { # determine if fishing takes place, but no MPA
    each.patch.pop <- array(0, c(NS.patches, EW.patches)) # create storage array
    for (i in 2:NUM.age.classes) {
      for (j in 1:NUM.sexes) {
        for (k in 1:NUM.genotypes) {
          each.patch.pop[, ] <- each.patch.pop[, ] + pop[, , i, j, k] # add population to storage array
        }
      }
    }
    mean.per.patch.pop <- mean(each.patch.pop) # calculate mean patch population
    ff <- mean.per.patch.pop * (1 / fished - 1) # distribute fishing pressure based on population in each patch
    patch.pop <- rowSums(pop[, , c(2, 3), , ], dims = 2) # get sum of rows in patch pop
    f <- patch.pop / (ff + patch.pop) # calculate patch fishing pressure
    for (i in 2:NUM.age.classes) {
      for (j in 1:NUM.sexes) {
        for (k in 1:NUM.genotypes) {
          survived <- Reshape(rbinom(NS.patches * EW.patches, pop[, , i, j, k], (1 - f)), NS.patches, EW.patches) # determine number of fish surviving fishing
          pop[, , i, j, k] <- survived
          fished.num[i, j, k] <- sum(each.patch.pop[, ]) - sum(survived) # save survivers in fished array
        }
      }
    }
  }
  if (gen > pre.reserve.gens + pre.fishing.gens) { # determine if there is an MPA
    reserve.area <- sum(reserve.patches[, , gen]) / (NS.patches * EW.patches) # calculate area of reserve
    buffer.area <- sum(buffer.patches) / (NS.patches * EW.patches) # calculate area of buffer region
    fished.adj <- (fished - (buffer.area * buffer.fished)) * 1 / (1 - (reserve.area + buffer.area)) # redistribute fishing pressure around reserve
    each.patch.pop <- array(0, c(NS.patches, EW.patches)) # create storage array
    for (i in 2:NUM.age.classes) {
      for (j in 1:NUM.sexes) {
        for (k in 1:NUM.genotypes) {
          each.patch.pop[, ] <- each.patch.pop[, ] + pop[, , i, j, k] # add population to storage array
        }
      }
    }
    each.patch.pop <- ifelse(reserve.patches[, , gen] == 1 | buffer.patches == 1, NaN, each.patch.pop) # remove population in reserve and buffer from average calculation
    mean.per.patch.pop <- mean(each.patch.pop, na.rm = TRUE) # calculate mean patch population
    ff <- mean.per.patch.pop * (1 / fished.adj - 1) # distribute fishing pressure based on population in each patch
    patch.pop <- rowSums(pop[, , c(2, 3), , ], dims = 2) # get sum of rows in patch pop
    patch.pop <- ifelse(reserve.patches[, , gen] == 1 | buffer.patches == 1, NaN, patch.pop) # remove fishing pressure from reserve and buffer
    f <- patch.pop / (ff + patch.pop) # calculate patch fishing pressure
    if (buffer.fished != 0) {
      f <- ifelse(buffer.patches == 1, buffer.fished, f) # if a buffer exists, remove fising from buffer
    }
    for (i in 2:NUM.age.classes) {
      for (j in 1:NUM.sexes) {
        for (k in 1:NUM.genotypes) {
          fished.array <- Reshape(rbinom(NS.patches * EW.patches, pop[, , i, j, k], (1 - f)), NS.patches, EW.patches) # determine number of fish surviving fishing
          pop[, , i, j, k] <- ifelse(is.na(fished.array[, ]), pop[, , i, j, k], fished.array[, ])
          fished.num[i, j, k] <- sum(each.patch.pop[, ], na.rm = TRUE) - sum(fished.array, na.rm = TRUE) # save survivers in fished array
        }
      }
    }
  }
  fished.list <- list("pop" = pop, "fish" = fished.num) # returned fish array with number of fish fished
  return(fished.list)
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
              
              #### Loop North/South
              # y_bool = ifelse((y <= patch.size*(NS.patches-lat)+patch.size/2) & (y >= patch.size/2-lat*patch.size), TRUE, FALSE)
              # y = ifelse((y_bool == FALSE) & (y > patch.size*(NS.patches-lat)+patch.size/2),(y - (patch.size * NS.patches)),y)
              # y = ifelse((y_bool == FALSE) & (y < patch.size/2-lat*patch.size),(y + (patch.size * NS.patches)),y)
              
              # convert movement distances into numbers of grid cells (assume fish start in centre of cell):
              x <- round(x) / patch.size # get number of patches
              y <- round(y) / patch.size # get number of patches
              xy <- as.data.frame(cbind(x, y)) # bind x and y together
              xy$count <- 1
              freq <- aggregate(count ~ x + y, data = xy, FUN = sum) # find number of fish that move each number of grids
              freq2D <- as.data.frame(array(0, c(length(unique(xy$y)), length(unique(xy$x)))))
              names(freq2D) <- sort(unique(xy$x))
              row.names(freq2D) <- sort(unique(xy$y)) # re frame data frame
              for (row in 1:length(freq$x)) {
                freq2D[as.character(freq$y[row]), as.character(freq$x[row])] <- freq$count[row] # assign frequencies to movement distance
              }
              # populate the move.array with movers (and stayers)
              for (xx in 1:length(unique(xy$x))) {
                for (yy in 1:length(unique(xy$y))) {
                  move.array[lat + as.numeric(row.names(freq2D)[yy]), lon + as.numeric(names(freq2D)[xx]), i, j, k] <- move.array[lat + as.numeric(row.names(freq2D)[yy]), lon + as.numeric(names(freq2D)[xx]), i, j, k] + freq2D[yy, xx] # move fish to new locations based on distances in movement dataframe
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
fished.array <- array(0, c(NUM.age.classes, NUM.sexes, NUM.genotypes, gens, reps))

start_time <- Sys.time()

SST.patches <- init_SST(years, "null") # assign either null, mean, enso, shock, or mean shock

for (rep in 1:reps) {
  print(rep)
  pop <- init()
  for (t in 1:gens) {
    output.array[, , , , , t, rep] <- pop
    pop <- spawn(pop)
    pop <- recruit(pop)
    if (t > pre.fishing.gens) {
      gen <- t
      fishing_result <- fishing(pop, gen)
      pop <- fishing_result$pop
      fished.array[, , , t, rep] <- fishing_result$fish
    }
    pop <- move(pop)
    print(t)
  }
  gc() # clear memory
}
gc()

end_time <- Sys.time()
end_time - start_time

# output array
save(output.array, file = here::here("model1.rda"))

# Output results into a dataframe
output_df <- data.frame() # create dataframe to hold results

# creates dataframe from array data
for (a in 1:reps) {
  for (b in 1:gens) {
    for (c in 1:NUM.genotypes) {
      for (d in 1:NUM.sexes) {
        for (e in 1:NUM.age.classes) {
          world_sub <- output.array[, , e, d, c, b, a] %>%
            as.data.frame()
          world_sub$rep <- paste0(a)
          world_sub$generation <- paste0(b)
          world_sub$genotype <- paste0(c)
          world_sub$sex <- paste0(d)
          world_sub$age <- paste0(e)
          world_sub$lat <- c(1:NS.patches)
          world_sub$max_temp <- SST.patches[1, 1, b]
          world_sub$mpa_temp <- SST.patches[50, 1, b]
          world_sub$min_temp <- SST.patches[100, 1, b]
          world_sub$fished <- fished.array[e, d, c, b, a]
          output_df <- bind_rows(output_df, world_sub)
        }
      }
    }
  }
}

# Wrangle dataframe into plottable format
output_df <- output_df %>%
  pivot_longer(V1:V20,
    names_to = "lon",
    values_to = "pop"
  ) %>%
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
  mutate(genotype = case_when( # assign real values to genotype
    genotype == 1 ~ "AA",
    genotype == 2 ~ "Aa",
    genotype == 3 ~ "aa"
  )) %>%
  mutate(genotype = as.factor(genotype)) %>% # turn genotype into factor
  mutate(sex = case_when( # assign real values to sex
    sex == 1 ~ "female",
    sex == 2 ~ "male"
  )) %>%
  mutate(sex = as.factor(sex)) %>% # turn sex into factor
  mutate(age = case_when( # assign real values to age
    age == 1 ~ "baby",
    age == 2 ~ "juvenile",
    age == 3 ~ "adult"
  )) %>%
  mutate(age = as.factor(age)) %>% # age as factor
  mutate(lat = as.numeric(lat)) %>%
  mutate(lon = as.numeric(lon)) %>%
  mutate(rep = as.numeric(rep)) %>%
  mutate(fished = as.numeric(fished)) %>%
  distinct()

# output CSV
write_csv(output_df, here::here("model1.csv"))

# Sum adults in population by location, generation, and genotype for each replicate
geno_sum2 <- output_df %>%
  filter(age == "adult") %>%
  group_by(lat, lon, generation, rep, genotype) %>%
  summarise(geno_pop_sum = sum(pop, na.rm = TRUE))

# average population across replicates
geno_mean2 <- geno_sum2 %>%
  group_by(lat, lon, generation, genotype) %>%
  summarise(
    geno_pop_mean = mean(geno_pop_sum, na.rm = TRUE),
    geno_pop_sd = sd(geno_pop_sum, na.rm = TRUE)
  ) %>%
  ungroup()

# sum adults in population by location and generation for each replicate
pop_sum2 <- output_df %>%
  filter(age == "adult") %>%
  group_by(lat, lon, generation, rep) %>%
  summarise(
    pop_sum = sum(pop, na.rm = TRUE),
    max_temp = max_temp, # get min temperature
    min_temp = min_temp, # get max temperature
    mpa_temp = mpa_temp, # get temperature within MPA
    fished_sum = sum(fished)
  )

# average population across replicates
pop_mean2 <- pop_sum2 %>%
  group_by(lat, lon, generation) %>%
  summarise(
    pop_mean = mean(pop_sum, na.rm = TRUE),
    pop_sd = sd(pop_sum, na.rm = TRUE),
    max_temp = max_temp,
    min_temp = min_temp,
    mpa_temp = mpa_temp,
    fished_mean = mean(fished_sum, na.rm = TRUE)
  ) %>%
  ungroup()

output_sum2 <- full_join(geno_mean2, pop_mean2) # join genotype based sum with full sum

output_sum2 <- output_sum2 %>%
  distinct() %>%
  mutate(
    freq = geno_pop_mean / pop_mean, # calculate frequency of each genotype
    mpa_size = "Large", # assign MPA size based on model_list.xlsx
    climate = "Mean Shock", # assign climate
    evolution = "Yes"
  ) # assign evolution or not

# output summary data to merge in 02_data_merge.Rmd
write_csv(output_sum2, here::here("summary_m1.csv"))
