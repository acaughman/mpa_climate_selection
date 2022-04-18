# Set Up ------------------------------------------------------------------

library(here)
library(tidyverse)
library(patchwork)
library(beepr)

set.seed(42)

# Fish Class --------------------------------------------------------------

Fish = setRefClass("Fish", fields = list(sex = "numeric", genotype = "numeric", age_class = "numeric", num = "numeric"))

# Parameters --------------------------------------------------------------

NUM.reps <- 1 # The number of replicate simulations to run
## 150 years total
NUM.gens.pre.fishing <- 25 # The number of generations before any fishery
NUM.gens.pre.reserve <- 50 # The number of generations of fishing before reserves are installed
NUM.gens.post.reserve <- 75 # The number of generations with the reserve installed
years = NUM.gens.pre.fishing+NUM.gens.pre.reserve+NUM.gens.post.reserve

NS.patches <- 5 # the number of patches on the north-south axis
EW.patches <- 5 # the number of patches on the east-west axis
patch.size <- 100 # the width and height of each grid cell in nautical miles (COULD BE METERS?)
## View the "world" coordinates:
view.world <- array(seq(1,NS.patches*EW.patches),c(NS.patches,EW.patches))
view.world

init.a <- 0.3 # The initial frequency of the low movement allele

sb <- 0.37 # survival proportion for babies
s <- 0.37 # survival proportion
dd <- 0.0005 # density dependence of baby survival 
fecundity <- 1500 # The number of babies produced, on average, by each adult female each year.
maturity.age <- 1.5 # The average age at which individuals mature (i.e., the age at which 50% of individuals are mature)
fished.factor <- 0.7
#fished <- fished.factor*(1-s) # Fishing mortalty: the proportion of adults that get fished per year
fished <- fished.factor
buffer.fished <- 0 #buffer fishing pressure (lower than total = buffer zone, higher than total = fishing the line)
reserves.at <- c(12,13,17,18) # This determines which patches are marine reserves. Should be a list: e.g., for one reserve, c(369,370,371,372,389,390,391,392,409,410,411,412,429,430,431,432)
buffer.at <- c()
bold.mover.distance <- 200 # Individuals with AA genotype move this distance on average every year, in nautical miles
lazy.mover.distance <- 100 # Individuals with aa genotype move this distance on average every year, in nautical miles
Dominance.coefficient <- 0.5 # Dominance coefficient
Heritability.index <- 2 # Influences stochastic variation in movement distance. High numbers decrease variation by reducing the variance around the phenotypic mean in a negative binomial distribution. The phenotypic mean is determined by the genotype.
opt.temp = 25 #optimal temperature of species
temp.range = 5 #thermal breath of species

NUM.age.classes <- 3 #babies, juvenile, adult
NUM.sexes <- 2 #female male
NUM.genotypes <- 3 #AA,Aa,aa

world = matrix(list(), nrow = NS.patches, ncol = EW.patches)

# Reserve and Buffer Function ---------------------------------------------

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

# Init Function -----------------------------------------------------------

## world[[1]][[1]]$genotype to access peices

init <- function() {
  init.AA <- list(Fish(sex = 1, age_class = 1, genotype = 1, num = round(200*(1-init.a)^2)), 
                  Fish(sex = 1, age_class = 2, genotype = 1, num = round(200*(1-init.a)^2)),
                  Fish(sex = 1, age_class = 3, genotype = 1, num = round(200*(1-init.a)^2)),
                  Fish(sex = 2, age_class = 1, genotype = 1, num = round(200*(1-init.a)^2)),
                  Fish(sex = 2, age_class = 2, genotype = 1, num = round(200*(1-init.a)^2)),
                  Fish(sex = 2, age_class = 3, genotype = 1, num = round(200*(1-init.a)^2)))     
  init.Aa <- list(Fish(sex = 1, age_class = 1, genotype = 2, num = round(200*2*(init.a)*(1-init.a))), 
                  Fish(sex = 1, age_class = 2, genotype = 2, num = round(200*2*(init.a)*(1-init.a))),
                  Fish(sex = 1, age_class = 3, genotype = 2, num = round(200*2*(init.a)*(1-init.a))),
                  Fish(sex = 2, age_class = 1, genotype = 2, num = round(200*2*(init.a)*(1-init.a))),
                  Fish(sex = 2, age_class = 2, genotype = 2, num = round(200*2*(init.a)*(1-init.a))),
                  Fish(sex = 2, age_class = 3, genotype = 2, num = round(200*2*(init.a)*(1-init.a))))                         
  init.aa <- list(Fish(sex = 1, age_class = 1, genotype = 3, num = round(200*(init.a)^2)),
                  Fish(sex = 1, age_class = 2, genotype = 3, num = round(200*(init.a)^2)),
                  Fish(sex = 1, age_class = 3, genotype = 3, num = round(200*(init.a)^2)),
                  Fish(sex = 2, age_class = 1, genotype = 3, num = round(200*(init.a)^2)),
                  Fish(sex = 2, age_class = 2, genotype = 3, num = round(200*(init.a)^2)),
                  Fish(sex = 2, age_class = 3, genotype = 3, num = round(200*(init.a)^2)))                        
  pop <- world
  for(lat in 1:NS.patches) {
    for(lon in 1:EW.patches) {
      pop[lat,lon] <- list(init.AA)
      pop[lat,lon] <- list(init.Aa)
      pop[lat,lon] <- list(init.aa)
    }
  }
  return(pop)
}


# Sea Surface Temperature Init --------------------------------------------

init_SST <- function(years) {
  
  ### UNCOMMENT FOR CONSISTENT SST
  SST.patches <- array(25, c(NS.patches, EW.patches, years))
  
  ### UNCOMMENT FOR CONSTANT MEAN SHIFT SST
  # SST.patches <- array(0, c(NS.patches, EW.patches, years))
  # start_SST = opt.temp + NS.patches*0.1
  # 
  # for (i in 1:years) {
  #   SST = start_SST
  #   if (SST > 35) {
  #     SST = 35
  #   }
  #   for (lat in 1:NS.patches) {
  #     SST.patches[lat,,i] = SST
  #     SST = SST - 0.1
  #   }
  #   start_SST = start_SST + 0.018
  # }
  
  
  ### UNCOMMENT FOR LOW VARIABLE MEAN SST
  # SST.patches <- array(0, c(NS.patches, EW.patches, years))
  # start_SST = opt.temp + NS.patches*0.1
  # 
  # for (i in 1:years) {
  #   SST = start_SST
  #   if (SST > 35) {
  #     SST = 35
  #   }
  #   for (lat in 1:NS.patches) {
  #     SST.patches[lat,,i] = SST
  #     SST = SST - 0.1
  #   }
  #   start_SST = start_SST + rnorm(1, mean = 0.018, sd = 0.01)
  # }
  
  
  ### UNCOMMENT FOR ENSO VARIABLE MEAN SST
  # SST.patches <- array(0, c(NS.patches, EW.patches, years))
  # start_SST = opt.temp + NS.patches*0.1
  # 
  # for (i in 1:years) {
  #   SST = start_SST
  #   if (SST > 35) {
  #     SST = 35
  #   }
  #   for (lat in 1:NS.patches) {
  #     SST.patches[lat,,i] = SST
  #     SST = SST - 0.1
  #   }
  #   start_SST = start_SST + rnorm(1, mean = 0.018, sd = 0.1)
  # }
  # #
  # 
  return(SST.patches) ### DO NOT COMMENT OUT
}


# Spawn -------------------------------------------------------------------

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





# The Simulation ----------------------------------------------------------
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
