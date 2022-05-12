library(here)
library(tidyverse)
library(patchwork)
library(pracma)
library(beepr)
library(reticulate)

seed = 42
addTaskCallback(function(...) {set.seed(seed);TRUE})

############################################################################
## BUFFER CODE PLAY ##########################################################

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
    patch.pop <- rowSums(pop[,,c(2,3),,], dims = 2)
    f <- patch.pop/(ff+patch.pop)
    for(i in 2:NUM.age.classes) {
      for(j in 1:NUM.sexes) {
        for(k in 1:NUM.genotypes) {
          pop[,,i,j,k] <- Reshape(rbinom(NS.patches * EW.patches,pop[,,i,j,k],(1-f)), NS.patches, EW.patches)
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
    each.patch.pop = ifelse(reserve.patches == 1, NaN, each.patch.pop)
    mean.per.patch.pop <- mean(each.patch.pop,na.rm=TRUE)
    ff <- mean.per.patch.pop*(1/fished.adj-1)
    patch.pop <- rowSums(pop[,,c(2,3),,], dims=2)
    patch.pop = ifelse(reserve.patches == 1, NaN, patch.pop)
    f <- patch.pop/(ff+patch.pop)
    for(i in 2:NUM.age.classes) {
      for(j in 1:NUM.sexes) {
        for(k in 1:NUM.genotypes) {
          fished.array = Reshape(rbinom(NS.patches * EW.patches,pop[,,i,j,k],(1-f)), NS.patches, EW.patches)
          pop[,,i,j,k] <- ifelse(is.na(fished.array[,]),pop[,,i,j,k],fished.array[,])
        }
      }
    }
  }
  return(pop)
}

############################################################################
## THE SIMULATION ##########################################################

reps <- NUM.reps

pre.fishing.gens <- NUM.gens.pre.fishing
pre.reserve.gens <- NUM.gens.pre.reserve
post.reserve.gens <- NUM.gens.post.reserve
gens <- pre.fishing.gens+pre.reserve.gens+post.reserve.gens

output.array <- array(0 ,c(NS.patches, EW.patches, NUM.age.classes, NUM.sexes, NUM.genotypes, gens, reps))

## NEW SIM
set.seed(seed)
start_time <- Sys.time()

for(rep in 1:reps) {
  #print(rep)
  pop <- init()
  SST.patches <- init_SST(gens)
  for(t in 1:gens) {
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
  gc() #clear memory
}
gc()

end_time <- Sys.time()
end_time - start_time # 30.24764


### OLD SIM
set.seed(seed)
start_time <- Sys.time()

for(rep in 1:reps) {
  #print(rep)
  pop2 <- init()
  SST.patches <- init_SST(gens)
  for(t in 1:gens) {
    output.array[,,,,,t,rep] <- pop2
    pop2 <- spawn2(pop2)
    pop2 <- recruit2(pop2)
    if(t > pre.fishing.gens) {
      gen <- t
      pop2 <- fishing2(pop2,gen)
    }
    pop2 <- move(pop2)
    print(t)
  }
  gc() #clear memory
}
gc()

end_time <- Sys.time()
end_time - start_time #41.98303 secs

set.seed(seed)
pop = init()
pop = spawn(pop)
pop = recruit(pop)
pop = fishing(pop, 83)

set.seed(seed)
pop2 = init()
pop2 = spawn2(pop2)
pop2 = recruit2(pop2)
pop2 = fishing2(pop2, 83)


mean(pop - pop2)
