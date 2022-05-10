library(here)
library(tidyverse)
library(patchwork)
library(pracma)
library(beepr)
library(reticulate)

seed = 42
addTaskCallback(function(...) {set.seed(seed);TRUE})



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
