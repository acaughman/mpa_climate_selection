library(pracma)

##fishing


fishing2 <- function(pop,gen) {
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

##move

pop = init()
pop2 = init()
SST.patches = init_SST(years = 84)
pop = spawn(pop)
pop2 = spawn(pop2)
pop = recruit(pop) 
pop2 = recruit(pop2)
pop = fishing(pop,20) 
pop2 = fishing2(pop2, 20)
pop
pop2
pop - pop2
mean(pop - pop2)
