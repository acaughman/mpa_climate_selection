library(pracma)

##recruit

calc_temp_mortality2 <- function(SST, opt.temp, temp.range, s) {
  nat.m = array(0, c(nrow(SST), ncol(SST)))
  m = 1 - exp((-(SST[,] - opt.temp)^2)/(temp.range^2)) # temperature based mortality function from Walsworth et al.
  m = 1 - m
  for (i in 1:nrow(m)) {
    for (j in 1:ncol(m)) {
      if(m[i,j] > s) {
        nat.m[i,j] = s
      } else if (m[i,j] < s) {
        nat.m[i,j] = m[i,j]
      }
    }
  }
  return(nat.m)
}

survival_b2 <- function(num, SST) {
  s = calc_temp_mortality2(SST, opt.temp, temp.range, sb)
  dd <- dd # density dependence of survival
  result <- s/(1 + dd * num)
  return(result)
}

survival2 <- function(SST) {
  result = calc_temp_mortality2(SST, opt.temp, temp.range, s)
  return(result)
}

p <- 1/(maturity.age)

recruit2 <- function(pop) {
  recruit.array <- world
  SST = SST.patches[,, t]
  for(i in 1:NUM.age.classes) {
    if(i == 1) {
      # Some babies survive and recruit to juvenile age class
      s1 <- survival_b2(sum(pop[,,i,,]), SST)
      s <- survival2(SST)
      for(j in 1:NUM.sexes) {
        for(k in 1:NUM.genotypes) {
          recruit.array[,,i+1,j,k] <- recruit.array[,,i+1,j,k] + Reshape(rbinom(NS.patches * EW.patches,pop[,,i,j,k],s1), NS.patches, EW.patches)
        }
      }
    } 
    if(i == 2) {
      # Some juveniles survive
      for(j in 1:NUM.sexes) {
        for(k in 1:NUM.genotypes) {
          juvies.surviving <- rbinom(NS.patches * EW.patches,pop[,,i,j,k],s)
          # Some juveniles recruit to adult age class
          juvies.recruiting <- rbinom(NS.patches * EW.patches,juvies.surviving,p)
          recruit.array[,,i+1,j,k] <- recruit.array[,,i+1,j,k] + juvies.recruiting
          # The rest of the juveniles remain in the juvenile age class
          recruit.array[,,i,j,k] <- recruit.array[,,i,j,k] + juvies.surviving-juvies.recruiting
        }
      }
    }
    if(i == 3) {
      # Some adults survive
      for(j in 1:NUM.sexes) {
        for(k in 1:NUM.genotypes) {
          recruit.array[,,i,j,k] <- recruit.array[,,i,j,k] + Reshape(rbinom(NS.patches * EW.patches,pop[l,,i,j,k],s), NS.patches, EW.patches)
        }
      }
    }
  }
  return(recruit.array)
}

##fishing


##move

pop = init()
pop2 = init()
SST.patches = init_SST(years = 84)
pop = spawn(pop)
pop2 = spawn(pop2)
pop = recruit(pop)
pop2 = recruit2(pop2)
pop
pop2
