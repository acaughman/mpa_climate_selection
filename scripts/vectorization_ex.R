world <- array(25, c(8, 8))
world2 <- array(25, c(8, 8))

for (i in 1:8) {
  for (j in 1:8) {
    world[i,j] = rbinom(1,world[i,j],(1-f))
  }
}

world2[,] = rbinom(8*8,world2[,],(1-f))
