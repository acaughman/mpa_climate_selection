library(pracma)


##move

move2 <- function(pop) {
  
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

pre.fishing.gens <- 5
pre.reserve.gens <- 10
post.reserve.gens <- 15
gens <- pre.fishing.gens+pre.reserve.gens+post.reserve.gens

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

for(rep in 1:reps) {
  #print(rep)
  pop2 <- init()
  SST.patches <- init_SST(gens)
  for(t in 1:gens) {
    #print(SST.patches[,,t])
    output.array[,,,,,t,rep] <- pop2
    pop <- spawn(pop2)
    pop <- recruit(pop2)
    if(t > pre.fishing.gens) {
      gen <- t
      pop2 <- fishing(pop,gen)
    }
    pop2 <- move2(pop2)
    print(t)
  }
}
# pop
# pop2
# pop - pop2
mean(pop - pop2)

