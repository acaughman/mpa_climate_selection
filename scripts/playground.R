# setup -------------------------------------------------------------------

#run final model parameters and functions first

addTaskCallback(function(...) {set.seed(42);TRUE})
options(warn=-1)
options(dplyr.summarise.inform = FALSE)

reps = 1

pre.fishing.gens <- 1
pre.reserve.gens <- 2
post.reserve.gens <- 3
gens <- 5

output.array <- array(0 ,c(NS.patches, EW.patches, NUM.age.classes, NUM.sexes, NUM.genotypes, gens, reps))
output.array.ifelse <- array(0 ,c(NS.patches, EW.patches, NUM.age.classes, NUM.sexes, NUM.genotypes, gens, reps))

# current work ------------------------------------------------------------

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
                    y[m] <- y[m] - (patch.size * NS.patches)
                    # distance penalty for hitting an edge
                    #y[m] <- y[m] + 1
                  }
                  if(y[m] < patch.size/2-lat*patch.size) {
                    y[m] <- y[m] + (patch.size * NS.patches)
                    # distance penalty for hitting an edge
                    #y[m] <- y[m] - 1
                  }
                }
              }
              # convert movement distances into numbers of grid cells (assume fish start in centre of cell):
              xy <- as.data.frame(cbind(round(x),round(y))) %>% 
                rename(x = V1) %>% 
                rename(y = V2)
              freq = xy %>% 
                group_by(x,y) %>% 
                summarize(count = n())
              freq2D = as.data.frame(array(0,c(length(unique(xy$y)), length(unique(xy$x)))))
              names(freq2D) <- sort(unique(xy$x))
              row.names(freq2D) <- sort(unique(xy$y))
              for(row in 1:length(freq$x)) {
                freq2D[as.character(freq$y[row]), as.character(freq$x[row])] = freq$count[row]
              }
              # populate the move.array with movers (and stayers)
              for(xx in 1:length(unique(xy$x))) {
                for(yy in 1:length(unique(xy$y))) {
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

move_ifelse <- function(pop) {
  
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
              dist <- ifelse(dist > EW.patches, EW.patches, dist)
              # determine the direction of each move
              theta <- runif(pop[lat,lon,i,j,k],0,2*pi)
              # bias this movement in the north-south direction (along coasts) if this is a great white shark simulation (otherwise, comment out the next three lines):
              f.adj <- function(x, u) x-cos(x)*sin(x) - u
              my.uniroot <- function(x) uniroot(f.adj, c(0, 2*pi), tol = 0.0001, u = x)$root
              theta <- vapply(theta, my.uniroot, numeric(1))
              # convert direction and distance into a distance in the x-direction (longitude)
              x <- cos(theta)*dist
              #bounce off edges
              x_bool = ifelse((x <= patch.size*(EW.patches-lon)+patch.size/2) & (x >= patch.size/2-lon*patch.size), TRUE, FALSE)
              x = ifelse((x_bool == FALSE) & (x > patch.size*(EW.patches-lon)+patch.size/2),(-(x-2*(patch.size*(EW.patches-lon)+patch.size/2))),x)
              x = ifelse((x_bool == FALSE) & (x < patch.size/2-lon*patch.size),(-(x-2*(patch.size/2-lon*patch.size))),x)
              # convert direction and distance into a distance in the y-direction (latitude)
              y <- sin(theta)*dist
              #bounce off edges
              y_bool = ifelse((y <= patch.size*(NS.patches-lat)+patch.size/2) & (y >= patch.size/2-lat*patch.size), TRUE, FALSE)
              y = ifelse((y_bool == FALSE) & (y > patch.size*(NS.patches-lat)+patch.size/2),(y - (patch.size * NS.patches)),y)
              y = ifelse((y_bool == FALSE) & (y < patch.size/2-lat*patch.size),(y + (patch.size * NS.patches)),y)
              # convert movement distances into numbers of grid cells (assume fish start in centre of cell):
              x = round(x)
              y = round(y)
              xy <- as.data.frame(cbind(x,y))
              freq <<- xy %>% 
                group_by(x,y) %>% 
                summarize(count = n())
              freq2D = as.data.frame(array(0,c(length(unique(xy$y)), length(unique(xy$x)))))
              names(freq2D) <- sort(unique(xy$x))
              row.names(freq2D) <- sort(unique(xy$y))
              for(row in 1:length(freq$x)) {
                freq2D[as.character(freq$y[row]), as.character(freq$x[row])] = freq$count[row]
              }
              # populate the move.array with movers (and stayers)
              for(xx in 1:length(unique(xy$x))) {
                for(yy in 1:length(unique(xy$y))) {
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
  gc()
}

# current simulation ------------------------------------------------------

start_time <- Sys.time()

for(rep in 1:reps) {
  print(rep)
  pop <- init()
  #save(SST.patches, file = here::here("data", "null.rda"))
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


# working simulation ------------------------------------------------------

start_time_ifelse <- Sys.time()

for(rep in 1:reps) {
  print(rep)
  pop <- init()
  #save(SST.patches, file = here::here("data", "null.rda"))
  for(t in 1:gens) {
    output.array[,,,,,t,rep] <- pop
    pop <- spawn(pop)
    pop <- recruit(pop)
    if(t > pre.fishing.gens) {
      gen <- t
      pop <- fishing(pop,gen)
    }
    pop <- move_ifelse(pop)
    print(t)
  }
  gc() #clear memory
}
gc()

end_time_ifelse <- Sys.time()

# sim comparison ----------------------------------------------------------

end_time - start_time
end_time_ifelse - start_time_ifelse

# Old Move 7/23/2022 ------------------------------------------------------

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
                    y[m] <- y[m] - (patch.size * NS.patches)
                    # distance penalty for hitting an edge
                    #y[m] <- y[m] + 1
                  }
                  if(y[m] < patch.size/2-lat*patch.size) {
                    y[m] <- y[m] + (patch.size * NS.patches)
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

