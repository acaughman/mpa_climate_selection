library(tidyverse)

init <- function() {
  init.AA <- round(50*(1-.3)^2)
  init.Aa <- round(50*2*(.3)*(1-.3))
  init.aa <- round(50*(.3)^2)
  pop <- world
  pop[,,,,1] <- init.AA
  pop[,,,,2] <- init.Aa
  pop[,,,,3] <- init.aa
  return(pop)
}

patch.size = 100
EW.patches = 10
NS.patches = 20

world = array(0, c(NS.patches, EW.patches, 3, 2, 3))

pop = init()

h <- .5 # Dominance coefficient
Hom.A.movers <- 400 # Individuals with AA genotype move this distance on average, in nautical miles
Hom.a.movers <- 300 # Individuals with aa genotype move this distance on average, in nautical miles
Het.movers <- min(Hom.A.movers,Hom.a.movers) + h * abs(Hom.A.movers - Hom.a.movers)
Herit <- 2 # Influences heritability of movement. High numbers increase heritability by reducing the variance around the phenotypic mean. The phenotypic mean is determined by the genotype.

move.array <- array(0, c(NS.patches, EW.patches, 3, 2, 3))
move.array2 <- array(0, c(NS.patches, EW.patches, 3, 2, 3))

for(lat in 5:6) {
  for(lon in 1:2) {
    for(i in 3:3) {
      for(j in 2:2) {
        for(k in 1:1) {
          if(k == 1) { mean.dist <- Hom.A.movers }
          if(k == 2) { mean.dist <- Het.movers }
          if(k == 3) { mean.dist <- Hom.a.movers }
          if(pop[lat,lon,i,j,k] > 0) {
            # movers are subtracted from the present grid cell
            move.array[lat,lon,i,j,k] <- move.array[lat,lon,i,j,k] - pop[lat,lon,i,j,k]
            move.array2[lat,lon,i,j,k] <- move.array2[lat,lon,i,j,k] - pop[lat,lon,i,j,k]
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
            x.round = round(x/patch.size)
            y.round = round(y/patch.size)
            xy = as.data.frame(cbind(x.round,y.round))
            xy_sum = xy %>%
              group_by(x.round, y.round) %>%
              summarise(sum = n())
            xy_pivot = xy_sum %>%
              ungroup() %>%
              pivot_wider(values_from = sum, names_from = x.round)
            final_xy = xy_pivot %>%
              select(-y.round)
            final_xy[is.na(final_xy)] <- 0
            final_xy = as.data.frame(final_xy)
            rownames(final_xy) = sort(unique(xy_sum$y.round))
            # populate the move.array with movers (and stayers)
            for(xx in 1:length(unique(xy_sum$x.round))) {
              for(yy in 1:length(unique(xy_sum$y.round))) {
                move.array[lat+as.numeric(row.names(final_xy)[yy]),lon+as.numeric(names(final_xy)[xx]),i,j,k] <- move.array[lat+as.numeric(row.names(final_xy)[yy]),lon+as.numeric(names(final_xy)[xx]),i,j,k] + final_xy[yy,xx]
              }
            }
            # convert movement distances into numbers of grid cells (assume fish start in centre of cell):
            xy2 <- as.data.frame(cbind(x,y))
            hx <- hist(xy2$x, breaks = seq(round(min(x),-2)-patch.size/2,round(max(x),-2)+patch.size/2,patch.size), plot = FALSE)
            xbreaks <- hx$breaks
            xmids <- hx$mids
            hy <- hist(xy2$y, breaks = seq(round(min(y),-2)-patch.size/2,round(max(y),-2)+patch.size/2,patch.size), plot = FALSE)
            ybreaks <- hy$breaks
            ymids <- hy$mids
            freq <-  as.data.frame(table(cut(xy2[,2], ybreaks,labels=ymids),cut(xy2[,1], xbreaks,labels=xmids)))
            freq2D <- as.data.frame(array(0,c(length(ymids),length(xmids))))
            names(freq2D) <- xmids/patch.size
            row.names(freq2D) <- ymids/patch.size
            freq2D[cbind(freq[,1], freq[,2])] <- freq[,3]
            # populate the move.array with movers (and stayers)
            for(xx in 1:length(xmids)) {
              for(yy in 1:length(ymids)) {
                move.array2[lat+as.numeric(row.names(freq2D)[yy]),lon+as.numeric(names(freq2D)[xx]),i,j,k] <- move.array2[lat+as.numeric(row.names(freq2D)[yy]),lon+as.numeric(names(freq2D)[xx]),i,j,k] + freq2D[yy,xx]
              }
            }
          }
        }
      }
    }
  }
}
