# setup -------------------------------------------------------------------

# run final model parameters and functions first

addTaskCallback(function(...) {
  set.seed(42)
  TRUE
})
options(warn = -1)
options(dplyr.summarise.inform = FALSE)

reps <- 1

pre.fishing.gens <- 1
pre.reserve.gens <- 2
post.reserve.gens <- 3
gens <- 5

output.array <- array(0, c(NS.patches, EW.patches, NUM.age.classes, NUM.sexes, NUM.genotypes, gens, reps))
output.array2 <- array(0, c(NS.patches, EW.patches, NUM.age.classes, NUM.sexes, NUM.genotypes, gens, reps))

# current simulation ------------------------------------------------------

start_time <- Sys.time()

for (rep in 1:reps) {
  print(rep)
  pop <- init()
  # save(SST.patches, file = here::here("data", "null.rda"))
  for (t in 1:gens) {
    output.array[, , , , , t, rep] <- pop
    pop <- spawn(pop)
    pop <- recruit(pop)
    if (t > pre.fishing.gens) {
      gen <- t
      pop <- fishing(pop, gen)
    }
    pop <- move(pop)
    print(t)
  }
  gc() # clear memory
}
gc()

end_time <- Sys.time()


# working simulation ------------------------------------------------------

start_time2 <- Sys.time()

for (rep in 1:reps) {
  print(rep)
  pop <- init()
  # save(SST.patches, file = here::here("data", "null.rda"))
  for (t in 1:gens) {
    output.array2[, , , , , t, rep] <- pop
    pop <- spawn(pop)
    pop <- recruit(pop)
    if (t > pre.fishing.gens) {
      gen <- t
      pop <- fishing(pop, gen)
    }
    pop <- move2(pop)
    print(t)
  }
  gc() # clear memory
}
gc()

end_time2 <- Sys.time()

# sim comparison ----------------------------------------------------------

end_time - start_time
end_time2 - start_time2

# current work ------------------------------------------------------------



# replace tidyverse -------------------------------------------------------

move2 <- function(pop) {
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
            if (pop[lat, lon, i, j, k] > 0) {
              # movers are subtracted from the present grid cell
              move.array[lat, lon, i, j, k] <- move.array[lat, lon, i, j, k] - pop[lat, lon, i, j, k]
              # determine the distribution of movement distances in nautical miles:
              dist <- rnbinom(pop[lat, lon, i, j, k], mu = mean.dist, size = Herit)
              dist <- ifelse(dist > EW.patches, EW.patches, dist)
              # determine the direction of each move
              theta <- runif(pop[lat, lon, i, j, k], 0, 2 * pi)
              # bias this movement in the north-south direction (along coasts) if this is a great white shark simulation (otherwise, comment out the next three lines):
              f.adj <- function(x, u) x - cos(x) * sin(x) - u
              my.uniroot <- function(x) uniroot(f.adj, c(0, 2 * pi), tol = 0.0001, u = x)$root
              theta <- vapply(theta, my.uniroot, numeric(1))
              # convert direction and distance into a distance in the x-direction (longitude)
              x <- cos(theta) * dist
              # bounce off edges
              x_bool <- ifelse((x <= patch.size * (EW.patches - lon) + patch.size / 2) & (x >= patch.size / 2 - lon * patch.size), TRUE, FALSE)
              x <- ifelse((x_bool == FALSE) & (x > patch.size * (EW.patches - lon) + patch.size / 2), (-(x - 2 * (patch.size * (EW.patches - lon) + patch.size / 2))), x)
              x <- ifelse((x_bool == FALSE) & (x < patch.size / 2 - lon * patch.size), (-(x - 2 * (patch.size / 2 - lon * patch.size))), x)
              # convert direction and distance into a distance in the y-direction (latitude)
              y <- sin(theta) * dist
              # bounce off edges
              y_bool <- ifelse((y <= patch.size * (NS.patches - lat) + patch.size / 2) & (y >= patch.size / 2 - lat * patch.size), TRUE, FALSE)
              y <- ifelse((y_bool == FALSE) & (y > patch.size * (NS.patches - lat) + patch.size / 2), (y - (patch.size * NS.patches)), y)
              y <- ifelse((y_bool == FALSE) & (y < patch.size / 2 - lat * patch.size), (y + (patch.size * NS.patches)), y)
              # convert movement distances into numbers of grid cells (assume fish start in centre of cell):
              x <- round(x)
              y <- round(y)
              xy <- as.data.frame(cbind(x, y))
              xy$count <- 1
              freq <- aggregate(count ~ x + y, data = xy, FUN = sum)
              freq2D <- as.data.frame(array(0, c(length(unique(xy$y)), length(unique(xy$x)))))
              names(freq2D) <- sort(unique(xy$x))
              row.names(freq2D) <- sort(unique(xy$y))
              for (row in 1:length(freq$x)) {
                freq2D[as.character(freq$y[row]), as.character(freq$x[row])] <- freq$count[row]
              }
              # populate the move.array with movers (and stayers)
              for (xx in 1:length(unique(xy$x))) {
                for (yy in 1:length(unique(xy$y))) {
                  move.array[lat + as.numeric(row.names(freq2D)[yy]), lon + as.numeric(names(freq2D)[xx]), i, j, k] <- move.array[lat + as.numeric(row.names(freq2D)[yy]), lon + as.numeric(names(freq2D)[xx]), i, j, k] + freq2D[yy, xx]
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
            if (pop[lat, lon, i, j, k] > 0) {
              # movers are subtracted from the present grid cell
              move.array[lat, lon, i, j, k] <- move.array[lat, lon, i, j, k] - pop[lat, lon, i, j, k]
              # determine the distribution of movement distances in nautical miles:
              dist <- rnbinom(pop[lat, lon, i, j, k], mu = mean.dist, size = Herit)
              dist <- ifelse(dist > EW.patches, EW.patches, dist)
              # determine the direction of each move
              theta <- runif(pop[lat, lon, i, j, k], 0, 2 * pi)
              # bias this movement in the north-south direction (along coasts) if this is a great white shark simulation (otherwise, comment out the next three lines):
              f.adj <- function(x, u) x - cos(x) * sin(x) - u
              my.uniroot <- function(x) uniroot(f.adj, c(0, 2 * pi), tol = 0.0001, u = x)$root
              theta <- vapply(theta, my.uniroot, numeric(1))
              # convert direction and distance into a distance in the x-direction (longitude)
              x <- cos(theta) * dist
              # bounce off edges
              x_bool <- ifelse((x <= patch.size * (EW.patches - lon) + patch.size / 2) & (x >= patch.size / 2 - lon * patch.size), TRUE, FALSE)
              x <- ifelse((x_bool == FALSE) & (x > patch.size * (EW.patches - lon) + patch.size / 2), (-(x - 2 * (patch.size * (EW.patches - lon) + patch.size / 2))), x)
              x <- ifelse((x_bool == FALSE) & (x < patch.size / 2 - lon * patch.size), (-(x - 2 * (patch.size / 2 - lon * patch.size))), x)
              # convert direction and distance into a distance in the y-direction (latitude)
              y <- sin(theta) * dist
              # bounce off edges
              y_bool <- ifelse((y <= patch.size * (NS.patches - lat) + patch.size / 2) & (y >= patch.size / 2 - lat * patch.size), TRUE, FALSE)
              y <- ifelse((y_bool == FALSE) & (y > patch.size * (NS.patches - lat) + patch.size / 2), (y - (patch.size * NS.patches)), y)
              y <- ifelse((y_bool == FALSE) & (y < patch.size / 2 - lat * patch.size), (y + (patch.size * NS.patches)), y)
              # convert movement distances into numbers of grid cells (assume fish start in centre of cell):
              x <- round(x)
              y <- round(y)
              xy <- as.data.frame(cbind(x, y))
              freq <- xy %>%
                group_by(x, y) %>%
                summarize(count = n())
              freq2D <- as.data.frame(array(0, c(length(unique(xy$y)), length(unique(xy$x)))))
              names(freq2D) <- sort(unique(xy$x))
              row.names(freq2D) <- sort(unique(xy$y))
              for (row in 1:length(freq$x)) {
                freq2D[as.character(freq$y[row]), as.character(freq$x[row])] <- freq$count[row]
              }
              # populate the move.array with movers (and stayers)
              for (xx in 1:length(unique(xy$x))) {
                for (yy in 1:length(unique(xy$y))) {
                  move.array[lat + as.numeric(row.names(freq2D)[yy]), lon + as.numeric(names(freq2D)[xx]), i, j, k] <- move.array[lat + as.numeric(row.names(freq2D)[yy]), lon + as.numeric(names(freq2D)[xx]), i, j, k] + freq2D[yy, xx]
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

# profvis -----------------------------------------------------------------

profvis::profvis({
  gens <- 1

  for (rep in 1:reps) {
    print(rep)
    pop <- init()
    # save(SST.patches, file = here::here("data", "enso.rda"))
    for (t in 1:gens) {
      output.array[, , , , , t, rep] <- pop
      pop <- spawn(pop)
      pop <- recruit(pop)
      if (t > pre.fishing.gens) {
        gen <- t
        pop <- fishing(pop, gen)
      }
      pop <- move(pop)
      print(t)
    }
    gc() # clear memory
  }
})

# unvectorized move -------------------------------------------------------

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
            if (pop[lat, lon, i, j, k] > 0) {
              # movers are subtracted from the present grid cell
              move.array[lat, lon, i, j, k] <- move.array[lat, lon, i, j, k] - pop[lat, lon, i, j, k]
              # determine the distribution of movement distances in nautical miles:
              dist <- rnbinom(pop[lat, lon, i, j, k], mu = mean.dist, size = Herit)
              # determine the direction of each move
              theta <- runif(pop[lat, lon, i, j, k], 0, 2 * pi)
              # bias this movement in the north-south direction (along coasts) if this is a great white shark simulation (otherwise, comment out the next three lines):
              f.adj <- function(x, u) x - cos(x) * sin(x) - u
              my.uniroot <- function(x) uniroot(f.adj, c(0, 2 * pi), tol = 0.0001, u = x)$root
              theta <- vapply(theta, my.uniroot, numeric(1))
              # convert direction and distance into a distance in the x-direction (longitude)
              x <- cos(theta) * dist
              # bounce off edges (assume fish start in centre of cell)
              for (m in 1:length(x)) {
                miss_x.edges <- FALSE
                while (miss_x.edges == FALSE) {
                  if (x[m] <= patch.size * (EW.patches - lon) + patch.size / 2) {
                    if (x[m] >= patch.size / 2 - lon * patch.size) {
                      miss_x.edges <- TRUE
                    }
                  }
                  if (x[m] > patch.size * (EW.patches - lon) + patch.size / 2) {
                    x[m] <- -(x[m] - 2 * (patch.size * (EW.patches - lon) + patch.size / 2))
                    # distance penalty for hitting an edge
                    # x[m] <- x[m] + 1
                  }
                  if (x[m] < patch.size / 2 - lon * patch.size) {
                    x[m] <- -(x[m] - 2 * (patch.size / 2 - lon * patch.size))
                    # distance penalty for hitting an edge
                    # x[m] <- x[m] - 1
                  }
                }
              }
              # convert direction and distance into a distance in the y-direction (latitude)
              y <- sin(theta) * dist
              # bounce off edges (assume fish start in centre of cell)
              for (m in 1:length(y)) {
                miss_y.edges <- FALSE
                while (miss_y.edges == FALSE) {
                  if (y[m] <= patch.size * (NS.patches - lat) + patch.size / 2) {
                    if (y[m] >= patch.size / 2 - lat * patch.size) {
                      miss_y.edges <- TRUE
                    }
                  }
                  if (y[m] > patch.size * (NS.patches - lat) + patch.size / 2) {
                    y[m] <- y[m] - (patch.size * NS.patches)
                    # distance penalty for hitting an edge
                    # y[m] <- y[m] + 1
                  }
                  if (y[m] < patch.size / 2 - lat * patch.size) {
                    y[m] <- y[m] + (patch.size * NS.patches)
                    # distance penalty for hitting an edge
                    # y[m] <- y[m] - 1
                  }
                }
              }
              # convert movement distances into numbers of grid cells (assume fish start in centre of cell):
              xy <- as.data.frame(cbind(round(x), round(y)))
              freq <- xy %>%
                group_by(V1, V2) %>%
                summarize(count = n())
              freq2D <- as.data.frame(array(0, c(length(unique(xy$y)), length(unique(xy$x)))))
              names(freq2D) <- sort(unique(xy$x))
              row.names(freq2D) <- sort(unique(xy$y))
              for (row in 1:length(freq$x)) {
                freq2D[as.character(freq$y[row]), as.character(freq$x[row])] <- freq$count[row]
              }
              # populate the move.array with movers (and stayers)
              for (xx in 1:length(unique(xy$x))) {
                for (yy in 1:length(unique(xy$y))) {
                  move.array[lat + as.numeric(row.names(freq2D)[yy]), lon + as.numeric(names(freq2D)[xx]), i, j, k] <- move.array[lat + as.numeric(row.names(freq2D)[yy]), lon + as.numeric(names(freq2D)[xx]), i, j, k] + freq2D[yy, xx]
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
}

move_ifelse <- function(pop) {
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
            if (pop[lat, lon, i, j, k] > 0) {
              # movers are subtracted from the present grid cell
              move.array[lat, lon, i, j, k] <- move.array[lat, lon, i, j, k] - pop[lat, lon, i, j, k]
              # determine the distribution of movement distances in nautical miles:
              dist <- rnbinom(pop[lat, lon, i, j, k], mu = mean.dist, size = Herit)
              dist <- ifelse(dist > EW.patches, EW.patches, dist)
              # determine the direction of each move
              theta <- runif(pop[lat, lon, i, j, k], 0, 2 * pi)
              # bias this movement in the north-south direction (along coasts) if this is a great white shark simulation (otherwise, comment out the next three lines):
              f.adj <- function(x, u) x - cos(x) * sin(x) - u
              my.uniroot <- function(x) uniroot(f.adj, c(0, 2 * pi), tol = 0.0001, u = x)$root
              theta <- vapply(theta, my.uniroot, numeric(1))
              # convert direction and distance into a distance in the x-direction (longitude)
              x <- cos(theta) * dist
              # bounce off edges
              x_bool <- ifelse((x <= patch.size * (EW.patches - lon) + patch.size / 2) & (x >= patch.size / 2 - lon * patch.size), TRUE, FALSE)
              x <- ifelse((x_bool == FALSE) & (x > patch.size * (EW.patches - lon) + patch.size / 2), (-(x - 2 * (patch.size * (EW.patches - lon) + patch.size / 2))), x)
              x <- ifelse((x_bool == FALSE) & (x < patch.size / 2 - lon * patch.size), (-(x - 2 * (patch.size / 2 - lon * patch.size))), x)
              # convert direction and distance into a distance in the y-direction (latitude)
              y <- sin(theta) * dist
              # bounce off edges
              y_bool <- ifelse((y <= patch.size * (NS.patches - lat) + patch.size / 2) & (y >= patch.size / 2 - lat * patch.size), TRUE, FALSE)
              y <- ifelse((y_bool == FALSE) & (y > patch.size * (NS.patches - lat) + patch.size / 2), (y - (patch.size * NS.patches)), y)
              y <- ifelse((y_bool == FALSE) & (y < patch.size / 2 - lat * patch.size), (y + (patch.size * NS.patches)), y)
              # convert movement distances into numbers of grid cells (assume fish start in centre of cell):
              x <- round(x)
              y <- round(y)
              xy <- as.data.frame(cbind(x, y))
              freq <<- xy %>%
                group_by(x, y) %>%
                summarize(count = n())
              freq2D <- as.data.frame(array(0, c(length(unique(xy$y)), length(unique(xy$x)))))
              names(freq2D) <- sort(unique(xy$x))
              row.names(freq2D) <- sort(unique(xy$y))
              for (row in 1:length(freq$x)) {
                freq2D[as.character(freq$y[row]), as.character(freq$x[row])] <- freq$count[row]
              }
              # populate the move.array with movers (and stayers)
              for (xx in 1:length(unique(xy$x))) {
                for (yy in 1:length(unique(xy$y))) {
                  move.array[lat + as.numeric(row.names(freq2D)[yy]), lon + as.numeric(names(freq2D)[xx]), i, j, k] <- move.array[lat + as.numeric(row.names(freq2D)[yy]), lon + as.numeric(names(freq2D)[xx]), i, j, k] + freq2D[yy, xx]
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


# plasticity move ---------------------------------------------------------


move <- function(pop, gen) {
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
            if (reserve.patches[, , gen] == 1) {
              mean.dist <- mean.dist - 1
            }
            if (pop[lat, lon, i, j, k] > 0) {
              # movers are subtracted from the present grid cell
              move.array[lat, lon, i, j, k] <- move.array[lat, lon, i, j, k] - pop[lat, lon, i, j, k]
              # determine the distribution of movement distances in nautical miles:
              dist <- rnbinom(pop[lat, lon, i, j, k], mu = mean.dist, size = Herit)
              # determine the direction of each move
              theta <- runif(pop[lat, lon, i, j, k], 0, 2 * pi)
              # bias this movement in the north-south direction (along coasts) if this is a great white shark simulation (otherwise, comment out the next three lines):
              f.adj <- function(x, u) x - cos(x) * sin(x) - u
              my.uniroot <- function(x) uniroot(f.adj, c(0, 2 * pi), tol = 0.0001, u = x)$root
              theta <- vapply(theta, my.uniroot, numeric(1))
              # convert direction and distance into a distance in the x-direction (longitude)
              x <- cos(theta) * dist
              # bounce off edges (assume fish start in centre of cell)
              for (m in 1:length(x)) {
                miss_x.edges <- FALSE
                while (miss_x.edges == FALSE) {
                  if (x[m] <= patch.size * (EW.patches - lon) + patch.size / 2) {
                    if (x[m] >= patch.size / 2 - lon * patch.size) {
                      miss_x.edges <- TRUE
                    }
                  }
                  if (x[m] > patch.size * (EW.patches - lon) + patch.size / 2) {
                    x[m] <- -(x[m] - 2 * (patch.size * (EW.patches - lon) + patch.size / 2))
                    # distance penalty for hitting an edge
                    # x[m] <- x[m] + 1
                  }
                  if (x[m] < patch.size / 2 - lon * patch.size) {
                    x[m] <- -(x[m] - 2 * (patch.size / 2 - lon * patch.size))
                    # distance penalty for hitting an edge
                    # x[m] <- x[m] - 1
                  }
                }
              }
              # convert direction and distance into a distance in the y-direction (latitude)
              y <- sin(theta) * dist
              # bounce off edges (assume fish start in centre of cell)
              for (m in 1:length(y)) {
                miss_y.edges <- FALSE
                while (miss_y.edges == FALSE) {
                  if (y[m] <= patch.size * (NS.patches - lat) + patch.size / 2) {
                    if (y[m] >= patch.size / 2 - lat * patch.size) {
                      miss_y.edges <- TRUE
                    }
                  }
                  if (y[m] > patch.size * (NS.patches - lat) + patch.size / 2) {
                    y[m] <- y[m] - (patch.size * NS.patches)
                    # distance penalty for hitting an edge
                    # y[m] <- y[m] + 1
                  }
                  if (y[m] < patch.size / 2 - lat * patch.size) {
                    y[m] <- y[m] + (patch.size * NS.patches)
                    # distance penalty for hitting an edge
                    # y[m] <- y[m] - 1
                  }
                }
              }
              # convert movement distances into numbers of grid cells (assume fish start in centre of cell):
              xy <- as.data.frame(cbind(x, y))
              hx <- hist(xy$x, breaks = seq(round(min(x), -2) - patch.size / 2, round(max(x), -2) + patch.size / 2, patch.size), plot = FALSE)
              xbreaks <- hx$breaks
              xmids <- hx$mids
              hy <- hist(xy$y, breaks = seq(round(min(y), -2) - patch.size / 2, round(max(y), -2) + patch.size / 2, patch.size), plot = FALSE)
              ybreaks <- hy$breaks
              ymids <- hy$mids
              freq <- as.data.frame(table(cut(xy[, 2], ybreaks, labels = ymids), cut(xy[, 1], xbreaks, labels = xmids)))
              freq2D <- as.data.frame(array(0, c(length(ymids), length(xmids))))
              names(freq2D) <- xmids / patch.size
              row.names(freq2D) <- ymids / patch.size
              freq2D[cbind(freq[, 1], freq[, 2])] <- freq[, 3]
              # populate the move.array with movers (and stayers)
              for (xx in 1:length(xmids)) {
                for (yy in 1:length(ymids)) {
                  move.array[lat + as.numeric(row.names(freq2D)[yy]), lon + as.numeric(names(freq2D)[xx]), i, j, k] <- move.array[lat + as.numeric(row.names(freq2D)[yy]), lon + as.numeric(names(freq2D)[xx]), i, j, k] + freq2D[yy, xx]
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
}



# Old Move 7/23/2022 ------------------------------------------------------

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
            if (pop[lat, lon, i, j, k] > 0) {
              # movers are subtracted from the present grid cell
              move.array[lat, lon, i, j, k] <- move.array[lat, lon, i, j, k] - pop[lat, lon, i, j, k]
              # determine the distribution of movement distances in nautical miles:
              dist <- rnbinom(pop[lat, lon, i, j, k], mu = mean.dist, size = Herit)
              # determine the direction of each move
              theta <- runif(pop[lat, lon, i, j, k], 0, 2 * pi)
              # bias this movement in the north-south direction (along coasts) if this is a great white shark simulation (otherwise, comment out the next three lines):
              f.adj <- function(x, u) x - cos(x) * sin(x) - u
              my.uniroot <- function(x) uniroot(f.adj, c(0, 2 * pi), tol = 0.0001, u = x)$root
              theta <- vapply(theta, my.uniroot, numeric(1))
              # convert direction and distance into a distance in the x-direction (longitude)
              x <- cos(theta) * dist
              # bounce off edges (assume fish start in centre of cell)
              for (m in 1:length(x)) {
                miss_x.edges <- FALSE
                while (miss_x.edges == FALSE) {
                  if (x[m] <= patch.size * (EW.patches - lon) + patch.size / 2) {
                    if (x[m] >= patch.size / 2 - lon * patch.size) {
                      miss_x.edges <- TRUE
                    }
                  }
                  if (x[m] > patch.size * (EW.patches - lon) + patch.size / 2) {
                    x[m] <- -(x[m] - 2 * (patch.size * (EW.patches - lon) + patch.size / 2))
                    # distance penalty for hitting an edge
                    # x[m] <- x[m] + 1
                  }
                  if (x[m] < patch.size / 2 - lon * patch.size) {
                    x[m] <- -(x[m] - 2 * (patch.size / 2 - lon * patch.size))
                    # distance penalty for hitting an edge
                    # x[m] <- x[m] - 1
                  }
                }
              }
              # convert direction and distance into a distance in the y-direction (latitude)
              y <- sin(theta) * dist
              # bounce off edges (assume fish start in centre of cell)
              for (m in 1:length(y)) {
                miss_y.edges <- FALSE
                while (miss_y.edges == FALSE) {
                  if (y[m] <= patch.size * (NS.patches - lat) + patch.size / 2) {
                    if (y[m] >= patch.size / 2 - lat * patch.size) {
                      miss_y.edges <- TRUE
                    }
                  }
                  if (y[m] > patch.size * (NS.patches - lat) + patch.size / 2) {
                    y[m] <- y[m] - (patch.size * NS.patches)
                    # distance penalty for hitting an edge
                    # y[m] <- y[m] + 1
                  }
                  if (y[m] < patch.size / 2 - lat * patch.size) {
                    y[m] <- y[m] + (patch.size * NS.patches)
                    # distance penalty for hitting an edge
                    # y[m] <- y[m] - 1
                  }
                }
              }
              # convert movement distances into numbers of grid cells (assume fish start in centre of cell):
              xy <- as.data.frame(cbind(x, y))
              hx <- hist(xy$x, breaks = seq(round(min(x), -2) - patch.size / 2, round(max(x), -2) + patch.size / 2, patch.size), plot = FALSE)
              xbreaks <- hx$breaks
              xmids <- hx$mids
              hy <- hist(xy$y, breaks = seq(round(min(y), -2) - patch.size / 2, round(max(y), -2) + patch.size / 2, patch.size), plot = FALSE)
              ybreaks <- hy$breaks
              ymids <- hy$mids
              freq <- as.data.frame(table(cut(xy[, 2], ybreaks, labels = ymids), cut(xy[, 1], xbreaks, labels = xmids)))
              freq2D <- as.data.frame(array(0, c(length(ymids), length(xmids))))
              names(freq2D) <- xmids / patch.size
              row.names(freq2D) <- ymids / patch.size
              freq2D[cbind(freq[, 1], freq[, 2])] <- freq[, 3]
              # populate the move.array with movers (and stayers)
              for (xx in 1:length(xmids)) {
                for (yy in 1:length(ymids)) {
                  move.array[lat + as.numeric(row.names(freq2D)[yy]), lon + as.numeric(names(freq2D)[xx]), i, j, k] <- move.array[lat + as.numeric(row.names(freq2D)[yy]), lon + as.numeric(names(freq2D)[xx]), i, j, k] + freq2D[yy, xx]
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
}


# older functions ---------------------------------------------------------


calc_temp_mortality2 <- function(SST, opt.temp, temp.range, s) {
  m <- 1 - exp((-(SST - opt.temp)^2) / (temp.range^2)) # temperature based mortality function from Walsworth et al.
  m <- 1 - m
  if (m > s) {
    nat.m <- s
  } else if (m < s) {
    nat.m <- m
  }
  return(nat.m)
}

############################################################################
## This function determines density dependent survival proportion for babies

survival_b2 <- function(num, SST) {
  s <- calc_temp_mortality2(SST, opt.temp, temp.range, sb)
  dd <- dd # density dependence of survival
  result <- s / (1 + dd * num)
  return(result)
}

############################################################################
## This function determines density dependent survival proportion for juveniles and adults

survival2 <- function(SST) {
  result <- calc_temp_mortality2(SST, opt.temp, temp.range, s)
  return(result)
}

recruit2 <- function(pop) {
  recruit.array <- world
  for (lat in 1:NS.patches) {
    for (lon in 1:EW.patches) {
      SST <- SST.patches[lat, lon, t]
      for (i in 1:NUM.age.classes) {
        if (i == 1) {
          # Some babies survive and recruit to juvenile age class
          s1 <- survival_b2(sum(pop[lat, lon, i, , ]), SST)
          s <- survival2(SST)
          for (j in 1:NUM.sexes) {
            for (k in 1:NUM.genotypes) {
              if (pop[lat, lon, i, j, k] > 0) {
                recruit.array[lat, lon, i + 1, j, k] <- recruit.array[lat, lon, i + 1, j, k] + rbinom(1, pop[lat, lon, i, j, k], s1)
              }
            }
          }
        }
        if (i == 2) {
          # Some juveniles survive
          for (j in 1:NUM.sexes) {
            for (k in 1:NUM.genotypes) {
              if (pop[lat, lon, i, j, k] > 0) {
                juvies.surviving <- rbinom(1, pop[lat, lon, i, j, k], s)
                # Some juveniles recruit to adult age class
                juvies.recruiting <- rbinom(1, juvies.surviving, p)
                recruit.array[lat, lon, i + 1, j, k] <- recruit.array[lat, lon, i + 1, j, k] + juvies.recruiting
                # The rest of the juveniles remain in the juvenile age class
                recruit.array[lat, lon, i, j, k] <- recruit.array[lat, lon, i, j, k] + juvies.surviving - juvies.recruiting
              }
            }
          }
        }
        if (i == 3) {
          # Some adults survive
          for (j in 1:NUM.sexes) {
            for (k in 1:NUM.genotypes) {
              if (pop[lat, lon, i, j, k] > 0) {
                recruit.array[lat, lon, i, j, k] <- recruit.array[lat, lon, i, j, k] + rbinom(1, pop[lat, lon, i, j, k], s)
              }
            }
          }
        }
      }
    }
  }
  return(recruit.array)
}

fishing2 <- function(pop, gen) {
  if (gen <= pre.reserve.gens + pre.fishing.gens) {
    each.patch.pop <- array(0, c(NS.patches, EW.patches))
    for (i in 2:NUM.age.classes) {
      for (j in 1:NUM.sexes) {
        for (k in 1:NUM.genotypes) {
          each.patch.pop[, ] <- each.patch.pop[, ] + pop[, , i, j, k]
        }
      }
    }
    mean.per.patch.pop <- mean(each.patch.pop)
    ff <- mean.per.patch.pop * (1 / fished - 1)
    if (mean.per.patch.pop > 0) {
      for (lat in 1:NS.patches) {
        for (lon in 1:EW.patches) {
          patch.pop <- sum(pop[lat, lon, c(2, 3), , ])
          if (patch.pop > 0) {
            f <- patch.pop / (ff + patch.pop)
            for (i in 2:NUM.age.classes) {
              for (j in 1:NUM.sexes) {
                for (k in 1:NUM.genotypes) {
                  pop[lat, lon, i, j, k] <- rbinom(1, pop[lat, lon, i, j, k], (1 - f))
                }
              }
            }
          }
        }
      }
    }
  }
  if (gen > pre.reserve.gens + pre.fishing.gens) {
    reserve.area <- sum(reserve.patches) / (NS.patches * EW.patches)
    fished.adj <- fished * 1 / (1 - reserve.area)
    each.patch.pop <- array(0, c(NS.patches, EW.patches))
    for (i in 2:NUM.age.classes) {
      for (j in 1:NUM.sexes) {
        for (k in 1:NUM.genotypes) {
          each.patch.pop[, ] <- each.patch.pop[, ] + pop[, , i, j, k]
        }
      }
    }
    for (lat in 1:NS.patches) {
      for (lon in 1:EW.patches) {
        if (reserve.patches[lat, lon] == 1) {
          each.patch.pop[lat, lon] <- NaN
        }
      }
    }
    mean.per.patch.pop <- mean(each.patch.pop, na.rm = TRUE)
    ff <- mean.per.patch.pop * (1 / fished.adj - 1)
    if (mean.per.patch.pop > 0) {
      for (lat in 1:NS.patches) {
        for (lon in 1:EW.patches) {
          if (reserve.patches[lat, lon] == 0) {
            patch.pop <- sum(pop[lat, lon, c(2, 3), , ])
            if (patch.pop > 0) {
              f <- patch.pop / (ff + patch.pop)
              for (i in 2:NUM.age.classes) {
                for (j in 1:NUM.sexes) {
                  for (k in 1:NUM.genotypes) {
                    pop[lat, lon, i, j, k] <- rbinom(1, pop[lat, lon, i, j, k], (1 - f))
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

spawn2 <- function(pop) {
  fec <- fecundity

  for (lat in 1:NS.patches) {
    for (lon in 1:EW.patches) {
      num.females <- sum(pop[lat, lon, 3, 1, ])
      num.males <- sum(pop[lat, lon, 3, 2, ])
      # Spawning only occurs if there is at least one males and one females in the patch
      if (num.females > 0 && num.males > 0) {
        # All females produce the same mean number of eggs
        NUM.A.eggs <- rpois(1, fec * pop[lat, lon, 3, 1, 1] + fec * pop[lat, lon, 3, 1, 2] / 2)
        NUM.a.eggs <- rpois(1, fec * pop[lat, lon, 3, 1, 3] + fec * pop[lat, lon, 3, 1, 2] / 2)
        # Males produce sperm in proportion to their genotypes
        freq.A.sperm <- pop[lat, lon, 3, 2, 1] / num.males + (pop[lat, lon, 3, 2, 2] / num.males) / 2
        freq.a.sperm <- pop[lat, lon, 3, 2, 3] / num.males + (pop[lat, lon, 3, 2, 2] / num.males) / 2
        # Sperm fertilize eggs in proportion to sperm genotype frequencies
        AA <- rbinom(1, NUM.A.eggs, freq.A.sperm)
        aa <- rbinom(1, NUM.a.eggs, freq.a.sperm)
        Aa <- NUM.A.eggs + NUM.a.eggs - AA - aa
        # Divide zygotes 50:50 among the sexes
        AA.f <- rbinom(1, AA, 0.5)
        AA.m <- AA - AA.f
        Aa.f <- rbinom(1, Aa, 0.5)
        Aa.m <- Aa - Aa.f
        aa.f <- rbinom(1, aa, 0.5)
        aa.m <- aa - aa.f
        # Female babies
        pop[lat, lon, 1, 1, 1] <- pop[lat, lon, 1, 1, 1] + AA.f
        pop[lat, lon, 1, 1, 2] <- pop[lat, lon, 1, 1, 2] + Aa.f
        pop[lat, lon, 1, 1, 3] <- pop[lat, lon, 1, 1, 3] + aa.f
        # Male babies
        pop[lat, lon, 1, 2, 1] <- pop[lat, lon, 1, 2, 1] + AA.m
        pop[lat, lon, 1, 2, 2] <- pop[lat, lon, 1, 2, 2] + Aa.m
        pop[lat, lon, 1, 2, 3] <- pop[lat, lon, 1, 2, 3] + aa.m
      }
    }
  }
  return(pop)
}

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
            if (pop[lat, lon, i, j, k] > 0) {
              # movers are subtracted from the present grid cell
              move.array[lat, lon, i, j, k] <- move.array[lat, lon, i, j, k] - pop[lat, lon, i, j, k]
              # determine the distribution of movement distances in nautical miles:
              dist <- rnbinom(pop[lat, lon, i, j, k], mu = mean.dist, size = Herit)
              # determine the direction of each move
              theta <- runif(pop[lat, lon, i, j, k], 0, 2 * pi)
              # bias this movement in the north-south direction (along coasts) if this is a great white shark simulation (otherwise, comment out the next three lines):
              # f.adj <- function(x, u) x-cos(x)*sin(x) - u
              # my.uniroot <- function(x) uniroot(f.adj, c(0, 2*pi), tol = 0.0001, u = x)$root
              # theta <- vapply(theta, my.uniroot, numeric(1))
              # convert direction and distance into a distance in the x-direction (longitude)
              x <- cos(theta) * dist
              # bounce off edges (assume fish start in centre of cell)
              for (m in 1:length(x)) {
                miss_x.edges <- FALSE
                while (miss_x.edges == FALSE) {
                  if (x[m] <= patch.size * (EW.patches - lon) + patch.size / 2) {
                    if (x[m] >= patch.size / 2 - lon * patch.size) {
                      miss_x.edges <- TRUE
                    }
                  }
                  if (x[m] > patch.size * (EW.patches - lon) + patch.size / 2) {
                    x[m] <- -(x[m] - 2 * (patch.size * (EW.patches - lon) + patch.size / 2))
                    # distance penalty for hitting an edge
                    # x[m] <- x[m] + 1
                  }
                  if (x[m] < patch.size / 2 - lon * patch.size) {
                    x[m] <- -(x[m] - 2 * (patch.size / 2 - lon * patch.size))
                    # distance penalty for hitting an edge
                    # x[m] <- x[m] - 1
                  }
                }
              }
              # convert direction and distance into a distance in the y-direction (latitude)
              y <- sin(theta) * dist
              # bounce off edges (assume fish start in centre of cell)
              for (m in 1:length(y)) {
                miss_y.edges <- FALSE
                while (miss_y.edges == FALSE) {
                  if (y[m] <= patch.size * (NS.patches - lat) + patch.size / 2) {
                    if (y[m] >= patch.size / 2 - lat * patch.size) {
                      miss_y.edges <- TRUE
                    }
                  }
                  if (y[m] > patch.size * (NS.patches - lat) + patch.size / 2) {
                    y[m] <- -(y[m] - 2 * (patch.size * (NS.patches - lat) + patch.size / 2))
                    # distance penalty for hitting an edge
                    # y[m] <- y[m] + 1
                  }
                  if (y[m] < patch.size / 2 - lat * patch.size) {
                    y[m] <- -(y[m] - 2 * (patch.size / 2 - lat * patch.size))
                    # distance penalty for hitting an edge
                    # y[m] <- y[m] - 1
                  }
                }
              }
              # convert movement distances into numbers of grid cells (assume fish start in centre of cell):
              xy <- as.data.frame(cbind(x, y))
              hx <- hist(xy$x, breaks = seq(round(min(x), -2) - patch.size / 2, round(max(x), -2) + patch.size / 2, patch.size), plot = FALSE)
              xbreaks <- hx$breaks
              xmids <- hx$mids
              hy <- hist(xy$y, breaks = seq(round(min(y), -2) - patch.size / 2, round(max(y), -2) + patch.size / 2, patch.size), plot = FALSE)
              ybreaks <- hy$breaks
              ymids <- hy$mids
              freq <- as.data.frame(table(cut(xy[, 2], ybreaks, labels = ymids), cut(xy[, 1], xbreaks, labels = xmids)))
              freq2D <- as.data.frame(array(0, c(length(ymids), length(xmids))))
              names(freq2D) <- xmids / patch.size
              row.names(freq2D) <- ymids / patch.size
              freq2D[cbind(freq[, 1], freq[, 2])] <- freq[, 3]
              # populate the move.array with movers (and stayers)
              for (xx in 1:length(xmids)) {
                for (yy in 1:length(ymids)) {
                  move.array[lat + as.numeric(row.names(freq2D)[yy]), lon + as.numeric(names(freq2D)[xx]), i, j, k] <- move.array[lat + as.numeric(row.names(freq2D)[yy]), lon + as.numeric(names(freq2D)[xx]), i, j, k] + freq2D[yy, xx]
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
}
