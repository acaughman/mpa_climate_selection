# Setup -------------------------------------------------------------------

library(RangeShiftR)
library(raster)
library(RColorBrewer)
library(rasterVis)
library(latticeExtra)
library(viridis)
library(grid)
library(gridExtra)

# relative path from working directory:
dirpath = "Tutorial_04/"

dir.create(paste0(dirpath,"Inputs"), showWarnings = TRUE)
dir.create(paste0(dirpath,"Outputs"), showWarnings = TRUE)
dir.create(paste0(dirpath,"Output_Maps"), showWarnings = TRUE)


# Simulation Parameters ---------------------------------------------------

#create landscape
land <- ArtificialLandscape(propSuit = 0.3,
                            K_or_DensDep = 100,
                            dimX = 50,
                            dimY = 800,
                            continuous = FALSE)

#simulation parameters
sim_0 <- Simulation(Simulation = 1,
                    Replicates = 5,
                    Years = 1300,
                    Gradient = 1,
                    GradSteep = 0.02,
                    Optimum = 100,
                    f = 0.0,
                    Shifting = TRUE,
                    ShiftRate = 1,
                    ShiftStart = 500,
                    ShiftEnd = 800,
                    OutIntPop = 0,
                    OutIntInd = 50,
                    OutIntRange = 5,
                    OutIntTraitRow = 50)

#demography parameters
demo <- Demography(Rmax = 4.0)

#genetic parameters
gene <- Genetics(Architecture = 0, 
                 NLoci = 3,
                 ProbMutn = 0.001,
                 MutationSD = 1.0,
                 ProbCross = 0.3,
                 AlleleSD = 0.1)

#initialization parameter
init <- Initialise(InitDens = 0,
                   maxY = 200)

#dispersal transfer parameter
trans_b <- DispersalKernel(IndVar = TRUE,
                           Distances = matrix(c(250, 50), ncol = 2),
                           TraitScaleFactor = 50)

#dispersal parameter
disp_b <- Dispersal(Emigration = Emigration(EmigProb = 0.1), 
                    Transfer   = trans_b, 
                    Settlement = Settlement() )

#simulation master
s <- RSsim(land = land, demog = demo, dispersal = disp_b, simul = sim_0, gene = gene, init = init)

# Evolution of dispersal distance -----------------------------------------

#run the simulation
RunRS(s, dirpath)


# Simulation Plots --------------------------------------------------------

# load 'traits by rows' output file
trait_ts <- read.table(paste0(dirpath,"Outputs/Batch1_Sim1_Land1_TraitsXrow.txt"), header = T)
# plot the time series for the replicate chosen above:
trait_ts_07 <- subset.data.frame(trait_ts, Rep==7)

plot(NULL, type = "n", ylab = "Mean dispersal distance [m]", xlab = "y coordinate", xlim=c(0, 500), ylim=c(150,450), main = "Trait time-series (Replicate 7)")
years <-  c(500, 650, 800, 1200)
cols <- hcl.colors(length(years), palette = "Dark 3", alpha = 1, rev = FALSE)
leg.txt <- c()
for(i in 1:length(years)) {
  trait_ts_07_yr <- subset.data.frame(trait_ts_07, Year==years[i])
  polygon(c(trait_ts_07_yr$y,rev(trait_ts_07_yr$y)),
          c((trait_ts_07_yr$mean_distI+trait_ts_07_yr$std_distI), rev(pmax(0,trait_ts_07_yr$mean_distI-trait_ts_07_yr$std_distI))),
          border=NA, col='grey80')
  lines(trait_ts_07_yr$y, trait_ts_07_yr$mean_distI, type = "l", lwd = 1, col = cols[i])
  leg.txt <- c(leg.txt, paste("Year", years[i]))
}
legend("bottomright", leg.txt, col = cols, lwd = 1.5)