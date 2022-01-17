# Setup -------------------------------------------------------------------

library(RangeShiftR)
library(raster)
library(RColorBrewer)
library(rasterVis)
library(latticeExtra)
library(viridis)
library(grid)
library(tidyverse)
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
                    Replicates = 6,
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
trait_ts <- read.table(paste0(dirpath,"Outputs/Batch1_Sim1_Land1_TraitsXrow.txt"), header = T) %>% 
  as.data.frame() %>% 
  mutate(Rep = as.factor(Rep)) %>% 
  filter(Year > 400 & Year < 900) %>% 
  mutate(Year = as.factor(Year))

ggplot(trait_ts, aes(y, mean_distI, color = Year)) + theme_bw() +
  geom_point() +
  facet_grid(~Rep)

