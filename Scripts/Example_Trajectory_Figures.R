#Script for studying the evolution of ballistic length scales
# in pred/prey systems via simulation. This script is used to generate example figures
# of the movement used in the simulation study

#Written by Michael Noonan

#Last updated: Mar 24th 2022

#----------------------------------------------------------------------
# Preamble
#----------------------------------------------------------------------

#Set the working directory
setwd("~/Dropbox (Personal)/UBC/Projects/BallisticMotion")

#Set the random seed
set.seed(1)

#Import necessary packages
library(extraDistr)
library(parallel)

#Source the functions
source("Scripts/Functions.R")


#----------------------------------------------------------------------
# Figure showing an example prey trajectory and foraging patches
#----------------------------------------------------------------------

#Predator mass (g)
mass_pred <- 40000

#Prey mass (g)
mass_prey <- prey.mass(mass_pred)

#"Lifespan" and sampling interval for the simulations
t <- sampling(mass_prey, crossings = 30)


# Build the raster of food patches for prey to feed on
FOOD <- patches(mass_pred,
                width = round(sqrt(pred.SIG(mass_prey))/10),
                pred = T)

#Energetic value
((10^(0.774 + 0.727*log10(mass_prey)))^1.22)/150
# Prey movement parameters
prey_tau_p <- prey.tau_p(mass_prey, variance = TRUE)
prey_tau_v <- prey.tau_v(mass_prey, variance = TRUE)
prey_sig <- prey.SIG(mass_prey)

#Define the movement model based on the parameter values
PREY_mod <- ctmm(tau = c(prey_tau_p,prey_tau_v),
                 mu = c(0,0),
                 sigma = prey_sig)

#Simulate the movement from the defined model
PREY_tracks <- simulate(PREY_mod,t = t)

#For defining the colours of the trajectories
COLS <- viridis::viridis(10)


png(filename="Results/Prey_Movement.png",
    width = 6.86, height = 6, units = "in",
    res = 600)
plot(PREY_tracks,
     type = "l",
     col = COLS[5])
plot(rasterToPolygons(FOOD),
     add=TRUE,
     border='black',
     lwd=0.5,
     xlim = c(-500,500),
     ylim = c(-500,500))
plot(PREY_tracks,
     type = "l",
     add=TRUE,
     col = COLS[5])
#title("Simulated movement and foraging patches for a kg prey species")
dev.off()