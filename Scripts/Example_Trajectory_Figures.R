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


#----------------------------------------------------------------------
# Figure showing example predator and prey trajectories
#----------------------------------------------------------------------

#Predator mass (g)
mass_pred <- 40000

#Prey mass (g)
mass_prey <- prey.mass(mass_pred)

#"Lifespan" and sampling interval for the simulations
t <- sampling(mass_prey, crossings = 30)


#Generate the prey movement models
#Generate the HR centres of the prey
CENTRES <- rbvpois(n = 10,
                   a = pred.SIG(mass_pred)*.75,
                   b = pred.SIG(mass_pred)*.75,
                   c = 0)
CENTRES <- scale(CENTRES, scale = FALSE)

PREY_mods <- list()
for(i in 1:10){
    # Prey movement parameters
    prey_tau_p <- mean(prey_details[which(prey_details$generation == max(prey_details$generation)),"tau_p"])
    prey_tau_v <- mean(prey_details[which(prey_details$generation == max(prey_details$generation)),"tau_v"])
    prey_sig <- mean(prey_details[which(prey_details$generation == max(prey_details$generation)),"sig"])
    
    PREY_mods[[i]] <- ctmm(tau = c(prey_tau_p,prey_tau_v),
                           mu = c(CENTRES[i,1],CENTRES[i,2]),
                           sigma = prey_sig)
} #Closes loop over n_prey


# Predator movement parameters
pred_tau_p <- mean(pred_details[which(pred_details$generation == max(pred_details$generation)),"tau_p"])
pred_tau_v <- mean(pred_details[which(pred_details$generation == max(pred_details$generation)),"tau_v"])
pred_sig <- mean(pred_details[which(pred_details$generation == max(pred_details$generation)),"sig"])

PRED_mod <- ctmm(tau = c(pred_tau_p,
                         pred_tau_v),
                 mu = c(0,0),
                 sigma = pred_sig)

#Simulate the prey movement
#Parallelised to speed up run times
PREY_tracks <- lapply(PREY_mods,
                      FUN = simulate,
                      t = t)

#Simulate the predator movement
PRED_tracks <- simulate(PRED_mod,t = t)

COLS <- viridis::viridis(10)

tracks <- 
    ggplot() +
    geom_path(aes(y=PREY_tracks[[1]]$y, x=PREY_tracks[[1]]$x), color =COLS[1], size = 0.1, alpha = 0.3) +
    geom_path(aes(y=PREY_tracks[[2]]$y, x=PREY_tracks[[2]]$x), color =COLS[2], size = 0.1, alpha = 0.3) +
    geom_path(aes(y=PREY_tracks[[3]]$y, x=PREY_tracks[[3]]$x), color =COLS[3], size = 0.1, alpha = 0.3) +
    geom_path(aes(y=PREY_tracks[[4]]$y, x=PREY_tracks[[4]]$x), color =COLS[4], size = 0.1, alpha = 0.3) +
    geom_path(aes(y=PREY_tracks[[5]]$y, x=PREY_tracks[[5]]$x), color =COLS[5], size = 0.1, alpha = 0.3) +
    geom_path(aes(y=PREY_tracks[[6]]$y, x=PREY_tracks[[6]]$x), color =COLS[6], size = 0.1, alpha = 0.3) +
    geom_path(aes(y=PREY_tracks[[7]]$y, x=PREY_tracks[[7]]$x), color =COLS[7], size = 0.1, alpha = 0.3) +
    geom_path(aes(y=PREY_tracks[[8]]$y, x=PREY_tracks[[8]]$x), color =COLS[8], size = 0.1, alpha = 0.3) +
    geom_path(aes(y=PREY_tracks[[9]]$y, x=PREY_tracks[[9]]$x), color =COLS[9], size = 0.1, alpha = 0.3) +
    geom_path(aes(y=PREY_tracks[[10]]$y, x=PREY_tracks[[10]]$x), color =COLS[10], size = 0.1, alpha = 0.3) +
    geom_path(aes(y=PRED_tracks$y, x=PRED_tracks$x), color = "black", size = 0.1, alpha = 0.8) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_text(size=8, family = "serif"),
          axis.title.x = element_text(size=8, family = "serif"),
          axis.text.y = element_text(size=6, family = "serif"),
          axis.text.x  = element_text(size=6, family = "serif"),
          plot.title = element_text(hjust = -0.05, size = 12, family = "serif"),
          legend.position = "none",
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent", color = NA)) +
    ylab("Y (m)") +
    xlab("X (m)")


ggsave(tracks,
       width = 3.23, height = 3, units = "in",
       dpi = 600,
       bg = "transparent",
       file="Results/Sim_Example.png")