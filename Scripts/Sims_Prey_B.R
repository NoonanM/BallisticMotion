#Script for studying the evolution of ballistic length scales
# in pred/prey systems via simulation

#Written by Michael Noonan

#Last updated: Jan 22nd 2022

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
# Set up the global parameters for the simulation
#----------------------------------------------------------------------

#Predator mass (g)
mass_pred <- 15000

#Prey mass (g)
mass_prey <- prey.mass(mass_pred)

#"Lifespan" and sampling interval for the simulations
t <- sampling(mass_prey)

#Number of preds & prey in each "arena"
n_prey <- 10
n_pred <- 1

#Number of "arenas" (becomes important when there are predators)
REPS <- 50

#Number of generations
GENS <- 250

# Build the raster of food patches for prey to feed on
FOOD <- patches(mass_pred, width = 25, pred = T)


#Lists for storing results and drawing params
prey_res <- list()
prey_details <- list()

#----------------------------------------------------------------------
# Run the simulation
#----------------------------------------------------------------------

#Loop over the number of generations
for(G in 1:GENS){
  
  #Empty lists for storing the results of the current generation
  prey <- list()
  
  for(R in 1:REPS){
    
    #Generate the prey movement models
    #If the first gen, generate movement parameters from the mass functions
    if(G == 1){
      
      PREY_mods <- list()
      for(i in 1:n_prey){
        # Prey movement parameters
        prey_tau_p <- prey.tau_p(mass_prey, variance = TRUE)
        prey_tau_v <- prey.tau_v(mass_prey, variance = TRUE)
        prey_sig <- prey.SIG(mass_prey)
        prey_lv <- sqrt((prey_tau_v/prey_tau_p)*prey_sig)
        
        PREY_mods[[i]] <- ctmm(tau = c(prey_tau_p,prey_tau_v),
                               mu = c(0,0),
                               sigma = prey_sig)
      } #Closes loop over n_prey
    } #Closes if gen 1 ask 
    
    
    #If any other generation, draw movement from pool of offspring params
    if(G != 1){
      
      PREY_mods <- list()
      for(i in 1:n_prey){
        # Prey movement parameters
        prey_tau_p <- sample(PREY_tau_p, 1) + rnorm(1, 0, 10) #Add some 'mutation' based variance
        prey_tau_p <- ctmm:::clamp(prey_tau_p, min = 0.1, max = Inf) #Clamp the minimum to 0
        prey_tau_v <- sample(PREY_tau_v, 1) + rnorm(1, 0, 2)  #Add some 'mutation' based variance
        prey_tau_v <- ctmm:::clamp(prey_tau_v, min = 0.1, max = Inf) #Clamp the minimum to 0
        prey_sig <- sample(PREY_sig, 1)
        prey_lv <- sqrt((prey_tau_v/prey_tau_p)*prey_sig)
        
        PREY_mods[[i]] <- ctmm(tau = c(prey_tau_p,prey_tau_v),
                               mu = c(0,0),
                               sigma = prey_sig)
      } #Closes loop over n_prey
    } #Closes if not gen 1 ask 
    
    
    #Simulate the prey movement
    # PREY_tracks <- list()
    # for(i in 1:n_prey){
    #   PREY_tracks[[i]] <- simulate(PREY_mods[[i]], t = t)
    # }
    
    #Parallelised version to speed up run times
    PREY_tracks <- mclapply(PREY_mods,
                            FUN = simulate,
                            t = t,
                            mc.cores = 6)
    
    # Calculate prey benefits
    benefits_prey <- vector()
    for(i in 1:n_prey){
      benefits_prey[i] <- grazing(PREY_tracks[[i]], FOOD)
      
    }
    
    #Calculate prey fitness
    offspring_prey <- prey.fitness(benefits_prey,
                                   mass_prey,
                                   models = PREY_mods,
                                   calories = 12)
    
    
    #Get the values of the prey movement model parameters
    prey_lvs <- vector()
    prey_TAU_V <- vector() 
    prey_TAU_P <- vector()
    prey_SIGMA <- vector()
    prey_SPEED <- vector()
    for(i in 1:n_prey){
      prey_TAU_V[i] <- PREY_mods[[i]]$tau["velocity"]
      prey_TAU_P[i] <- PREY_mods[[i]]$tau["position"]
      prey_SIGMA[i] <- ctmm:::area.covm(PREY_mods[[i]]$sigma)
      prey_SPEED[i] <- if(nrow(summary(PREY_mods[[i]], units = FALSE)$CI)==4){summary(PREY_mods[[i]], units = FALSE)$CI[4,2]} else{Inf}
      prey_lvs[i] <- sqrt((prey_TAU_V[i]/prey_TAU_P[i])*prey_SIGMA[i])
    }
    
    #Summarise the results of the prey 
    prey[[R]] <- data.frame(generation = G,
                            tau_p = prey_TAU_P,
                            tau_v = prey_TAU_V,
                            sig = prey_SIGMA,
                            speed = prey_SPEED,
                            lv = prey_lvs,
                            patches = benefits_prey,
                            offspring = offspring_prey)
  }
  
  #Compile the results from the generation
  prey <- do.call(rbind, prey)
  
  # Save the results
  # prey
  prey_res[[G]] <- data.frame(generation = G,
                              lv = mean(prey$lv),
                              var = var(prey$lv))
  
  prey_details[[G]] <- prey
  
  
  #Set up the parameters for the next generation based on
  #Fitness of current generation
  PREY_tau_p <- vector()
  PREY_tau_v <- vector()
  PREY_sig <- vector()
  for(i in 1:nrow(prey)){
    if(prey[i,"offspring"] >0){
      PREY_tau_p <- c(PREY_tau_p,
                      rep(prey[i,"tau_p"], prey[i,"offspring"]))
      
      PREY_tau_v <- c(PREY_tau_v,
                      rep(prey[i,"tau_v"], prey[i,"offspring"]))
      
      PREY_sig <- c(PREY_sig,
                    rep(prey[i,"sig"], prey[i,"offspring"]))
      
    } #Closes the if statement
  } #Closes the loop over the number of prey
  
  print(G)
  save(prey_res, file = 'Results/lv_Evo_Prey_Ben.Rda')
  save(prey_details, file = 'Results/lv_Evo_Prey_Ben_details.Rda')
}

#Compile the results
load('Results/lv_Evo_Prey_Ben.Rda')
load('Results/lv_Evo_Prey_Ben_details.Rda')

prey_res <- do.call(rbind, prey_res)
prey_details <- do.call(rbind, prey_details)

#Compare with starting values
PREY_LV <- prey_res$lv[1] #(sqrt((tau_v/tau_p)*sig))

#par(mfrow = c(1,2))
plot(lv/PREY_LV ~ generation,
     data = prey_res,
     type = "l",
     col = "red",
     ylim = c(0,1.2),
     ylab = "Relative change in lv",
     xlab = "Generation")
abline(h = 1, col = 'grey30', lty = 'dashed')

polygon(c(prey_res$generation, rev(prey_res$generation)), c((prey_res$lv - prey_res$var)/PREY_LV, rev((prey_res$lv + prey_res$var)/PREY_LV)),
        col = adjustcolor("red",alpha.f=0.5), border=NA)
lines(lv/PREY_LV ~ generation,
      data = prey_res,
      col = "red")


plot(offspring ~ lv,
     data = prey_details,
     #type = "l",
     col = "red",
     #ylim = c(0,7),
     xlab = "Ballistic lengthscale (m)",
     ylab = "Offspring")
abline(lm(offspring ~ lv,
          data = prey_details), lty = 'dashed')


plot(offspring ~ speed,
     data = prey_details,
     #type = "l",
     col = "red",
     #ylim = c(0,7),
     xlab = "Movement speed (m/s)",
     ylab = "Offspring")
abline(lm(offspring ~ speed,
          data = prey_details), lty = 'dashed')


plot(speed ~ lv,
     data = prey_details,
     #type = "l",
     col = "red",
     #ylim = c(0,7),
     xlab = "Ballistic lengthscale (m)",
     ylab = "Movement speed (m/s)")


plot(patches ~ generation,
     data = prey_details,
     #type = "l",
     col = "red",
     #ylim = c(0,7),
     xlab = "Generation",
     ylab = "Patches")
abline(lm(patches ~ generation,
          data = prey_details), lty = 'dashed')


plot(lv ~ generation,
     data = prey_details,
     #type = "l",
     col = "red",
     #ylim = c(0,7),
     xlab = "Generation",
     ylab = "Ballistic lengthscale (m)")
abline(lm(lv ~ generation,
          data = prey_details), lty = 'dashed')

plot(tau_v ~ generation,
     data = prey_details,
     #type = "l",
     col = "red",
     #ylim = c(0,7),
     xlab = "Generation",
     ylab = "Directional persistance (sec)")
abline(lm(tau_v ~ generation,
          data = prey_details), lty = 'dashed')


plot(tau_p ~ generation,
     data = prey_details,
     #type = "l",
     col = "red",
     #ylim = c(0,7),
     xlab = "Generation",
     ylab = "Range crossing time (sec)")
abline(lm(tau_p ~ generation,
          data = prey_details), lty = 'dashed')

plot(speed ~ generation,
     data = prey_details,
     #type = "l",
     col = "red",
     #ylim = c(0,7),
     xlab = "Generation",
     ylab = "Movement speed (m/s)")
abline(lm(speed ~ generation,
          data = prey_details), lty = 'dashed')


