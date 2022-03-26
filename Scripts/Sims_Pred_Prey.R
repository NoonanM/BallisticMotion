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
t <- sampling(mass_prey, crossings = 30)

#Number of preds & prey in each "arena"
n_prey <- 10
n_pred <- 1

#Number of "arenas"
REPS <- 50

#Number of generations
GENS <- 2000

# Build the raster of food patches for prey to feed on
FOOD <- patches(mass_pred, width = 25, pred = T)

#Lists for storing results and drawing params
prey_res <- list()
prey_details <- list()

pred_res <- list()
pred_details <- list()

#----------------------------------------------------------------------
# Run the simulation
#----------------------------------------------------------------------

#Loop over the number of generations
for(G in 1702:GENS){
  
  #Empty lists for storing the results of the current generation
  prey <- list()
  pred <- list()
  
  for(R in 1:REPS){
    
    #Generate the prey movement models
    #If the first gen, generate movement parameters from the mass functions
    if(G == 1){
      #Generate the HR centres of the prey
      CENTRES <- rbvpois(n = n_prey,
                         a = pred.SIG(mass_pred)*.75,
                         b = pred.SIG(mass_pred)*.75,
                         c = 0)
      CENTRES <- scale(CENTRES, scale = FALSE)
      
      PREY_mods <- list()
      for(i in 1:n_prey){
        # Prey movement parameters
        prey_tau_p <- prey.tau_p(mass_prey, variance = TRUE)
        prey_tau_v <- prey.tau_v(mass_prey, variance = TRUE)
        prey_sig <- prey.SIG(mass_prey)
        prey_lv <- sqrt((prey_tau_v/prey_tau_p)*prey_sig)
        
        PREY_mods[[i]] <- ctmm(tau = c(prey_tau_p,prey_tau_v),
                               mu = c(CENTRES[i,1],CENTRES[i,2]),
                               sigma = prey_sig)
      } #Closes loop over n_prey
      
      PRED_mods <- list()
      for(i in 1:n_pred){
        # Predator movement parameters
        pred_tau_p <- pred.tau_p(mass_pred, variance = TRUE)
        pred_tau_v <- pred.tau_v(mass_pred, variance = TRUE)
        pred_sig <- pred.SIG(mass_pred)
        pred_lv <- sqrt((pred_tau_v/pred_tau_p)*pred_sig)
        
        PRED_mods[[i]] <- ctmm(tau = c(pred_tau_p,
                                       pred_tau_v),
                               mu = c(0,0),
                               sigma = pred_sig)
      } # Closes loop over n_pred
      
    } # Closes if gen 1 ask 
    
    
    # If any other generation, draw movement from pool of offspring params
    if(G != 1){
      # Generate the HR centres of the prey
      CENTRES <- rbvpois(n = n_prey,
                         a = pred.SIG(mass_pred)*.75,
                         b = pred.SIG(mass_pred)*.75,
                         c = 0)
      CENTRES <- scale(CENTRES, scale = FALSE)
      
      PREY_mods <- list()
      for(i in 1:n_prey){
        # Prey movement parameters
        prey_tau_p <- sample(PREY_tau_p, 1) + rnorm(1, 0, 10) #Add some 'mutation' based variance
        prey_tau_p <- ctmm:::clamp(prey_tau_p, min = 0.1, max = Inf) #Clamp the minimum to 0
        prey_tau_v <- sample(PREY_tau_v, 1) + rnorm(1, 0, 2)  #Add some 'mutation' based variance
        prey_tau_v <- ctmm:::clamp(prey_tau_v, min = 0.1, max = Inf) #Clamp the minimum to 0
        prey_sig <- sample(PREY_sig, 1)
        prey_lv <- sqrt((prey_tau_v/prey_tau_p)*prey_sig)
        
        # Define the movement models
        PREY_mods[[i]] <- ctmm(tau = c(prey_tau_p,prey_tau_v),
                               mu = c(CENTRES[i,1],CENTRES[i,2]),
                               sigma = prey_sig)
      } # Closes loop over n_prey
      
      PRED_mods <- list()
      for(i in 1:n_pred){
        # Predator movement parameters
        pred_tau_p <- sample(PRED_tau_p, 1) + rnorm(1, 0, 15) #Add some 'mutation' based variance
        pred_tau_p <- ctmm:::clamp(pred_tau_p, min = 0.1, max = Inf) #Clamp the minimum to 0
        pred_tau_v <- sample(PRED_tau_v, 1) + rnorm(1, 0, 5)  #Add some 'mutation' based variance
        pred_tau_v <- ctmm:::clamp(pred_tau_v, min = 0.1, max = Inf) #Clamp the minimum to 0
        pred_sig <- sample(PRED_sig, 1)
        pred_lv <- sqrt((pred_tau_v/pred_tau_p)*pred_sig)
        
        PRED_mods[[i]] <- ctmm(tau = c(pred_tau_p,
                                       pred_tau_v),
                               mu = c(0,0),
                               sigma = pred_sig)
      } #Closes loop over n_pred
    } #Closes if not gen 1 ask 
    
    
    #Simulate the prey movement
    #Parallelised to speed up run times
    PREY_tracks <- mclapply(PREY_mods,
                            FUN = simulate,
                            t = t,
                            mc.cores = 2)
    
    #Simulate the predator movement
    PRED_tracks <- list()
    for(i in 1:n_pred){
      PRED_tracks[[i]] <- simulate(PRED_mods[[i]],t = t)
    }
    
    # Calculate prey benefits
    benefits_prey <- vector()
    for(i in 1:n_prey){
      benefits_prey[i] <- grazing(PREY_tracks[[i]], FOOD)
    }
    
    #Count the encounters (only setup for single predator/arena)
    encounters <- encounter(prey.tracks = PREY_tracks,
                            pred.tracks = PRED_tracks,
                            range = 25)
    
    #Calculate prey fitness
    offspring_prey <- prey.fitness(benefits = benefits_prey,
                                   costs = encounters,
                                   mass_prey,
                                   crossings = 30,
                                   models = PREY_mods,
                                   calories = 15)
    
    #Calculate predator fitness
    offspring_pred <- pred.fitness(encounters = encounters,
                                   mass = mass_pred,
                                   models = PRED_mods,
                                   time = t)
    
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
                            encounter = encounters,
                            offspring = offspring_prey)
    
    
    #Get the values of the predator movement model parameters
    pred_lvs <- vector()
    pred_TAU_V <- vector() 
    pred_TAU_P <- vector()
    pred_SIGMA <- vector()
    pred_SPEED <- vector()
    for(i in 1:n_pred){
      pred_TAU_V[i] <- PRED_mods[[i]]$tau["velocity"]
      pred_TAU_P[i] <- PRED_mods[[i]]$tau["position"]
      pred_SIGMA[i] <- ctmm:::area.covm(PRED_mods[[i]]$sigma)
      pred_SPEED[i] <- if(nrow(summary(PRED_mods[[i]], units = FALSE)$CI)==4){summary(PRED_mods[[i]], units = FALSE)$CI[4,2]} else{Inf}
      pred_lvs[i] <- sqrt((pred_TAU_V[i]/pred_TAU_P[i])*pred_SIGMA[i])
    }
    
    #Summarise the results of the prey 
    pred[[R]] <- data.frame(generation = G,
                            tau_p = pred_TAU_P,
                            tau_v = pred_TAU_V,
                            sig = pred_SIGMA,
                            speed = pred_SPEED,
                            lv = pred_lvs,
                            encounter = sum(encounters),
                            offspring = offspring_pred)
    
  }
  
  #Compile the results from the generation
  prey <- do.call(rbind, prey)
  pred <- do.call(rbind, pred)
  
  # Save the results
  # prey
  prey_res[[G]] <- data.frame(generation = G,
                              lv = mean(prey$lv),
                              var = var(prey$lv),
                              pred_lv = mean(pred$lv))
  
  prey_details[[G]] <- prey
  
  # Predator
  pred_res[[G]] <- data.frame(generation = G,
                              lv = mean(pred$lv),
                              var = var(pred$lv))
  
  pred_details[[G]] <- pred
  
  
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
  
  
  
  #Set up the parameters for the next generation of predators to be based on
  #Fitness of current generation
  PRED_tau_p <- vector()
  PRED_tau_v <- vector()
  PRED_sig <- vector()
  for(i in 1:nrow(pred)){
    if(pred[i,"offspring"] >0){
      PRED_tau_p <- c(PRED_tau_p,
                      rep(pred[i,"tau_p"], pred[i,"offspring"]))
      
      PRED_tau_v <- c(PRED_tau_v,
                      rep(pred[i,"tau_v"], pred[i,"offspring"]))
      
      PRED_sig <- c(PRED_sig,
                    rep(pred[i,"sig"], pred[i,"offspring"]))
      
    } #Closes the if statement
  } #Closes the loop over the number of pred
  
  
  
  
  #Save prey results
  save(prey_res, file = 'Results/lv_Evo_Full_Prey.Rda')
  save(prey_details, file = 'Results/lv_Evo_Full_Prey_details.Rda')
  
  #Save predator results
  save(pred_res, file = 'Results/lv_Evo_Full_Pred.Rda')
  save(pred_details, file = 'Results/lv_Evo_Full_Pred_details.Rda')
  
  # Progress report
  print(G)
}

#Compile the results
load('Results/lv_Evo_Full_Prey.Rda')
load('Results/lv_Evo_Full_Prey_details.Rda')

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


