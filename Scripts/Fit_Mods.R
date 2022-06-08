
# This script is a function for fitting movement models to relocation data
# It uses the best fit model to then estimate a an AKDE HR
# The fitted model, and UD are then exported for future use, and the area estimate +/- CIs is returned
# Finally, a 4 panel figure is returned for diagnostic purposes


#It requires as inputs: i) tracking data on 1 individual, as a telemetry object
#                       ii) A path to the location where the movement models will be saved
#                       iv) A path to the location where the figures will be saved



#Written by Michael Noonan

#Last updated: Mar 26th 2022



library(ctmm)

#This function is used to fit movement models to gps data and calculate an akde home range estimate
#It requires ctmm

CTMM_FIT <- function(cilla,
                     Model_path,
                     Fig_Path,
                     error = FALSE,
                     binomial = NA){
  
  message("Fitting movement model for: ", unlist(cilla@info[1]))
  
  #First fit the movement model using the usual workflow
  #Generate the variogram
  vg.cilla <- variogram(cilla)
  
  #fit the variogram
  GUESS <- ctmm.guess(cilla, interactive= FALSE)
  
  if(error){GUESS$error <- TRUE}
  
  
  #Fit the movement models
  cilla.mods <- tryCatch(
    {
      cilla.mods <- ctmm.select(cilla,
                                CTMM = GUESS,
                                cores = -1,
                                trace = TRUE)
    }, error=function(err) {
      message(paste("Model fitting using pHREML and optim failed, defaulting to ML"))
      
      #If pREML and/or the optimiser fail, this will try it again
      #using the standard ML without any optimisation
      
      cilla.mods <- ctmm.select(cilla,
                                CTMM = GUESS,
                                method = "ML")
      
      #The return here is essential, otherwise results aren't returned and it fails
      return(cilla.mods)
    }
  )
  
  #Save the best fit model
  #First get the best fit model if more than 1 have been fit
  if(class(cilla.mods) == "list") {FIT <- cilla.mods[[1]]} else {
    
    FIT <- cilla.mods
  }
  
  #Assign the file path
  ctmm.path <- file.path(Model_path, paste("Fits_", cilla@info[1], ".rda", sep = ""))
  
  
  #And save
  save(FIT, file = ctmm.path)
  
  
  
  
  
  
  ############################################
  #Plot all the results and save them as a pdf
  ############################################
  
  message("\n", "Saving the figures")
  
  #Assign the file path and name for saving the results
  fig.path <- file.path(Fig_Path,
                        paste("ctmm_", cilla@info[1], ".png", sep = ""))
  
  #Save the graphic device's output as a pdf
  png(file=fig.path,
      type="cairo",
      units = "in",
      width = 6.81, height = 8,
      pointsize = 10,
      res = 600) #
  
  #Set the par to plot all on same screen
  par(mfrow=c(2,2),
      mgp = c(1.5, 0.5, 0),
      oma=c(0,0,0,0),
      mar=c(3,3,2,2),
      cex.lab=1.5,
      cex.main = 2) 
  
  
  
  #Plot the zoomed in variogram 
  tryCatch(
    {
      plot(vg.cilla, CTMM=FIT,family="serif", fraction = 0.01) 
      title(main = "a)", family = "serif", adj = 0)
      
    }, error=function(err) {
      
      plot(vg.cilla, CTMM=FIT,family="serif", fraction = 0.05) 
      title(main = "a)", family = "serif",  adj = 0)
      
    }
    
  )
  
  
  
  #Plot the variogram of the full time series
  plot(vg.cilla, CTMM=FIT,family="serif") 
  title(main = "b)", family = "serif",  adj = 0)
  
  
  #Plot a check for outliers
  OUTLIERS <- outlie(cilla, family="serif")
  title(main = "c)", family = "serif", adj = 0)
  
  
  #plot(OUTLIERS, family="serif")
  #title(main = "d)", family = "serif",  adj = 0)
  
  
  #Plot the AKDE range estimate, with the relocation data, coloured by time
  #Create a function that scales colours between red and blue
  rbPal <- colorRampPalette(c('#FF0000','#046C9A'))
  #Then create a variable that scales from red to blue between the two times
  cilla$Col <- rbPal(nrow(cilla))[as.numeric(cut(cilla$t,breaks = nrow(cilla)))]
  
  
  #Plot of the range estimate
  plot(cilla,
       col.grid = NA,
       family = "serif",
       pch = 20,
       cex = 0.5,
       col = cilla$Col,
       labels=FALSE,
       error = F)
  
  title(main = "d)", family = "serif",  adj = 0)
  
  
  
  dev.off()
  
  
  
  ##########################################################################################
  ##########################################################################################
  # Export all of the results
  ##########################################################################################
  ##########################################################################################
  
  #Get basic stats on the dataset
  res <- data.frame(binomial = binomial)
  res$ID <- cilla@info$identity
  res$Year <- median(as.numeric(format(as.Date(cilla$timestamp),"%Y")))
  res$Lat <- median(cilla$latitude)
  res$Long <- median(cilla$longitude)
  res$Frequency <- median(diff(cilla$t))/60
  res$Duration <- (tail(cilla$t,1) - head(cilla$t,1))/60/60/24
  res$n <- nrow(cilla)
  res$n_area <- summary(FIT)$DOF["area"]
  
  #Get tau_p
  res$tau_p <- if(nrow(summary(FIT, units = FALSE)$CI) >= 2){summary(FIT, units = FALSE)$CI[2,2]} else {NA}
  res$tau_p_min <- if(nrow(summary(FIT, units = FALSE)$CI) >= 2){summary(FIT, units = FALSE)$CI[2,1]} else {NA}
  res$tau_p_max <- if(nrow(summary(FIT, units = FALSE)$CI) >= 2){summary(FIT, units = FALSE)$CI[2,3]} else {NA}
  if(grepl("OUf", summary(FIT)$name, fixed = TRUE)){res$tau_p_var <- FIT$COV["tau","tau"]} else if(nrow(summary(FIT, units = FALSE)$CI) >= 2){
    res$tau_p_var <- FIT$COV["tau position","tau position"]} else {NA}
  
  
  
  #Get tau_v
  res$tau_v <- if(grepl("OUF", summary(FIT)$name, fixed = TRUE)){summary(FIT, units = FALSE)$CI[3,2]} else {NA}
  res$tau_v_min <- if(grepl("OUF", summary(FIT)$name, fixed = TRUE)){summary(FIT, units = FALSE)$CI[3,1]} else {NA}
  res$tau_v_max <- if(grepl("OUF", summary(FIT)$name, fixed = TRUE)){summary(FIT, units = FALSE)$CI[3,3]} else {NA}
  res$tau_v_var <- if(grepl("OUF", summary(FIT)$name, fixed = TRUE)){FIT$COV["tau velocity","tau velocity"]} else {NA}
  
  
  #Get spatial variance (sigma)
  res$sigma <- ctmm:::area.covm(FIT$sigma)
  
  #Get ballistic length scale
  res$l_v <- sqrt((res$tau_v/res$tau_p)*res$sigma)
  
  #Return all of the results
  return(res)
  
}

