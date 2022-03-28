library(ctmm)


source("Scripts/Fit_Mods.R")


FILES <- list.files("~/Dropbox (Personal)/MultiSpecies_Data/Mammals",
                    pattern = ".Rda",
                    full.names = TRUE)


#For testing purposes, use a lion that runs through quickly
#i <- 51
#j <- 1


#Then walk through each of them sequentially
for(i in 6:length(FILES)){
  
  #Load in the tracking data
  DATA <- get(load(FILES[i]))
  
  #Store the species binomial
  BINOMIAL <- gsub("/Users/michaelnoonan/Dropbox (Personal)/MultiSpecies_Data/Mammals/", "", gsub(".Rda", "", FILES[i]), fixed = TRUE)
  
  
  message("Working on the tracking data on ", BINOMIAL)
  
  #Drop an outlier for the lion dataset
  if(BINOMIAL == "Panthera_leo"){DATA <- DATA[-which(DATA$individual.local.identifier == "Diana" & DATA$timestamp == "2005-10-12 00:00:47.000"),]}
  
  if(BINOMIAL == "Brachylagus_idahoensis"){
    
    DATA$gps.vdop <- 1/(DATA$gps.satellite.count-2)
    
    DATA$height.above.msl <- NULL
    
    #Convert to a telemetry object
    DATA <- as.telemetry(DATA)
    
    UERE <- uere.fit(DATA[1:12])
    
    DATA <- DATA[13:29]
    
    uere(DATA) <- UERE
    
  } else {
    
    
    #Convert to a telemetry object
    
    #if(i == 3 || i == 26 || i == 36){
    #  DATA <- as.telemetry(DATA[,c(3:5, 10)])
    #} else {
    
    #Convert to a telemetry object
    DATA <- as.telemetry(DATA)
    #}
  }
  
  
  #Create the paths to the file directory
  dir.create(paste("~/Dropbox (Personal)/UBC/Projects/BallisticMotion/Results/Empirical_Results/", BINOMIAL, sep = ""))
  dir.create(paste("~/Dropbox (Personal)/UBC/Projects/BallisticMotion/Results/Empirical_Results/", BINOMIAL, "/Fits", sep = ""))
  Model_path = paste("~/Dropbox (Personal)/UBC/Projects/BallisticMotion/Results/Empirical_Results/", BINOMIAL, "/Fits", sep = "")
  
  dir.create(paste("~/Dropbox (Personal)/UBC/Projects/BallisticMotion/Results/Empirical_Results/", BINOMIAL, "/Figs", sep = ""))
  Fig_Path = paste("~/Dropbox (Personal)/UBC/Projects/BallisticMotion/Results/Empirical_Results/", BINOMIAL, "/Figs", sep = "")
  
  
  #Then walk through each individual sequentially
  for(j in 1:length(DATA)){
    
    cat("Working on individual ", j, " of ", length(DATA), "\n")
    
    #Extract the current individual
    
    if(class(DATA) == "list") {cilla <- DATA[[j]]} else {
      
      cilla <- DATA
    }
    
    
    #Error handling
    error <- FALSE
    #Turn error on for the datasets with calibration data
    if(!is.infinite(uere(cilla)[[3]])){error <- TRUE}
    
    #Turn error on for the datasets with DOP values
    if("vx" %in% names(cilla)){error <- TRUE}
    if("HDOP" %in% names(cilla)){error <- TRUE}
    if("DOP" %in% names(cilla)){error <- TRUE}
    
    RESULTS <- tryCatch(
      
      {
        
        RESULTS <- CTMM_FIT(cilla,
                            Model_path = Model_path,
                            Fig_Path = Fig_Path,
                            error = error,
                            binomial = BINOMIAL)
        
        
      }, error=function(err) {
        message(cat("Model fitting failed, returning NAs \n"))
        
        RESULTS <- as.data.frame(t(unlist(c(BINOMIAL,
                                            cilla@info[1],
                                            rep(NA, 17)))))
        
        return(RESULTS)
      }
    )
    
    
    
    
    
    if(i ==1 && j == 1){
      
      write.table(RESULTS,
                  file = "~/Dropbox (Personal)/UBC/Projects/BallisticMotion/Results/lv_Results.csv",
                  row.names=FALSE,
                  col.names=TRUE,
                  sep=",",
                  append=TRUE)
      
    } else {
      
      write.table(RESULTS,
                  file = "~/Dropbox (Personal)/UBC/Projects/BallisticMotion/Results/lv_Results.csv",
                  row.names=FALSE,
                  col.names=FALSE,
                  sep=",",
                  append=TRUE)
      
    }
    
    
    
    
  }#Closes the loop that runs over the as.telemetry object
  
}#Closes the loop that runs over the individual datasets



