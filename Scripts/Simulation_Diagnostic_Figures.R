#load in the necessary packages
library(ggplot2)
library(gridExtra)

#Source the functions
source("Scripts/Functions.R")


#Load in the results for the 40kg predator simulations
load('Results/Sockeye_Simulations/lv_Evo_40000g_Prey_Details.Rda')
load('Results/Sockeye_Simulations/lv_Evo_40000g_Pred_Details.Rda')

load('Results/Sockeye_Simulations/lv_Evo_40000g_Prey.Rda')
load('Results/Sockeye_Simulations/lv_Evo_40000g_Pred.Rda')


prey_res <- do.call(rbind, prey_res)
pred_res <- do.call(rbind, pred_res)
prey_details <- do.call(rbind, prey_details)
pred_details <- do.call(rbind, pred_details)


#Calculate prey mass
prey_mass <- prey.mass(40000)


#----------------------------------------------------------------------
# Figure depicting movement speed versus ballistic length scale
#----------------------------------------------------------------------
a <- 
  ggplot(data=prey_details) +
  ggtitle("A") +
  geom_point(aes(y=speed, x=lv), color = "#3471bc", size = 0.2, alpha = 0.1, pch = 16) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans", face = "bold"),
        axis.title.x = element_text(size=8, family = "sans", face = "bold"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none") +
  xlab(expression(paste("Prey ", l[v]))) +
  ylab("Movement speed (m/s)")


#----------------------------------------------------------------------
# Figure depicting energetic costs versus ballistic length scale
#----------------------------------------------------------------------

#First calculate E_total
# Field metabolic rate (in kj/day) from Nagy 1987 https://doi.org/10.2307/1942620 
FMR <- 10^(0.774 + 0.727*log10(prey_mass))

# total lifespan in days (based on 30 range crossings)
lifespan <- round(prey.tau_p(prey_mass)*30) /60/60/24

# Metabolic cost of movement in watts/kg from Taylor et al. 1982 https://doi.org/10.1242/jeb.97.1.1 
E = 10.7*(prey_mass/1000)^(-0.316)*prey_details$speed + 6.03*(prey_mass/1000)^(-0.303)

#Convert to kJ/s
E <- (E * (prey_mass/1000))/1000

#Total energetic cost in kj as a function of FMR and movement speed
prey_details$E_total <- FMR * lifespan + E*prey.tau_p(prey_mass)*30

b <- 
  ggplot(data=prey_details) +
  ggtitle("B") +
  geom_point(aes(y=E_total, x=lv), color = "#3471bc", size = 0.2, alpha = 0.1, pch = 16) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans", face = "bold"),
        axis.title.x = element_text(size=8, family = "sans", face = "bold"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none") +
  xlab(expression(paste("Prey ", l[v]))) +
  ylab(expression(paste(E[total]," (kj)")))




#----------------------------------------------------------------------
# Figure depicting number of encounters versus ballistic length scale
#----------------------------------------------------------------------


c <- 
  ggplot(data=prey_details) +
  ggtitle("C") +
  geom_point(aes(y=patches, x=lv), color = "#3471bc", size = 0.2, alpha = 0.1, pch = 16) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans", face = "bold"),
        axis.title.x = element_text(size=8, family = "sans", face = "bold"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none") +
  xlab(expression(paste("Prey ", l[v]))) +
  ylab("Number of patches")



#----------------------------------------------------------------------
# Figure depicting number of encounters versus ballistic length scale
#----------------------------------------------------------------------


d <- 
  ggplot(data=prey_details) +
  ggtitle("D") +
  geom_point(aes(y=patches, x=speed), color = "#3471bc", size = 0.2, alpha = 0.1, pch = 16) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans", face = "bold"),
        axis.title.x = element_text(size=8, family = "sans", face = "bold"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none") +
  xlab("Movement speed (m/s)") +
  ylab("Number of patches") +
  scale_x_continuous(limits = c(0,5.8))



#----------------------------------------------------------------------
# Compile and save
#----------------------------------------------------------------------

plots <- grid.arrange(a,b,
                      c,d,
                      #top = "14kg Prey",
                      ncol = 2)

ggsave(plots,
       width = 6.86, height = 6, units = "in",
       dpi = 600,
       bg = "transparent",
       file="Results/Prey_Diagnostics.png")





#----------------------------------------------------------------------
# Repeat for the predator
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Figure depicting movement speed versus ballistic length scale
#----------------------------------------------------------------------
a <- 
  ggplot(data=pred_details) +
  ggtitle("A") +
  geom_point(aes(y=speed, x=lv), color = "#e6c141", size = 0.2, alpha = 0.1, pch = 16) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans", face = "bold"),
        axis.title.x = element_text(size=8, family = "sans", face = "bold"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none") +
  xlab(expression(paste("Predator ", l[v]))) +
  ylab("Movement speed (m/s)")


#----------------------------------------------------------------------
# Figure depicting energetic costs versus ballistic length scale
#----------------------------------------------------------------------

#First calculate E_total
# Field metabolic rate (in kj/day) from Nagy 1987 https://doi.org/10.2307/1942620 
FMR <- 10^(0.774 + 0.727*log10(40000))

# total lifespan in days (based on 30 range crossings)
lifespan <- round(prey.tau_p(prey_mass)*30) /60/60/24

# Metabolic cost of movement in watts/kg from Taylor et al. 1982 https://doi.org/10.1242/jeb.97.1.1 
E = 10.7*(40000/1000)^(-0.316)*pred_details$speed + 6.03*(40000/1000)^(-0.303)

#Convert to kJ/s
E <- (E * (40000/1000))/1000

#Total energetic cost in kj as a function of FMR and movement speed
pred_details$E_total <- FMR * lifespan + E*prey.tau_p(prey_mass)*30

b <- 
  ggplot(data=pred_details) +
  ggtitle("B") +
  geom_point(aes(y=E_total, x=lv), color = "#e6c141", size = 0.2, alpha = 0.1, pch = 16) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans", face = "bold"),
        axis.title.x = element_text(size=8, family = "sans", face = "bold"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none") +
  xlab(expression(paste("Predator ", l[v]))) +
  ylab(expression(paste(E[total]," (kj)")))




#----------------------------------------------------------------------
# Figure depicting number of encounters versus ballistic length scale
#----------------------------------------------------------------------

pred_details$lv2 <- round(pred_details$lv,1)

AGG <- aggregate(encounter ~ lv2,
                 data = pred_details, FUN = "mean")

c <- 
  ggplot(data=AGG) +
  ggtitle("C") +
  geom_point(aes(y=encounter, x=lv2), color = "#e6c141", size = 1, alpha = 0.3, pch = 16) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans", face = "bold"),
        axis.title.x = element_text(size=8, family = "sans", face = "bold"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none") +
  xlab(expression(paste("Predator ", l[v]))) +
  ylab("Number of prey encounters") #+
  #scale_x_log10(breaks = c(0.01,0.1,1,10,100,1000,10000), labels = c(0,0.1,1,10,100,1000,10000), expand = c(0,0.05))



#----------------------------------------------------------------------
# Figure depicting number of encounters versus ballistic length scale
#----------------------------------------------------------------------

pred_details$speed2 <- round(pred_details$speed,4)

AGG <- aggregate(encounter ~ speed2,
                 data = pred_details, FUN = "mean")

d <- 
  ggplot(data=AGG) +
  ggtitle("D") +
  geom_point(aes(y=encounter, x=speed2), color = "#e6c141", size = 1, alpha = 0.3, pch = 16) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans", face = "bold"),
        axis.title.x = element_text(size=8, family = "sans", face = "bold"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none") +
  xlab(expression(paste("Movement speed (m/s)"))) +
  ylab("Number of prey encounters")



#----------------------------------------------------------------------
# Compile and save
#----------------------------------------------------------------------

plots <- grid.arrange(a,b,
                      c,d,
                      #top = "40kg Predator",
                      ncol = 2)

ggsave(plots,
       width = 6.86, height = 6, units = "in",
       dpi = 600,
       bg = "transparent",
       file="Results/Predator_Diagnostics.png")
