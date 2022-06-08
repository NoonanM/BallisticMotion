#Set the working directory
setwd("~/Dropbox (Personal)/UBC/Projects/BallisticMotion")

#Import necessary packages
library(ggplot2)
library(gridExtra)
library(metafor)
library(rphylopic)
library(extraDistr)
library(cowplot)

#Source the functions
source("Scripts/Functions.R")


#Pull in the phylopics from the web (All are creative commons licensed)
wolf <- image_data("8cad2b22-30d3-4cbd-86a3-a6d2d004b201", size = 256)[[1]]
hare <- image_data("f69eb95b-3d0d-491d-9a7f-acddd419afed", size = 256)[[1]]
grass <- image_data("2af0a13e-69a8-4245-832e-ee3d981089b7", size = 256)[[1]]

#----------------------------------------------------------------------
# Predator and Prey lv across the mass spectrum (from simulations)
#----------------------------------------------------------------------

weights <- c(seq(5000, 100000, by = 5000))
setwd("Results/Sockeye_Simulations/Changing_Biomass/")

RESULTS <- list()
for(i in 1:length(weights)){
  FILES <- list.files(pattern = paste("_",format(weights[i],scientific = FALSE),"g", sep = ""))
  
  #Load the results
  for(j in 1:length(FILES)){load(FILES[j])}
  
  #Convert to data frames
  prey_res <- do.call(rbind, prey_res)
  pred_res <- do.call(rbind, pred_res)
  prey_details <- do.call(rbind, prey_details)
  pred_details <- do.call(rbind, pred_details)
  
  prey_lv <- prey_res[which(prey_res$generation %in% (nrow(prey_res)-50):nrow(prey_res)),"lv"]
  prey_var <- prey_res[which(prey_res$generation %in% (nrow(prey_res)-50):nrow(prey_res)),"var"]
  
  pred_lv <- pred_res[which(pred_res$generation %in% (nrow(pred_res)-50):nrow(pred_res)),"lv"]
  pred_var <- pred_res[which(pred_res$generation %in% (nrow(pred_res)-50):nrow(pred_res)),"var"]
  
  CIs <- ctmm:::F.CI(mean(pred_lv),
                     var(pred_details[which(pred_details$generation %in% (nrow(pred_res)-50):nrow(pred_res)),"lv"]),
                     mean(1/prey_lv),
                     var(1/prey_details[which(prey_details$generation %in% (nrow(prey_res)-50):nrow(prey_res)),"lv"]),
                     level=0.95)
  
  RESULTS[[i]]<-
    data.frame(pred_mass = weights[i], 
               prey_mass = prey.mass(weights[i]),
               prey_lv = mean(prey_lv),
               prey_var = var(prey_details[which(prey_details$generation %in% (nrow(prey_res)-50):nrow(prey_res)),"lv"]),
               
               pred_lv = mean(pred_lv),
               pred_var = var(pred_details[which(pred_details$generation %in% (nrow(pred_res)-50):nrow(pred_res)),"lv"]),
               
               ratio_est = CIs[2],
               ratio_low = CIs[1],
               ratio_high = CIs[3],
               
               
               ratio = mean(pred_lv/prey_lv),
               ratio_var = var(pred_lv/prey_lv))
  
  # 
  # #Generate figures to check for convergence
  # a <- 
  #   ggplot() +
  #   geom_point(data=prey_details,
  #               aes(y=lv, x=generation),
  #              color = "grey70",
  #               alpha = 0.6, size = 0.4, stroke = 0, shape=16) +
  #   geom_line(data=prey_res, aes(y=lv, x=generation), color = "#046C9A", size = 0.2) +
  #   theme_bw() +
  #   theme(panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         axis.title.y = element_text(size=8, family = "sans"),
  #         axis.title.x = element_text(size=8, family = "sans"),
  #         axis.text.y = element_text(size=6, family = "sans"),
  #         axis.text.x  = element_text(size=6, family = "sans"),
  #         plot.title = element_text(hjust = -0.05, size = 12, family = "sans"),
  #         legend.position = "none") +
  #   ylab("Prey lv") + 
  #   xlab(expression(paste("Generation"))) +
  #   scale_x_continuous(limits = c(0,300))
  # 
  # b <- 
  #   ggplot() +
  #   geom_point(data=pred_details,
  #              aes(y=lv, x=generation),
  #              color = "grey70",
  #              alpha = 0.8, size = 0.4, stroke = 0, shape=16)  + 
  #   geom_line(data=pred_res, aes(y=lv, x=generation), color = "#FF0000", size = 0.4) +
  #   theme_bw() +
  #   theme(panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(),
  #         axis.title.y = element_text(size=8, family = "sans"),
  #         axis.title.x = element_text(size=8, family = "sans"),
  #         axis.text.y = element_text(size=6, family = "sans"),
  #         axis.text.x  = element_text(size=6, family = "sans"),
  #         plot.title = element_text(hjust = -0.05, size = 12, family = "sans"),
  #         legend.position = "none") +
  #   ylab("Predator lv") + 
  #   xlab(expression(paste("Generation"))) +
  # scale_x_continuous(limits = c(0,300))
  # 
  # path <- paste("~/Dropbox (Personal)/UBC/Projects/BallisticMotion/Results/Sim_",weights[i],".png", sep = "")
  # 
  # Pair  <- grid.arrange(a,b, ncol = 2)
  # 
  # ggsave(Pair,
  #        width = 6.86, height = 3, units = "in",
  #        dpi = 600,
  #        bg = "transparent",
  #        file=path)
} 

RESULTS <- do.call(rbind, RESULTS)

#RESULTS <- RESULTS[-which(RESULTS$ratio_est > 30),] #Temporarily remove runs that didn't converge

#Get CIs on the mean ratio
test <- t.test(RESULTS$ratio_est)
summary(lm(ratio_est ~ pred_mass, data = RESULTS))
plot(ratio_est ~ pred_mass, data = RESULTS)
abline(lm(ratio_est ~ pred_mass, data = RESULTS))

#####################################################################


a <- 
  ggplot(data=RESULTS) +
  ggtitle("A") +
  geom_point(aes(x=prey_mass/1000, y=prey_lv, size = prey_var), color = "#3471bc", alpha = 0.7,stroke = 0,shape=16) +
  geom_point(aes(x=pred_mass/1000, y=pred_lv, size = pred_var), color = "#e6c141", alpha = 0.7,stroke = 0,shape=16)  +
  geom_smooth(aes(x=prey_mass/1000, y=prey_lv), method = "lm", color = "#3471bc", fill = "#3471bc", se = F, size = 0.5) +
  geom_smooth(aes(x=pred_mass/1000, y=pred_lv), method = "lm", color = "#e6c141", fill = "#e6c141", se = F, size = 0.5) +
  
  add_phylopic(hare, alpha = 1, x = 0.5, y = 0.5, ysize = 0.3, color = "#3471bc") +
  add_phylopic(wolf, alpha = 1, x = 1.1, y = 2.1, ysize = 0.3, color = "#e6c141") +
  
  ylab("Ballistic length scale (m)") +
  xlab("Body mass (kg)") +
  
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = "sans", face = "bold"),
        axis.title.x = element_text(size=10, family = "sans", face = "bold"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  
  scale_x_log10(breaks = c(0.01,0.1,1,10,100,1000,10000),
                labels = c(0,0.1,1,10,100,1000,10000),
                expand = c(0,0),
                limits = c(2,125)) +
  scale_y_log10(breaks = c(0.01,0.1,1,10,100,1000,10000),
                labels = c(0,0.1,1,10,100,1000,10000),
                limits = c(0.5,1000)) +
  annotation_logticks(outside = TRUE,
                      size = 0.3,
                      short = unit(0.05, "cm"),
                      mid = unit(0.05, "cm"),
                      long = unit(0.1, "cm")) +
  coord_cartesian(clip = "off")


b <- 
  ggplot(data=RESULTS) +
  ggtitle("B") +
  geom_rect(aes(xmin = 2, xmax = Inf,
                ymin = test$conf.int[1],
                ymax = test$conf.int[2]),
            alpha = 0.1,
            fill = "grey80") +
  geom_hline(yintercept = median(RESULTS$ratio_est), linetype = "dashed", colour = "grey30", size = 0.5, alpha = 0.7) +
  geom_segment(aes(x=pred_mass/1000, xend=pred_mass/1000,
                   y=ratio_low,
                   yend=ratio_high),
               color = "#3c7a47",
               alpha = 0.5,size = 0.3) +
  geom_point(aes(y=ratio_est, x=pred_mass/1000), color = "black", fill = NA, alpha = 1, shape = 21, stroke = 0.2, size = 2) +
  geom_point(aes(y=ratio_est, x=pred_mass/1000), color = "#3c7a47", alpha = 0.8, stroke = 0, shape=16, size = 2) +
  scale_x_log10(breaks = c(0.01,0.1,1,10,100,1000,10000),
                labels = c(0,0.1,1,10,100,1000,10000),
                expand = c(0,0),
                limits = c(2,125)) +
  scale_y_log10(breaks = c(0.01,0.1,1,10,100,1000,10000),
                labels = c(0,0.1,1,10,100,1000,10000),
                limits = c(0.5,280)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=12, family = "sans", face = "bold"),
        axis.title.x = element_text(size=10, family = "sans", face = "bold"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  #ylab(expression(paste("Predator ", l[v], " : ", "Prey ", l[v]))) +
  ylab("Predator:Prey") + 
  xlab("Predator mass (kg)") +
  annotation_logticks(outside = TRUE,
                      size = 0.3,
                      short = unit(0.05, "cm"),
                      mid = unit(0.05, "cm"),
                      long = unit(0.1, "cm")) +
  coord_cartesian(clip = "off")


TOP <- grid.arrange(a,b, ncol = 2)



#----------------------------------------------------------------------
# Diagram panel
#----------------------------------------------------------------------

setwd("~/Dropbox (Personal)/UBC/Projects/BallisticMotion")

#Load in the results for the 40kg predator simulations
load('Results/Sockeye_Simulations/Changing_Biomass/lv_Evo_40000g_Prey_Details_decline.Rda')
load('Results/Sockeye_Simulations/Changing_Biomass/lv_Evo_40000g_Pred_Details_decline.Rda')

load('Results/Sockeye_Simulations/Changing_Biomass/lv_Evo_40000g_Prey_decline.Rda')
load('Results/Sockeye_Simulations/Changing_Biomass/lv_Evo_40000g_Pred_decline.Rda')


prey_res <- do.call(rbind, prey_res)
pred_res <- do.call(rbind, pred_res)
prey_details <- do.call(rbind, prey_details)
pred_details <- do.call(rbind, pred_details)


#Mean lv
prey_mu <- mean(prey_details[which(prey_details$generation %in% (nrow(prey_res)-50):nrow(prey_res)),"lv"])
#var lv
prey_sig <- var(prey_details[which(prey_details$generation %in% (nrow(prey_res)-50):nrow(prey_res)),"lv"])

prey_details$lv2 <- round(prey_details$lv,2)

AGG <- aggregate(offspring ~ lv2,
                 data = prey_details, FUN = "mean")

FIT <- nls(offspring ~ a * (lv2 + c) * exp(-b * (lv2 + c)),
           start = list(a = 30,
                        b = 0.1,
                        c = -1),
           data = AGG)


ricker <- function(x) {
  coef(FIT)[1] * (x + coef(FIT)[3]) * exp(-coef(FIT)[2] * (x + coef(FIT)[3]))
}
x <- seq(0,150, 0.01)
y <- ricker(x)

c <- 
  ggplot() +
  ggtitle("C") +
  geom_point(data=AGG, aes(y=offspring, x=lv2), col = "grey70", alpha = 0.6, size = 0.6, stroke = 0, shape=16) +
  geom_line(aes(y=y, x=x), color = "black", size = 0.2) +
  geom_rect(aes(xmin=prey_mu - prey_sig, xmax=prey_mu + prey_sig, ymin=-Inf, ymax=Inf), fill = "#3471bc", alpha = 0.3) +
  geom_vline(xintercept = prey_mu, linetype = "dashed", color = "#3471bc", size = 0.3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans", face = "bold"),
        axis.title.x = element_text(size=8, family = "sans", face = "bold"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  ylab("Prey fitness") +
  xlab(expression(paste(l[v]))) + 
  scale_x_continuous(limits = c(0,20), expand = c(0,0.1)) + #c(0,round(max(AGG$lv2)*1.1))
  scale_y_continuous(limits = c(0,40), expand = c(0,0.5))


#----------------------------------------------------------------------
# Predator fitness vs. lv
#----------------------------------------------------------------------

#Mean lv
pred_mu <- mean(pred_details[which(pred_details$generation %in% (nrow(pred_res)-50):nrow(pred_res)),"lv"])
#var lv
pred_sig <- var(pred_details[which(pred_details$generation %in% (nrow(pred_res)-50):nrow(pred_res)),"lv"])

pred_details$lv2 <- round(pred_details$lv,0)

AGG <- aggregate(offspring ~ lv2,
                 data = pred_details, FUN = "mean")

FIT <- nls(offspring ~ a * (lv2 ) * exp(-b * (lv2 )),
           start = list(a = 0.5,
                        b = 0.01),
           data = AGG)


ricker <- function(x) {
  coef(FIT)[1] * (x ) * exp(-coef(FIT)[2] * (x ))
}
pred_x <- seq(0,450, 0.01)
pred_y <- ricker(pred_x)

d <- 
  ggplot() +
  ggtitle("D") +
  geom_point(data=AGG, aes(y=offspring, x=lv2), col = "grey70", alpha = 0.6, size = 0.6, stroke = 0, shape=16) +
  geom_line(aes(y=pred_y, x=pred_x), color = "black", size = 0.2) +
  geom_rect(aes(xmin=pred_mu - pred_sig, xmax=pred_mu + pred_sig, ymin=-Inf, ymax=Inf), fill = "#e6c141", alpha = 0.3) +
  geom_vline(xintercept = pred_mu, linetype = "dashed", color = "#e6c141", size = 0.3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans", face = "bold"),
        axis.title.x = element_text(size=8, family = "sans", face = "bold"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  ylab("Predator fitness") +
  xlab(expression(paste(l[v]))) + 
  scale_x_continuous(limits = c(0,400), expand = c(0,0.1)) + #c(round(min(AGG$lv2)*0.9),round(max(AGG$lv2)*1.1))
  scale_y_continuous(limits = c(0,37), expand = c(0,0.5))


# c <-
#   ggplot(data=AGG) +
#   ggtitle("C") +
#   theme_bw() +
#   
#   add_phylopic(grass, alpha = 1, x = 0.05, y = 0.55, ysize = 0.201) +
#   add_phylopic(grass, alpha = 1, x = 0.05, y = 0.55, ysize = 0.2, color = "#3c7a47") +
#   add_phylopic(grass, alpha = 1, x = 0.18, y = 0.45, ysize = 0.201) +
#   add_phylopic(grass, alpha = 1, x = 0.18, y = 0.45, ysize = 0.2, color = "#3c7a47") +
#   add_phylopic(grass, alpha = 1, x = 0.07, y = 0.3, ysize = 0.201) +
#   add_phylopic(grass, alpha = 1, x = 0.07, y = 0.3, ysize = 0.2, color = "#3c7a47") +
#   
#   add_phylopic(hare, alpha = 1, x = 0.5, y = 0.5, ysize = 0.302) +
#   add_phylopic(hare, alpha = 1, x = 0.5, y = 0.5, ysize = 0.3, color = "#3471bc") +
#   
#   add_phylopic(wolf, alpha = 1, x = 0.9, y = 0.5, ysize = 0.302) +
#   add_phylopic(wolf, alpha = 1, x = 0.9, y = 0.5, ysize = 0.3, color = "#e6c141") +
#   
#   geom_text(x=0.7, y=0.87, label="Diffusive", family = "sans", size = 2.5) +
#   geom_text(x=0.27, y=0.87, label="Ballistic", family = "sans", size = 2.5) +
#   geom_text(x=0.68, y=0.58, label="Ballistic", family = "sans", size = 2.5) +
#   geom_curve(aes(x = 0.4,
#                  y = 0.7,
#                  xend = 0.15,
#                  yend = 0.7),
#              color = "#3471bc",
#              arrow = arrow(length = unit(0.03, "npc"))) +
#   geom_curve(aes(x = 0.55,
#                  y = 0.7,
#                  xend = 0.85,
#                  yend = 0.7),
#              curvature = -0.5,
#              color = "#3471bc",
#              arrow = arrow(length = unit(0.03, "npc"))) +
#   geom_curve(aes(x = 0.8,
#                  y = 0.6,
#                  xend = 0.55,
#                  yend = 0.6),
#              curvature = 0.5,
#              color = "#e6c141",
#              arrow = arrow(length = unit(0.03, "npc"))) +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.title.y = element_blank(),
#         axis.title.x = element_blank(),
#         axis.text.y = element_blank(),
#         axis.text.x  = element_blank(),
#         axis.ticks = element_blank(),
#         plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
#         legend.position = "none",
#         panel.background = element_rect(fill = "transparent"),
#         plot.background = element_rect(fill = "transparent", color = NA)) +
#   scale_y_continuous(limits = c(0,1)) +
#   scale_x_continuous(limits = c(0,1))
# 

# c <-
#   ggdraw() +
#   draw_plot(c) +
#   draw_plot(inset_1, x = 0.2, y = .05, width = .3, height = .3) +
#   draw_plot(inset_2, x = 0.55, y = .05, width = .3, height = .3)


#----------------------------------------------------------------------
# Example of the simulated movement data
#----------------------------------------------------------------------

# 
# #Predator mass (g)
# mass_pred <- 40000
# 
# #Prey mass (g)
# mass_prey <- prey.mass(mass_pred)
# 
# #"Lifespan" and sampling interval for the simulations
# t <- sampling(mass_prey, crossings = 30)
# 
# 
# #Generate the prey movement models
# #Generate the HR centres of the prey
# CENTRES <- rbvpois(n = 10,
#                    a = pred.SIG(mass_pred)*.75,
#                    b = pred.SIG(mass_pred)*.75,
#                    c = 0)
# CENTRES <- scale(CENTRES, scale = FALSE)
# 
# PREY_mods <- list()
# for(i in 1:10){
#   # Prey movement parameters
#   prey_tau_p <- mean(prey_details[which(prey_details$generation == max(prey_details$generation)),"tau_p"])
#   prey_tau_v <- mean(prey_details[which(prey_details$generation == max(prey_details$generation)),"tau_v"])
#   prey_sig <- mean(prey_details[which(prey_details$generation == max(prey_details$generation)),"sig"])
#   
#   PREY_mods[[i]] <- ctmm(tau = c(prey_tau_p,prey_tau_v),
#                          mu = c(CENTRES[i,1],CENTRES[i,2]),
#                          sigma = prey_sig)
# } #Closes loop over n_prey
# 
# 
# # Predator movement parameters
# pred_tau_p <- mean(pred_details[which(pred_details$generation == max(pred_details$generation)),"tau_p"])
# pred_tau_v <- mean(pred_details[which(pred_details$generation == max(pred_details$generation)),"tau_v"])
# pred_sig <- mean(pred_details[which(pred_details$generation == max(pred_details$generation)),"sig"])
# 
# PRED_mod <- ctmm(tau = c(pred_tau_p,
#                          pred_tau_v),
#                  mu = c(0,0),
#                  sigma = pred_sig)
# 
# #Simulate the prey movement
# #Parallelised to speed up run times
# PREY_tracks <- lapply(PREY_mods,
#                       FUN = simulate,
#                       t = t)
# 
# #Simulate the predator movement
# PRED_tracks <- simulate(PRED_mod,t = t)
# 
# COLS <- viridis::viridis(10)
# 
# d <- 
#   ggplot() +
#   ggtitle("D") +
#   # geom_path(aes(y=PREY_tracks[[1]]$y, x=PREY_tracks[[1]]$x), color =COLS[1], size = 0.1, alpha = 0.3) +
#   # geom_path(aes(y=PREY_tracks[[2]]$y, x=PREY_tracks[[2]]$x), color =COLS[2], size = 0.1, alpha = 0.3) +
#   # geom_path(aes(y=PREY_tracks[[3]]$y, x=PREY_tracks[[3]]$x), color =COLS[3], size = 0.1, alpha = 0.3) +
#   # geom_path(aes(y=PREY_tracks[[4]]$y, x=PREY_tracks[[4]]$x), color =COLS[4], size = 0.1, alpha = 0.3) +
#   # geom_path(aes(y=PREY_tracks[[5]]$y, x=PREY_tracks[[5]]$x), color =COLS[5], size = 0.1, alpha = 0.3) +
#   # geom_path(aes(y=PREY_tracks[[6]]$y, x=PREY_tracks[[6]]$x), color =COLS[6], size = 0.1, alpha = 0.3) +
#   # geom_path(aes(y=PREY_tracks[[7]]$y, x=PREY_tracks[[7]]$x), color =COLS[7], size = 0.1, alpha = 0.3) +
#   # geom_path(aes(y=PREY_tracks[[8]]$y, x=PREY_tracks[[8]]$x), color =COLS[8], size = 0.1, alpha = 0.3) +
#   # geom_path(aes(y=PREY_tracks[[9]]$y, x=PREY_tracks[[9]]$x), color =COLS[9], size = 0.1, alpha = 0.3) +
#   # geom_path(aes(y=PREY_tracks[[10]]$y, x=PREY_tracks[[10]]$x), color =COLS[10], size = 0.1, alpha = 0.3) +
#   geom_path(aes(y=PREY_tracks[[1]]$y, x=PREY_tracks[[1]]$x), color = "#3471bc", size = 0.1, alpha = 0.2) +
#   geom_path(aes(y=PREY_tracks[[2]]$y, x=PREY_tracks[[2]]$x), color ="#3471bc", size = 0.1, alpha = 0.2) +
#   geom_path(aes(y=PREY_tracks[[3]]$y, x=PREY_tracks[[3]]$x), color ="#3471bc", size = 0.1, alpha = 0.2) +
#   geom_path(aes(y=PREY_tracks[[4]]$y, x=PREY_tracks[[4]]$x), color ="#3471bc", size = 0.1, alpha = 0.2) +
#   geom_path(aes(y=PREY_tracks[[5]]$y, x=PREY_tracks[[5]]$x), color ="#3471bc", size = 0.1, alpha = 0.2) +
#   geom_path(aes(y=PREY_tracks[[6]]$y, x=PREY_tracks[[6]]$x), color ="#3471bc", size = 0.1, alpha = 0.2) +
#   geom_path(aes(y=PREY_tracks[[7]]$y, x=PREY_tracks[[7]]$x), color ="#3471bc", size = 0.1, alpha = 0.2) +
#   geom_path(aes(y=PREY_tracks[[8]]$y, x=PREY_tracks[[8]]$x), color ="#3471bc", size = 0.1, alpha = 0.2) +
#   geom_path(aes(y=PREY_tracks[[9]]$y, x=PREY_tracks[[9]]$x), color ="#3471bc", size = 0.1, alpha = 0.2) +
#   geom_path(aes(y=PREY_tracks[[10]]$y, x=PREY_tracks[[10]]$x), color ="#3471bc", size = 0.1, alpha = 0.3) +
#   geom_path(aes(y=PRED_tracks$y, x=PRED_tracks$x), color = "#e6c141", size = 0.1, alpha = 0.8) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title.y = element_text(size=8, family = "sans", face = "bold"),
#         axis.title.x = element_text(size=8, family = "sans", face = "bold"),
#         axis.text.y = element_text(size=6, family = "sans"),
#         axis.text.x  = element_text(size=6, family = "sans"),
#         plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
#         legend.position = "none",
#         panel.background = element_rect(fill = "transparent"),
#         plot.background = element_rect(fill = "transparent", color = NA)) +
#   ylab("Y (m)") +
#   xlab("X (m)")
# 
# 

#----------------------------------------------------------------------
# Example of the results for a single run
#----------------------------------------------------------------------


CIs <- data.frame(low = NA, est = NA, high = NA, standard = NA)
for(i in 1:nrow(prey_res)){
  
  CIs[i,1:3] <- ctmm:::F.CI(mean(pred_details[which(pred_details$generation == i),"lv"]),
                            var(pred_details[which(pred_details$generation == i),"lv"]),
                            mean(1/prey_details[which(prey_details$generation == i),"lv"]),
                            var(1/prey_details[which(prey_details$generation == i),"lv"]),
                            level=0.5)
  
  CIs[i,4] <- mean(pred_details[which(pred_details$generation == i),"lv"])/mean(prey_details[which(prey_details$generation == i),"lv"])
}

CIs$generation <- prey_res$generation

CIs <- CIs[which(CIs$generation < 301),]

e <- 
  ggplot(data=CIs) +
  ggtitle("E") +
  geom_hline(yintercept = median(RESULTS$ratio_est), linetype = "dashed", colour = "grey30", size = 0.3, alpha = 0.7) +
  geom_line(aes(y=est, x=generation), color = "#3c7a47", size = 0.2) +
  #scale_y_log10(breaks = c(0.01,0.1,1,10,100,1000,10000), labels = c(0,0.1,1,10,100,1000,10000)) +
  #scale_x_log10(breaks = c(0.01,0.1,1,10,100,1000,10000), labels = c(0,0.1,1,10,100,1000,10000)) +
  geom_ribbon(aes(ymin = low, ymax = high, x=generation), fill = "#3c7a47", alpha = 0.3) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "sans", face = "bold"),
        axis.title.x = element_text(size=8, family = "sans", face = "bold"),
        axis.text.y = element_text(size=6, family = "sans"),
        axis.text.x  = element_text(size=6, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "none",
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  #ylab(expression(paste("Predator ", l[v], " : ", "Prey ", l[v]))) +
  ylab("Predator:Prey") + 
  xlab("Generation") +
  scale_x_continuous(expand = c(0,2))


BOT <- grid.arrange(c,d,e, ncol = 3)

plots <- grid.arrange(TOP, BOT, ncol = 1, heights=c(1.5,1))


#Save the figures
ggsave(plots,
       width = 6.86, height = 4.5, units = "in",
       dpi = 600,
       bg = "transparent",
       file="Results/lv_Scaling_Simulations_ChangingBiomass.png")






