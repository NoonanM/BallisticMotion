#Set the working directory
setwd("~/Dropbox (Personal)/UBC/Projects/BallisticMotion")

#Import necessary packages
library(ggplot2)
library(gridExtra)

#Load in the results
load('Scripts/Sockeye_Simulations/BallisticMotion/Results/lv_Evo_45000g_Prey_Details_test.Rda')
load('Scripts/Sockeye_Simulations/BallisticMotion/Results/lv_Evo_45000g_Pred_Details_test.Rda')

load('Scripts/Sockeye_Simulations/BallisticMotion/Results/lv_Evo_45000g_Prey_test.Rda')
load('Scripts/Sockeye_Simulations/BallisticMotion/Results/lv_Evo_45000g_Pred_test.Rda')

# load('Results/lv_Evo_45000g_Prey_Details_test.Rda')
# load('Results/lv_Evo_45000g_Pred_Details_test.Rda')
# 
# load('Results/lv_Evo_45000g_Prey_test.Rda')
# load('Results/lv_Evo_45000g_Pred_test.Rda')


prey_res <- do.call(rbind, prey_res)
pred_res <- do.call(rbind, pred_res)
prey_details <- do.call(rbind, prey_details)
pred_details <- do.call(rbind, pred_details)

CIs <- data.frame(low = NA, est = NA, high = NA)
for(i in 1:nrow(prey_res)){
  
  
  CIs[i,] <- ctmm:::F.CI(mean(pred_details[which(pred_details$generation == i),"lv"]),
                         mean(pred_details[which(pred_details$generation == i),"lv"]),
                         mean(1/prey_details[which(prey_details$generation == i),"lv"]),
                         var(1/prey_details[which(prey_details$generation == i),"lv"]),
                         level=0.95)
}

CIs$generation <- prey_res$generation

c <- 
  ggplot(data=CIs) +
  ggtitle("c)") +
  #geom_hline(yintercept = median(RESULTS$ratio_est), linetype = "dashed", colour = "grey30", size = 0.3, alpha = 0.7) +
  geom_line(aes(y=est, x=generation), color = "purple", size = 0.2) +
  #scale_y_log10(breaks = c(0.01,0.1,1,10,100,1000,10000), labels = c(0,0.1,1,10,100,1000,10000)) +
  #scale_x_log10(breaks = c(0.01,0.1,1,10,100,1000,10000), labels = c(0,0.1,1,10,100,1000,10000)) +
  geom_ribbon(aes(ymin = low, ymax = high, x=generation), fill = "purple", alpha = 0.3) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "serif"),
        axis.title.x = element_text(size=8, family = "serif"),
        axis.text.y = element_text(size=6, family = "serif"),
        axis.text.x  = element_text(size=6, family = "serif"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "serif"),
        legend.position = "none") +
  ylab(expression(paste("Predator ", l[v], " : ", "Prey ", l[v]))) +
  xlab(expression(paste("Generation")))



#----------------------------------------------------------------------
# Prey fitness vs. lv
#----------------------------------------------------------------------

#Mean lv
prey_mu <- mean(prey_details[which(prey_details$generation %in% (nrow(prey_res)-50):nrow(prey_res)),"lv"])
#var lv
prey_sig <- var(prey_details[which(prey_details$generation %in% (nrow(prey_res)-50):nrow(prey_res)),"lv"])

prey_details$lv2 <- round(prey_details$lv,2)

AGG <- aggregate(offspring ~ lv2,
                 data = prey_details, FUN = "mean")

FIT <- nls(offspring ~ a * (lv2 + c) * exp(-b * (lv2 + c)),
           start = list(a = 60,
                        b = 0.1,
                        c = -1),
           data = AGG[which(AGG$lv2 < 9),])


ricker <- function(x) {
  coef(FIT)[1] * (x - coef(FIT)[3]) * exp(-coef(FIT)[2] * (x - coef(FIT)[3]))
}
x <- seq(0,150, 0.01)
y <- ricker(x)

d <- 
  ggplot() +
  ggtitle("d)") +
  geom_point(data=AGG, aes(y=offspring, x=lv2), col = "grey70", alpha = 0.8, size = 0.6, stroke = 0, shape=16) +
  geom_line(aes(y=y, x=x), color = "black", size = 0.4) +
  geom_rect(aes(xmin=prey_mu - prey_sig, xmax=prey_mu + prey_sig, ymin=-Inf, ymax=Inf), fill = "#046C9A", alpha = 0.3) +
  geom_vline(xintercept = prey_mu, linetype = "dashed", color = "#046C9A", size = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "serif"),
        axis.title.x = element_text(size=8, family = "serif"),
        axis.text.y = element_text(size=6, family = "serif"),
        axis.text.x  = element_text(size=6, family = "serif"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "serif"),
        legend.position = "none") +
  ylab("Prey fitness") +
  xlab("Ballistic length scale (m)") +
  #scale_x_continuous(limits = c(0,round(max(AGG$lv2)*1.1)), expand = c(0,0.1)) +
  scale_x_continuous(limits = c(0,30), expand = c(0,0.1)) +
  scale_y_continuous(limits = c(0,ceiling(max(AGG$offspring))), expand = c(0,0.5))

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

FIT <- nls(offspring ~ a * (lv2) * exp(-b * (lv2)),
           start = list(a = 0.5,
                        b = 0.02),
           data = AGG)


ricker <- function(x) {
  coef(FIT)[1] * (x) * exp(-coef(FIT)[2] * (x))
}
pred_x <- seq(0,500, 0.01)
pred_y <- ricker(pred_x)

e <- 
  ggplot() +
  ggtitle("e)") +
  geom_point(data=AGG, aes(y=offspring, x=lv2), col = "grey70", alpha = 0.6, size = 1.6, stroke = 0, shape=16) +
  geom_line(aes(y=pred_y, x=pred_x), color = "black", size = 0.4) +
  geom_rect(aes(xmin=pred_mu - pred_sig, xmax=pred_mu + pred_sig, ymin=-Inf, ymax=Inf), fill = "#FF0000", alpha = 0.3) +
  geom_vline(xintercept = pred_mu, linetype = "dashed", color = "#FF0000", size = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "serif"),
        axis.title.x = element_text(size=8, family = "serif"),
        axis.text.y = element_text(size=6, family = "serif"),
        axis.text.x  = element_text(size=6, family = "serif"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "serif"),
        legend.position = "none") +
  ylab("Predator fitness") +
  xlab("Ballistic length scale (m)") +
  #scale_x_continuous(limits = c(round(min(AGG$lv2)*0.9),round(max(AGG$lv2)*1.1)), expand = c(0,0.1)) +
  scale_x_continuous(limits = c(0,500), expand = c(0,0.1)) +
  scale_y_continuous(limits = c(0,ceiling(max(AGG$offspring))), expand = c(0,0.5))


BOT <- grid.arrange(c,d,e, ncol = 3)



a <- 
  ggplot() +
  ggtitle("a)") +
  geom_line(data=prey_res, aes(y=lv, x=generation), color = "#046C9A", size = 0.2) +
  geom_ribbon(data=prey_res, aes(ymin = lv - var, ymax = lv + var, x=generation), fill = "#046C9A", alpha = 0.3) + 
  # geom_point(data=prey_details,
  #            aes(y=lv, x=generation),
  #            color = "#046C9A",
  #            alpha = 0.6, size = 0.4, stroke = 0, shape=16) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "serif"),
        axis.title.x = element_text(size=8, family = "serif"),
        axis.text.y = element_text(size=6, family = "serif"),
        axis.text.x  = element_text(size=6, family = "serif"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "serif"),
        legend.position = "none") +
  ylab("Prey lv") + 
  xlab(expression(paste("Generation"))) +
  coord_cartesian(ylim = c(0,20), expand = c(0,0.5))
  #scale_y_continuous(limits = c(0,20), expand = c(0,0.5))

b <- 
ggplot() +
  ggtitle("b)") +

   geom_point(data=pred_details,
              aes(y=lv, x=generation),
              color = "grey70",
              alpha = 0.8, size = 0.4, stroke = 0, shape=16)  + 
  geom_line(data=pred_res, aes(y=lv, x=generation), color = "#FF0000", size = 0.4) +
  #geom_ribbon(data=pred_res, aes(ymin = lv - sqrt(var), ymax = lv + sqrt(var), x=generation), fill = "#FF0000", alpha = 0.3) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=8, family = "serif"),
        axis.title.x = element_text(size=8, family = "serif"),
        axis.text.y = element_text(size=6, family = "serif"),
        axis.text.x  = element_text(size=6, family = "serif"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "serif"),
        legend.position = "none") +
  ylab("Predator lv") + 
  xlab(expression(paste("Generation"))) +
  coord_cartesian(ylim = c(0,300), expand = c(0,0.5))
  #scale_y_continuous(limits = c(0,300), expand = c(0,0.5))


TOP  <- grid.arrange(a,b, ncol = 2)
BOT <- grid.arrange(c,d,e, ncol = 3)

grid.arrange(TOP,BOT, ncol = 1)
