#Set the working directory
setwd("~/Dropbox (Personal)/UBC/Projects/BallisticMotion")


# Load in the prey results
load('Results/lv_Evo_Full_Prey.Rda')
load('Results/lv_Evo_Full_Prey_details.Rda')

prey_res <- do.call(rbind, prey_res)
prey_details <- do.call(rbind, prey_details)


# Load in the predator results
load('Results/lv_Evo_Full_Pred.Rda')
load('Results/lv_Evo_Full_Pred_details.Rda')

pred_res <- do.call(rbind, pred_res)
pred_details <- do.call(rbind, pred_details)

ratio <- pred_res$lv/prey_res$lv

png(filename="Results/Simulation_Results.png",
    width = 6.86, height = 6.2, units = "in",
    res = 600)
par(mar = c(4,
            4,
            0.5,
            0.7),
    mgp=c(2,0.6,0),
    family = "serif")   
layout(matrix(c(1,1,
                2,3,
                4,5), ncol = 2, byrow = TRUE))

lv_pred <- pred_details[which(pred_details$generation %in% (nrow(pred_res)-100):(nrow(pred_res)-100)),"lv"]
lv_prey <- prey_details[which(prey_details$generation %in% (nrow(prey_res)-100):(nrow(prey_res)-100)),"lv"]
mu <- median(lv_pred/lv_prey)

#par(mfrow = c(1,2))
plot(ratio ~ generation,
     data = prey_res,
     type = "l",
     col = "purple",
     ylim = c(2,52),
     xlim = c(60,1950),
     ylab = expression(paste(Predator[l[v]], ":", Prey[l[v]])),
     xlab = "Generation")
abline(h = mu, col = 'grey30', lty = 'dashed')
#abline(h = 1, col = 'grey30', lty = 'dashed') #Add empirical ratio?



#Compare prey lv with starting values
PREY_LV <- prey_res$lv[1] #(sqrt((tau_v/tau_p)*sig))

#par(mfrow = c(1,2))
plot(lv/PREY_LV ~ generation,
     data = prey_res,
     type = "l",
     col = "#046C9A",
     ylim = c(0,1.2),
     xlim = c(60,1950),
     ylab = expression(paste("Relative change in ", Prey[l[v]])),
     xlab = "Generation")
abline(h = 1, col = 'grey30', lty = 'dashed')

polygon(c(prey_res$generation, rev(prey_res$generation)), c((prey_res$lv - prey_res$var)/PREY_LV, rev((prey_res$lv + prey_res$var)/PREY_LV)),
        col = adjustcolor("#046C9A",alpha.f=0.5), border=NA)



#Compare prey lv with starting values
PRED_LV <- pred_res$lv[1] #(sqrt((tau_v/tau_p)*sig))

#par(mfrow = c(1,2))
plot(lv/PRED_LV ~ generation,
     data = pred_res,
     type = "l",
     col = "red",
     ylim = c(0,1.2),
     xlim = c(60,1950),
     ylab = expression(paste("Relative change in ", Predator[l[v]])),
     xlab = "Generation")
abline(h = 1, col = 'grey30', lty = 'dashed')

polygon(c(pred_res$generation, rev(pred_res$generation)), c((pred_res$lv - pred_res$var)/PRED_LV, rev((pred_res$lv + pred_res$var)/PRED_LV)),
        col = adjustcolor("red",alpha.f=0.5), border=NA)

#Prey l_v and fitness
DATA <- prey_details[which(prey_details$offspring != 0),]

prey_details$lv2 <- round(prey_details$lv, 2)
AGG <- aggregate(offspring ~ lv2,
                 data = prey_details, FUN = "mean")

FIT <- nls(offspring ~ a * (lv2 + c) * exp(-b * (lv2 + c)),
           start = list(a = 6,
                        b = 0.1,
                        c = -2),
           data = AGG)


ricker <- function(x) {
  coef(FIT)[1] * (x + coef(FIT)[3]) * exp(-coef(FIT)[2] * (x + coef(FIT)[3]))
}
x <- seq(0,130, 0.01)
y <- ricker(x)

#Mean lv
mu <- mean(prey_details[which(prey_details$generation %in% (nrow(prey_res)-100):(nrow(prey_res)-100)),"lv"])
#var lv
sig <- var(prey_details[which(prey_details$generation %in% (nrow(prey_res)-100):(nrow(prey_res)-100)),"lv"])


plot(offspring ~ lv2,
     xlim = c(0.3,8),
     ylim = c(0,6),
     xlab = "Ballistic lengthscale (m)",
     ylab = "Prey fitness",
     data = AGG,
     pch = 20,
     col = adjustcolor("grey70", alpha = 0.6))
points(y ~ x,
       type = "l",
       xlab = "Ballistic lengthscale (m)",
       ylab = "Prey fitness")
abline(v = mu,
       col = '#046C9A',
       lty = 'dashed')
rect(xleft = mu - sig, xright = mu + sig,
     ybottom = par("usr")[3], ytop = par("usr")[4], 
     border = NA,
     col = adjustcolor("#046C9A", alpha = 0.3))



#Mean lv
mu <- mean(pred_details[which(pred_details$generation %in% (nrow(pred_res)-500):(nrow(pred_res)-500)),"lv"])
#var lv
sig <- var(pred_details[which(pred_details$generation %in% (nrow(pred_res)-500):(nrow(pred_res)-500)),"lv"])

pred_details$lv2 <- round(pred_details$lv,1)
AGG <- aggregate(offspring ~ lv2,
                 data = pred_details, FUN = "mean")
AGG <- AGG[-which(AGG$offspring > 9),]

FIT <- nls(offspring ~ a * (lv2 + c) * exp(-b * (lv2 + c)),
           start = list(a = 0.1,
                        b = 0.01,
                        c = -20),
           data = AGG)

ricker <- function(x) {
  coef(FIT)[1] * (x + coef(FIT)[3]) * exp(-coef(FIT)[2] * (x + coef(FIT)[3]))
}
x <- seq(0,130, 0.01)
y <- ricker(x)

plot(offspring ~ lv2,
     xlim = c(0,125),
     #ylim = c(0,9),
     xlab = "Ballistic lengthscale (m)",
     ylab = "Predator fitness",
     data = AGG,
     pch = 20,
     col = adjustcolor("grey70", alpha = 0.9))
points(y ~ x,
       type = "l",
       xlab = "Ballistic lengthscale (m)",
       ylab = "Prey fitness")

abline(v = mu,
       col = 'red',
       lty = 'dashed')
rect(xleft = mu - sig, xright = mu + sig,
     ybottom = par("usr")[3], ytop = par("usr")[4], 
     border = NA,
     col = adjustcolor("red", alpha = 0.3))

dev.off()
