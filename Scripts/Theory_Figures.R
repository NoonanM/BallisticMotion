x <- seq(1, 6, by = 0.001)

y_prey <- 1/exp(x)
y_pred <- (0.1*(x) + 0.3)
ratio <- y_pred/y_prey



png("Results/Encounter_Theory.png",
    width=6.86,
    height=4,
    units = "in",
    res = 600) 

par(mfrow=c(3,2), mgp=c(0.2,0.2,0), tcl=-0.4, mar=c(2,10,1,1), font = 1, font.axis = 2, font.lab = 2)
plot(x = x,
     y = y_prey,
     type = "l",
     yaxt='n',
     xaxt='n',
     col = "#3471bc",
     xlab = "Generation",
     ylab = "Ballistic length scale",
     ylim = c(0,1),
     cex.lab = 0.6)

lines(x = x,
     y = y_pred,
     type = "l",
     yaxt='n',
     xaxt='n',
     col = "#e6c141",
     xlab = "Generation",
     ylab = "Ballistic length scale",
     cex.lab = 0.6)


plot(x = x,
     y = ratio,
     type = "l",
     yaxt='n',
     xaxt='n',
     col = "#3c7a47",
     xlab = "Generation",
     ylab = "Predator:Prey",
     cex.lab = 0.6)


y_prey <- (0.1*(x))
y_pred <- 0.1*(x) + 0.3
ratio <- y_pred/y_prey



plot(x = x,
     y = y_prey,
     type = "l",
     yaxt='n',
     xaxt='n',
     col = "#3471bc",
     xlab = "Generation",
     ylab = "Ballistic length scale",
     ylim = c(0,1),
     cex.lab = 0.6)

lines(x = x,
      y = y_pred,
      type = "l",
      yaxt='n',
      xaxt='n',
      col = "#e6c141",
      xlab = "Generation",
      ylab = "Ballistic length scale",
      cex.lab = 0.6)

plot(x = x,
     y = ratio,
     type = "l",
     yaxt='n',
     xaxt='n',
     col = "#3c7a47",
     xlab = "Generation",
     ylab = "Predator:Prey",
     cex.lab = 0.6)




x <- seq(1, 200, by = 0.001)
y_prey <- (0.2 * x)/(3 + x) #rep(0.3, length(x))
y_pred <- (0.3 * x)/(3 + x) #rep(0.2, length(x))
ratio <- y_pred/y_prey



plot(x = x,
     y = y_prey,
     type = "l",
     yaxt='n',
     xaxt='n',
     col = "#3471bc",
     xlab = "Generation",
     ylab = "Ballistic length scale",
     ylim = c(0,0.5),
     cex.lab = 0.6)

lines(x = x,
      y = y_pred,
      type = "l",
      yaxt='n',
      xaxt='n',
      col = "#e6c141",
      xlab = "Generation",
      ylab = "Ballistic length scale",
      cex.lab = 0.6)

plot(x = x,
     y = ratio,
     type = "l",
     yaxt='n',
     xaxt='n',
     col = "#3c7a47",
     xlab = "Generation",
     ylab = "Predator:Prey",
     cex.lab = 0.6)

dev.off()