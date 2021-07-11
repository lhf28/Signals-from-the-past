########################################################### SETUP ##################################################################
setwd("~/Desktop/DissertationR")

library(coala)
activate_msms(jar = "~/Desktop/DissertationR/msms/lib/msms.jar", priority = 500, download = FALSE)

list_simulators()

########################################################### EXPANSION POWER DIFFERENCES ##################################################################

### NEUTRAL POWERS
load("Tajneutralpwr1.RData")
load("Tajneutralpwr2.RData")
load("nuc_divneutralpwr1.RData")
load("nuc_divneutralpwr2.RData")

### EXTREME
load("TajDExExpanpwr1.RData")
load("TajDExExpanpwr2.RData")
load("nuc_divExExpanpwr1.RData")
load("nuc_divExExpanpwr2.RData")

### MILD
load("TajDMExpanpwr1.RData")
load("TajDMExpanpwr2.RData")
load("nuc_divMExpanpwr1.RData")
load("nuc_divMExpanpwr2.RData")


### EXTREME EXPANSION POWERS

### s = 0.04
load("Tajneutralpwr1.RData")
coeff1taj <- c(FPRs2[[1]][1],FPRs2[[2]][1],FPRs2[[3]][1],FPRs2[[4]][1],FPRs2[[5]][1],FPRs2[[6]][1])
load("nuc_divneutralpwr1.RData")
coeff1nucdiv <- c(FPRs6[[1]][1],FPRs6[[2]][1],FPRs6[[3]][1],FPRs6[[4]][1],FPRs6[[5]][1],FPRs6[[6]][1])
load("TajDExExpanpwr1.RData")
coeff2taj <- c(FPRs1[[1]][1],FPRs1[[2]][1],FPRs1[[3]][1],FPRs1[[4]][1],FPRs1[[5]][1],FPRs1[[6]][1])
load("nuc_divExExpanpwr1.RData")
coeff2nucdiv <- c(FPRs1[[1]][1],FPRs1[[2]][1],FPRs1[[3]][1],FPRs1[[4]][1],FPRs1[[5]][1],FPRs1[[6]][1])
plot(coeff1taj, type = "o", lty = 2, lwd = 2, col = "pink", ylim = c(0.0:1.0), axes = FALSE, ann = FALSE)
axis(1, at=1:6, lab = c("5 kya", "10 kya", "20 kya", "50 kya", "100 kya", "200 kya"))
axis(2, tick = , las = 1, at = c(0.0,0.2,0.4,0.6,0.8,1.0))
box()
lines(coeff2taj, type = "o", pch = 21, lty = "solid", lwd = 2, col = "pink")
lines(coeff1nucdiv, type = "o", pch = 21, lty = 2, lwd = 2, col = "lightblue")
lines(coeff2nucdiv, type = "o", pch = 21, lty = "solid", lwd = 2, col = "lightblue")
title(main = "s = 0.04")
title(xlab = "Time interval")
title(ylab = "Power")
legend("bottomright", cex = 0.9, c("Neutral Tajima's D", "Neutral nucleotide diversity (π)","Extreme expansion Tajima's D", "Extreme expansion nucleotide diversity (π)"), col = c("pink", "lightblue", "pink", "lightblue"), pch = 21, lty = c(2,2,1,1))
### s = 0.06
load("Tajneutralpwr2.RData")
coeff3taj <- c(FPRs3[[1]][1],FPRs3[[2]][1],FPRs3[[3]][1],FPRs3[[4]][1],FPRs3[[5]][1],FPRs3[[6]][1])
load("nuc_divneutralpwr2.RData")
coeff3nucdiv <- c(FPRs7[[1]][1],FPRs7[[2]][1],FPRs7[[3]][1],FPRs7[[4]][1],FPRs7[[5]][1],FPRs7[[6]][1])
load("TajDExExpanpwr2.RData")
coeff4taj <- c(FPRs2[[1]][1],FPRs2[[2]][1],FPRs2[[3]][1],FPRs2[[4]][1],FPRs2[[5]][1],FPRs2[[6]][1])
load("nuc_divExExpanpwr2.RData")
coeff4nucdiv <- c(FPRs2[[1]][1],FPRs2[[2]][1],FPRs2[[3]][1],FPRs2[[4]][1],FPRs2[[5]][1],FPRs2[[6]][1])
plot(coeff3taj, type = "o", lty = 2, lwd = 2, col = "pink", ylim = c(0.0:1.0), axes = FALSE, ann = FALSE)
axis(1, at=1:6, lab = c("5 kya", "10 kya", "20 kya", "50 kya", "100 kya", "200 kya"))
axis(2, tick = , las = 1, at = c(0.0,0.2,0.4,0.6,0.8,1.0))
box()
lines(coeff4taj, type = "o", pch = 21, lty = "solid", lwd = 2, col = "pink")
lines(coeff3nucdiv, type = "o", pch = 21, lty = 2, lwd = 2, col = "lightblue")
lines(coeff4nucdiv, type = "o", pch = 21, lty = "solid", lwd = 2, col = "lightblue")
title(main = "s = 0.06")
title(xlab = "Time interval")
title(ylab = "Power")
legend("bottomright", cex = 0.9, c("Neutral Tajima's D", "Neutral nucleotide diversity (π)","Extreme expansion Tajima's D", "Extreme expansion nucleotide diversity (π)"), col = c("pink", "lightblue", "pink", "lightblue"), pch = 21, lty = c(2,2,1,1))


### MILD BOTTLENECK POWERS

### s = 0.04
load("Tajneutralpwr1.RData")
coeff1taj <- c(FPRs2[[1]][1],FPRs2[[2]][1],FPRs2[[3]][1],FPRs2[[4]][1],FPRs2[[5]][1],FPRs2[[6]][1])
load("nuc_divneutralpwr1.RData")
coeff1nucdiv <- c(FPRs6[[1]][1],FPRs6[[2]][1],FPRs6[[3]][1],FPRs6[[4]][1],FPRs6[[5]][1],FPRs6[[6]][1])
load("TajDMExpanpwr1.RData")
coeff2taj <- c(FPRs3[[1]][1],FPRs3[[2]][1],FPRs3[[3]][1],FPRs3[[4]][1],FPRs3[[5]][1],FPRs3[[6]][1])
load("nuc_divMExpanpwr1.RData")
coeff2nucdiv <- c(FPRs3[[1]][1],FPRs3[[2]][1],FPRs3[[3]][1],FPRs3[[4]][1],FPRs3[[5]][1],FPRs3[[6]][1])
plot(coeff1taj, type = "o", lty = 2, lwd = 2, col = "pink", ylim = c(0.0:1.0), axes = FALSE, ann = FALSE)
axis(1, at=1:6, lab = c("5 kya", "10 kya", "20 kya", "50 kya", "100 kya", "200 kya"))
axis(2, tick = , las = 1, at = c(0.0,0.2,0.4,0.6,0.8,1.0))
box()
lines(coeff2taj, type = "o", pch = 21, lty = "solid", lwd = 2, col = "pink")
lines(coeff1nucdiv, type = "o", pch = 21, lty = 2, lwd = 2, col = "lightblue")
lines(coeff2nucdiv, type = "o", pch = 21, lty = "solid", lwd = 2, col = "lightblue")
title(main = "s = 0.04")
title(xlab = "Time interval")
title(ylab = "Power")
legend("bottomright", cex = 0.9, c("Neutral Tajima's D", "Neutral nucleotide diversity (π)","Extreme expansion Tajima's D", "Extreme expansion nucleotide diversity (π)"), col = c("pink", "lightblue", "pink", "lightblue"), pch = 21, lty = c(2,2,1,1))
### s = 0.06
load("Tajneutralpwr2.RData")
coeff3taj <- c(FPRs3[[1]][1],FPRs3[[2]][1],FPRs3[[3]][1],FPRs3[[4]][1],FPRs3[[5]][1],FPRs3[[6]][1])
load("nuc_divneutralpwr2.RData")
coeff3nucdiv <- c(FPRs7[[1]][1],FPRs7[[2]][1],FPRs7[[3]][1],FPRs7[[4]][1],FPRs7[[5]][1],FPRs7[[6]][1])
load("TajDMExpanpwr2.RData")
coeff4taj <- c(FPRs4[[1]][1],FPRs4[[2]][1],FPRs4[[3]][1],FPRs4[[4]][1],FPRs4[[5]][1],FPRs4[[6]][1])
load("nuc_divMExpanpwr2.RData")
coeff4nucdiv <- c(FPRs4[[1]][1],FPRs4[[2]][1],FPRs4[[3]][1],FPRs4[[4]][1],FPRs4[[5]][1],FPRs4[[6]][1])
plot(coeff3taj, type = "o", lty = 2, lwd = 2, col = "pink", ylim = c(0.0:1.0), axes = FALSE, ann = FALSE)
axis(1, at=1:6, lab = c("5 kya", "10 kya", "20 kya", "50 kya", "100 kya", "200 kya"))
axis(2, tick = , las = 1, at = c(0.0,0.2,0.4,0.6,0.8,1.0))
box()
lines(coeff4taj, type = "o", pch = 21, lty = "solid", lwd = 2, col = "pink")
lines(coeff3nucdiv, type = "o", pch = 21, lty = 2, lwd = 2, col = "lightblue")
lines(coeff4nucdiv, type = "o", pch = 21, lty = "solid", lwd = 2, col = "lightblue")
title(main = "s = 0.06")
title(xlab = "Time interval")
title(ylab = "Power")
legend("bottomright", cex = 0.9, c("Neutral Tajima's D", "Neutral nucleotide diversity (π)","Extreme expansion Tajima's D", "Extreme expansion nucleotide diversity (π)"), col = c("pink", "lightblue", "pink", "lightblue"), pch = 21, lty = c(2,2,1,1))
