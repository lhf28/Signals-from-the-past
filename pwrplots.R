setwd("/Users/lydiahfurness/Desktop/DissertationR")

############################################# POWER PLOTS #############################################

mytimes <- c(0.005,0.01,0.02,0.05,0.1,0.2)

load("Taj_Selsims_modified.RData")
load("Taj_neutral_sim_modified.RData")

nvector <- sim$taji_d
percentile <- quantile(nvector,0.01)

FPRs1 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector1 = Selsims[[1]][[a]][["taji_d"]]
  FPRs1[a] = length(vector1[vector1 <= percentile])/length(vector1)
}

FPRs2 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector2 = Selsims[[2]][[a]][["taji_d"]]
  FPRs2[a] = length(vector2[vector2 <= percentile])/length(vector2)
}

FPRs3 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector3 = Selsims[[3]][[a]][["taji_d"]]
  FPRs3[a] = length(vector3[vector3 <= percentile])/length(vector3)
}

FPRs4 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector4 = Selsims[[4]][[a]][["taji_d"]]
  FPRs4[a] = length(vector4[vector4 <= percentile])/length(vector4)
}


load("nuc_div_sim.RData")
load("nuc_div_Selsims.RData")

nvector <- sim$nuc_div
percentile <- quantile(nvector,0.01)

FPRs5 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector5 = Selsims[[1]][[a]][["nuc_div"]]
  FPRs5[a] = length(vector5[vector5 <= percentile])/length(vector5)
}

FPRs6 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector6 = Selsims[[2]][[a]][["nuc_div"]]
  FPRs6[a] = length(vector6[vector6 <= percentile])/length(vector6)
}

FPRs7 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector7 = Selsims[[3]][[a]][["nuc_div"]]
  FPRs7[a] = length(vector7[vector7 <= percentile])/length(vector7)
}

FPRs8 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector8 = Selsims[[4]][[a]][["nuc_div"]]
  FPRs8[a] = length(vector8[vector8 <= percentile])/length(vector8)
}


par(mfrow = c(2,2))
par(mar = c(4.0, 4.0, 2.0, 2.0))

### s = 0.01
coeff1taj <- c(FPRs1[[1]][1],FPRs1[[2]][1],FPRs1[[3]][1],FPRs1[[4]][1],FPRs1[[5]][1],FPRs1[[6]][1])
coeff1nucdiv <- c(FPRs5[[1]][1],FPRs5[[2]][1],FPRs5[[3]][1],FPRs5[[4]][1],FPRs5[[5]][1],FPRs5[[6]][1])
plot(coeff1taj, type = "o", lty = 2, lwd = 2, col = "pink", ylim = c(0.0:1.0), axes = FALSE, ann = FALSE)
axis(1, at=1:6, lab = c("5 kya", "10 kya", "20 kya", "50 kya", "100 kya", "200 kya"))
axis(2, tick = , las = 1, at = c(0.0,0.2,0.4,0.6,0.8,1.0))
box()
lines(coeff1nucdiv, type = "o", pch = 22, lty = 2, lwd = 2, col = "lightblue")
title(main = "s = 0.01")
title(xlab = "Time interval")
title(ylab = "Power")
legend("bottomright", c("Tajima's D", "Nucleotide diversity (π)"), col = c("pink", "lightblue"), cex = 0.9, pch = 21:22, lty = 2:2)
### s = 0.04
coeff2taj <- c(FPRs2[[1]][1],FPRs2[[2]][1],FPRs2[[3]][1],FPRs2[[4]][1],FPRs2[[5]][1],FPRs2[[6]][1])
coeff2nucdiv <- c(FPRs6[[1]][1],FPRs6[[2]][1],FPRs6[[3]][1],FPRs6[[4]][1],FPRs6[[5]][1],FPRs6[[6]][1])
plot(coeff2taj, type = "o", lty = 2, lwd = 2, col = "pink", ylim = c(0.0:1.0), axes = FALSE, ann = FALSE)
axis(1, at=1:6, lab = c("5 kya", "10 kya", "20 kya", "50 kya", "100 kya", "200 kya"))
axis(2, tick = , las = 1, at = c(0.0,0.2,0.4,0.6,0.8,1.0))
box()
lines(coeff2nucdiv, type = "o", pch = 22, lty = 2, lwd = 2, col = "lightblue")
title(main = "s = 0.04")
title(xlab = "Time interval")
title(ylab = "Power")
legend("bottomright", c("Tajima's D", "Nucleotide diversity (π)"), col = c("pink", "lightblue"), cex = 0.9, pch = 21:22, lty = 2:2)
### s = 0.06
coeff3taj <- c(FPRs3[[1]][1],FPRs3[[2]][1],FPRs3[[3]][1],FPRs3[[4]][1],FPRs3[[5]][1],FPRs3[[6]][1])
coeff3nucdiv <- c(FPRs7[[1]][1],FPRs7[[2]][1],FPRs7[[3]][1],FPRs7[[4]][1],FPRs7[[5]][1],FPRs7[[6]][1])
plot(coeff3taj, type = "o", lty = 2, lwd = 2, col = "pink", ylim = c(0.0:1.0), axes = FALSE, ann = FALSE)
axis(1, at=1:6, lab = c("5 kya", "10 kya", "20 kya", "50 kya", "100 kya", "200 kya"))
axis(2, tick = , las = 1, at = c(0.0,0.2,0.4,0.6,0.8,1.0))
box()
lines(coeff3nucdiv, type = "o", pch = 22, lty = 2, lwd = 2, col = "lightblue")
title(main = "s = 0.06")
title(xlab = "Time interval")
title(ylab = "Power")
legend("bottomright", c("Tajima's D", "Nucleotide diversity (π)"), col = c("pink", "lightblue"), cex = 0.9, pch = 21:22, lty = 2:2)
### s = 0.2
coeff4taj <- c(FPRs4[[1]][1],FPRs4[[2]][1],FPRs4[[3]][1],FPRs4[[4]][1],FPRs4[[5]][1],FPRs4[[6]][1])
coeff4nucdiv <- c(FPRs8[[1]][1],FPRs8[[2]][1],FPRs8[[3]][1],FPRs8[[4]][1],FPRs8[[5]][1],FPRs8[[6]][1])
plot(coeff4taj, type = "o", lty = 2, lwd = 2, col = "pink", ylim = c(0.0:1.0), axes = FALSE, ann = FALSE)
axis(1, at=1:6, lab = c("5 kya", "10 kya", "20 kya", "50 kya", "100 kya", "200 kya"))
axis(2, tick = , las = 1, at = c(0.0,0.2,0.4,0.6,0.8,1.0))
box()
lines(coeff4nucdiv, type = "o", pch = 22, lty = 2, lwd = 2, col = "lightblue")
title(main = "s = 0.2")
title(xlab = "Time interval")
title(ylab = "Power")
legend("bottomright", c("Tajima's D", "Nucleotide diversity (π)"), col = c("pink", "lightblue"), cex = 0.9, pch = 21:22, lty = 2:2)
