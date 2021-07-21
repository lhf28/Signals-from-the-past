########################################################### SETUP ###################################################################
setwd("~/Desktop/DissertationR")

library(coala)
activate_msms(jar = "~/Desktop/DissertationR/msms/lib/msms.jar", priority = 500, download = FALSE)

list_simulators()

############################################################ POWER ##################################################################

### TAJ D
load("Selsims1OoA.RData")
load("Selsims2OoA.RData")
load("Selsims3OoA.RData")
load("Selsims4OoA.RData")

### NUC_DIV
load("Selsims5OoA.RData")
load("Selsims6OoA.RData")
load("Selsims7OoA.RData")
load("Selsims8OoA.RData")


mytimes <- c(0.005,0.01,0.02,0.05,0.1,0.2)

load("TajOoAsim.RData")
load("nuc_divOoAsim.RData")
TajDnvector <- TajOoAsim$taji_d
Tajpercentile <- quantile(TajDnvector,0.01)
Nuc_divnvector <- nuc_divOoAsim$nuc_div
Nuc_divpercentile <- quantile(Nuc_divnvector,0.01)

FPRsTajAfrica1 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector1 = Selsims1[[1]][[a]][["taji_d"]]
  FPRsTajAfrica1[a] = length(vector1[vector1 <= Tajpercentile])/length(vector1)
}

FPRsTajAfrica2 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector2 = Selsims1[[2]][[a]][["taji_d"]]
  FPRsTajAfrica2[a] = length(vector2[vector2 <= Tajpercentile])/length(vector2)
}

FPRsTajEurope1 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector3 = Selsims2[[1]][[a]][["taji_d"]]
  FPRsTajEurope1[a] = length(vector3[vector3 <= Tajpercentile])/length(vector3)
}

FPRsTajEurope2 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector4 = Selsims2[[2]][[a]][["taji_d"]]
  FPRsTajEurope2[a] = length(vector4[vector4 <= Tajpercentile])/length(vector4)
}

FPRsTajAsia1 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector5 = Selsims3[[1]][[a]][["taji_d"]]
  FPRsTajAsia1[a] = length(vector5[vector5 <= Tajpercentile])/length(vector5)
}

FPRsTajAsia2 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector6 = Selsims3[[2]][[a]][["taji_d"]]
  FPRsTajAsia2[a] = length(vector6[vector6 <= Tajpercentile])/length(vector6)
}

FPRsTajMexAm1 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector7 = Selsims4[[1]][[a]][["taji_d"]]
  FPRsTajMexAm1[a] = length(vector7[vector7 <= Tajpercentile])/length(vector7)
}

FPRsTajMexAm2 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector8 = Selsims4[[2]][[a]][["taji_d"]]
  FPRsTajMexAm2[a] = length(vector8[vector8 <= Tajpercentile])/length(vector8)
}

###

FPRsnuc_divAfrica1 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector9 = Selsims5[[1]][[a]][["nuc_div"]]
  FPRsnuc_divAfrica1[a] = length(vector9[vector9 <= Nuc_divpercentile])/length(vector9)
}

FPRsnuc_divAfrica2 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector10 = Selsims5[[2]][[a]][["nuc_div"]]
  FPRsnuc_divAfrica2[a] = length(vector10[vector10 <= Nuc_divpercentile])/length(vector10)
}

FPRsnuc_divEurope1 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector11 = Selsims6[[1]][[a]][["nuc_div"]]
  FPRsnuc_divEurope1[a] = length(vector11[vector11 <= Nuc_divpercentile])/length(vector11)
}

FPRsnuc_divEurope2 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector12 = Selsims6[[2]][[a]][["nuc_div"]]
  FPRsnuc_divEurope2[a] = length(vector12[vector12 <= Nuc_divpercentile])/length(vector12)
}

FPRsnuc_divAsia1 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector13 = Selsims7[[1]][[a]][["nuc_div"]]
  FPRsnuc_divAsia1[a] = length(vector13[vector13 <= Nuc_divpercentile])/length(vector13)
}

FPRsnuc_divAsia2 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector14 = Selsims7[[2]][[a]][["nuc_div"]]
  FPRsnuc_divAsia2[a] = length(vector14[vector14 <= Nuc_divpercentile])/length(vector14)
}

FPRsnuc_divMexAm1 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector15 = Selsims8[[1]][[a]][["nuc_div"]]
  FPRsnuc_divMexAm1[a] = length(vector15[vector15 <= Nuc_divpercentile])/length(vector15)
}

FPRsnuc_divMexAm2 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector16 = Selsims8[[2]][[a]][["nuc_div"]]
  FPRsnuc_divMexAm2[a] = length(vector16[vector16 <= Nuc_divpercentile])/length(vector16)
}


### PLOTTING POWER

par(mfrow = c(2,2))
par(mar = c(4.0, 4.0, 2.0, 2.0))

# 0.04 Taj D
load("Tajneutralpwr1.RData")
neutral <- c(FPRs2[[1]][1],FPRs2[[2]][1],FPRs2[[3]][1],FPRs2[[4]][1],FPRs2[[5]][1],FPRs2[[6]][1])
Africa <- c(FPRsTajAfrica1[[1]][1],FPRsTajAfrica1[[2]][1],FPRsTajAfrica1[[3]][1],FPRsTajAfrica1[[4]][1],FPRsTajAfrica1[[5]][1],FPRsTajAfrica1[[6]][1])
Europe <- c(FPRsTajEurope1[[1]][1],FPRsTajEurope1[[2]][1],FPRsTajEurope1[[3]][1],FPRsTajEurope1[[4]][1],FPRsTajEurope1[[5]][1],FPRsTajEurope1[[6]][1])
Asia <- c(FPRsTajAsia1[[1]][1],FPRsTajAsia1[[2]][1],FPRsTajAsia1[[3]][1],FPRsTajAsia1[[4]][1],FPRsTajAsia1[[5]][1],FPRsTajAsia1[[6]][1])
MexAm <- c(FPRsTajMexAm1[[1]][1],FPRsTajMexAm1[[2]][1],FPRsTajMexAm1[[3]][1],FPRsTajMexAm1[[4]][1],FPRsTajMexAm1[[5]][1],FPRsTajMexAm1[[6]][1])
plot(neutral, type = "o", lty = 2, lwd = 2, col = "pink", ylim = c(0.0:1.0), axes = FALSE, ann = FALSE)
axis(1, at=1:6, lab = c("5 kya", "10 kya", "20 kya", "50 kya", "100 kya", "200 kya"))
axis(2, tick = , las = 1, at = c(0.0,0.2,0.4,0.6,0.8,1.0))
box()
lines(Africa, type = "o", pch = 21, lty = 2, lwd = 2, col = "red")
lines(Europe, type = "o", pch = 21, lty = 2, lwd = 2, col = "blue")
lines(Asia, type = "o", pch = 21, lty = 2, lwd = 2, col = "green")
lines(MexAm, type = "o", pch = 21, lty = 2, lwd = 2, col = "orange")
load("OoATajpwr1.RData")
Pooled <- c(FPRs1[[1]][1],FPRs1[[2]][1],FPRs1[[3]][1],FPRs1[[4]][1],FPRs1[[5]][1],FPRs1[[6]][1])
lines(Pooled, type = "o", pch = 21, lty = "solid", lwd = 2, col = "pink")
title(main = "s = 0.04")
title(xlab = "Time interval")
title(ylab = "Power")
legend("right", c("Neutral", "African", "European", "Asian", "Mexican American", "Pooled"), col = c("pink", "red", "blue", "green", "orange", "pink"), cex = 0.9, pch = 21, lty = c(2,2,2,2,2,1))

# 0.06 Taj D
load("Tajneutralpwr2.RData")
neutral <- c(FPRs3[[1]][1],FPRs3[[2]][1],FPRs3[[3]][1],FPRs3[[4]][1],FPRs3[[5]][1],FPRs3[[6]][1])
Africa <- c(FPRsTajAfrica2[[1]][1],FPRsTajAfrica2[[2]][1],FPRsTajAfrica2[[3]][1],FPRsTajAfrica2[[4]][1],FPRsTajAfrica2[[5]][1],FPRsTajAfrica2[[6]][1])
Europe <- c(FPRsTajEurope2[[1]][1],FPRsTajEurope2[[2]][1],FPRsTajEurope2[[3]][1],FPRsTajEurope2[[4]][1],FPRsTajEurope2[[5]][1],FPRsTajEurope2[[6]][1])
Asia <- c(FPRsTajAsia2[[1]][1],FPRsTajAsia2[[2]][1],FPRsTajAsia2[[3]][1],FPRsTajAsia2[[4]][1],FPRsTajAsia2[[5]][1],FPRsTajAsia2[[6]][1])
MexAm <- c(FPRsTajMexAm2[[1]][1],FPRsTajMexAm2[[2]][1],FPRsTajMexAm2[[3]][1],FPRsTajMexAm2[[4]][1],FPRsTajMexAm2[[5]][1],FPRsTajMexAm2[[6]][1])
plot(neutral, type = "o", lty = 2, lwd = 2, col = "pink", ylim = c(0.0:1.0), axes = FALSE, ann = FALSE)
axis(1, at=1:6, lab = c("5 kya", "10 kya", "20 kya", "50 kya", "100 kya", "200 kya"))
axis(2, tick = , las = 1, at = c(0.0,0.2,0.4,0.6,0.8,1.0))
box()
lines(Africa, type = "o", pch = 21, lty = 2, lwd = 2, col = "red")
lines(Europe, type = "o", pch = 21, lty = 2, lwd = 2, col = "blue")
lines(Asia, type = "o", pch = 21, lty = 2, lwd = 2, col = "green")
lines(MexAm, type = "o", pch = 21, lty = 2, lwd = 2, col = "orange")
load("OoATajpwr2.RData")
Pooled <- c(FPRs2[[1]][1],FPRs2[[2]][1],FPRs2[[3]][1],FPRs2[[4]][1],FPRs2[[5]][1],FPRs2[[6]][1])
lines(Pooled, type = "o", pch = 21, lty = "solid", lwd = 2, col = "pink")
title(main = "s = 0.06")
title(xlab = "Time interval")
title(ylab = "Power")
legend("right", c("Neutral", "African", "European", "Asian", "Mexican American", "Pooled"), col = c("pink", "red", "blue", "green", "orange", "pink"), cex = 0.9, pch = 21, lty = c(2,2,2,2,2,1))

# 0.04 Taj D
load("nuc_divneutralpwr1.RData")
neutral <- c(FPRs6[[1]][1],FPRs6[[2]][1],FPRs6[[3]][1],FPRs6[[4]][1],FPRs6[[5]][1],FPRs6[[6]][1])
Africa <- c(FPRsnuc_divAfrica1[[1]][1],FPRsnuc_divAfrica1[[2]][1],FPRsnuc_divAfrica1[[3]][1],FPRsnuc_divAfrica1[[4]][1],FPRsnuc_divAfrica1[[5]][1],FPRsnuc_divAfrica1[[6]][1])
Europe <- c(FPRsnuc_divEurope1[[1]][1],FPRsnuc_divEurope1[[2]][1],FPRsnuc_divEurope1[[3]][1],FPRsnuc_divEurope1[[4]][1],FPRsnuc_divEurope1[[5]][1],FPRsnuc_divEurope1[[6]][1])
Asia <- c(FPRsnuc_divAsia1[[1]][1],FPRsnuc_divAsia1[[2]][1],FPRsnuc_divAsia1[[3]][1],FPRsnuc_divAsia1[[4]][1],FPRsnuc_divAsia1[[5]][1],FPRsnuc_divAsia1[[6]][1])
MexAm <- c(FPRsnuc_divMexAm1[[1]][1],FPRsnuc_divMexAm1[[2]][1],FPRsnuc_divMexAm1[[3]][1],FPRsnuc_divMexAm1[[4]][1],FPRsnuc_divMexAm1[[5]][1],FPRsnuc_divMexAm1[[6]][1])
plot(neutral, type = "o", lty = 2, lwd = 2, col = "lightblue", ylim = c(0.0:1.0), axes = FALSE, ann = FALSE)
axis(1, at=1:6, lab = c("5 kya", "10 kya", "20 kya", "50 kya", "100 kya", "200 kya"))
axis(2, tick = , las = 1, at = c(0.0,0.2,0.4,0.6,0.8,1.0))
box()
lines(Africa, type = "o", pch = 21, lty = 2, lwd = 2, col = "red")
lines(Europe, type = "o", pch = 21, lty = 2, lwd = 2, col = "blue")
lines(Asia, type = "o", pch = 21, lty = 2, lwd = 2, col = "green")
lines(MexAm, type = "o", pch = 21, lty = 2, lwd = 2, col = "orange")
load("OoAnuc_divpwr1.RData")
Pooled <- c(FPRs3[[1]][1],FPRs3[[2]][1],FPRs3[[3]][1],FPRs3[[4]][1],FPRs3[[5]][1],FPRs3[[6]][1])
lines(Pooled, type = "o", pch = 21, lty = "solid", lwd = 2, col = "lightblue")
title(main = "s = 0.04")
title(xlab = "Time interval")
title(ylab = "Power")
legend("bottomright", c("Neutral", "African", "European", "Asian", "Mexican American", "Pooled"), col = c("lightblue", "red", "blue", "green", "orange", "lightblue"), cex = 0.9, pch = 21, lty = c(2,2,2,2,2,1))

# 0.06 nuc_div
load("nuc_divneutralpwr2.RData")
neutral <- c(FPRs7[[1]][1],FPRs7[[2]][1],FPRs7[[3]][1],FPRs7[[4]][1],FPRs7[[5]][1],FPRs7[[6]][1])
Africa <- c(FPRsnuc_divAfrica2[[1]][1],FPRsnuc_divAfrica2[[2]][1],FPRsnuc_divAfrica2[[3]][1],FPRsnuc_divAfrica2[[4]][1],FPRsnuc_divAfrica2[[5]][1],FPRsnuc_divAfrica2[[6]][1])
Europe <- c(FPRsnuc_divEurope2[[1]][1],FPRsnuc_divEurope2[[2]][1],FPRsnuc_divEurope2[[3]][1],FPRsnuc_divEurope2[[4]][1],FPRsnuc_divEurope2[[5]][1],FPRsnuc_divEurope2[[6]][1])
Asia <- c(FPRsnuc_divAsia2[[1]][1],FPRsnuc_divAsia2[[2]][1],FPRsnuc_divAsia2[[3]][1],FPRsnuc_divAsia2[[4]][1],FPRsnuc_divAsia2[[5]][1],FPRsnuc_divAsia2[[6]][1])
MexAm <- c(FPRsnuc_divMexAm2[[1]][1],FPRsnuc_divMexAm2[[2]][1],FPRsnuc_divMexAm2[[3]][1],FPRsnuc_divMexAm2[[4]][1],FPRsnuc_divMexAm2[[5]][1],FPRsnuc_divMexAm2[[6]][1])
plot(neutral, type = "o", lty = 2, lwd = 2, col = "lightblue", ylim = c(0.0:1.0), axes = FALSE, ann = FALSE)
axis(1, at=1:6, lab = c("5 kya", "10 kya", "20 kya", "50 kya", "100 kya", "200 kya"))
axis(2, tick = , las = 1, at = c(0.0,0.2,0.4,0.6,0.8,1.0))
box()
lines(Africa, type = "o", pch = 21, lty = 2, lwd = 2, col = "red")
lines(Europe, type = "o", pch = 21, lty = 2, lwd = 2, col = "blue")
lines(Asia, type = "o", pch = 21, lty = 2, lwd = 2, col = "green")
lines(MexAm, type = "o", pch = 21, lty = 2, lwd = 2, col = "orange")
load("OoAnuc_divpwr2.RData")
Pooled <- c(FPRs4[[1]][1],FPRs4[[2]][1],FPRs4[[3]][1],FPRs4[[4]][1],FPRs4[[5]][1],FPRs4[[6]][1])
lines(Pooled, type = "o", pch = 21, lty = "solid", lwd = 2, col = "lightblue")
title(main = "s = 0.06")
title(xlab = "Time interval")
title(ylab = "Power")
legend("bottomright", c("Neutral", "African", "European", "Asian", "Mexican American", "Pooled"), col = c("lightblue", "red", "blue", "green", "orange", "lightblue"), cex = 0.9, pch = 21, lty = c(2,2,2,2,2,1))

