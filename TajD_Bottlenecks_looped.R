########################################################### SETUP ##################################################################
setwd("~/Desktop/DissertationR")

library(coala)
activate_msms(jar = "~/Desktop/DissertationR/msms/lib/msms.jar", priority = 500, download = FALSE)

list_simulators()
library(BSDA)

####################################################### NEUTRAL MODELS ############################################################ 
load("Taj_model.RData")
load("Taj_neutral_sim_modified.RData")

extrememodel <- model + feat_size_change(0.1, time = 0.005) + feat_size_change(1.0, time = 0)
check_model(extrememodel)
extremesim <- simulate(extrememodel)
save(extremesim, file = "extremesim.RData")

mildmodel <- model + feat_size_change(0.5, time = 0.005) + feat_size_change(1.0, time = 0)
check_model(mildmodel)
mildsim <- simulate(mildmodel)
save(mildsim, file = "mildsim.RData")

par(mfrow=c(1,1))
plot(density(sim$taji_d,na.rm=TRUE, adjust=1.5), lwd= 2, col = "black", main = "Distributions of Tajima's D under population bottleneck scenarios", xlab = "Tajima's D value", xlim = c(-2,2), ylim = c(0.0,2.3))
lines(density(extremesim$taji_d, na.rm=TRUE, adjust=1.5), lty= 3, lwd= 2, col = "green")
lines(density(mildsim$taji_d, na.rm=TRUE, adjust=1.5), lty= 3, lwd= 2, col = "blue")
colours=c("black","green","blue")
legend("topleft", 
       legend = c("Neutral","Extreme bottleneck","Mild bottleneck"), 
       fill = colours
)

### save plot

### STATS ###
neutral <- sim$taji_d
extreme <- extremesim$taji_d
mild <- mildsim$taji_d
qqnorm(neutral)
qqnorm(extreme)
qqnorm(mild)
sigma.y <- sd(neutral)
sigma.x1 <- sd(extreme)
sigma.x2 <- sd(mild)
z.test(extreme,neutral,"two.sided",sigma.x=sigma.x1,sigma.y = sigma.y)
z.test(mild,neutral,"two.sided",sigma.x=sigma.x2,sigma.y = sigma.y)


####################################################### SELECTION MODELS ############################################################
### extreme bottleneck situation - loss of 90% 5kya, regain Ne at time=0

mytimes <- c(0.005,0.01,0.02,0.05,0.1,0.2)

myselection <- vector(mode="list",length=2)
myselection[[1]] <- c(800,400) 
myselection[[2]] <- c(1200,600)

ExBotSelsims <- vector(mode="list", length=length(myselection))

for (a in 1:length(myselection)) {
  for (b in 1:length(mytimes)) {
    print(paste0("ExBotSelsims",a,"-",b," AA: ",myselection[[a]][1]," Aa: ",myselection[[a]][2]," time: ",mytimes[b]))
    SelMod1 <- extrememodel + feat_selection(
      strength_AA = myselection[[a]][1], 
      strength_Aa = myselection[[a]][2], 
      strength_aa = 0,
      strength_A = NULL,
      population = "all",
      time = mytimes[b], 
      start = TRUE,
      start_frequency = 5e-04, 
      Ne = 10000,
      position = 0.5,
      force_keep = TRUE,
      locus_group = "all")
    check_model(SelMod1)
    ExBotSelsims[[a]][[b]] <- simulate(SelMod1) 
  }
}

save(ExBotSelsims, file = "ExBotSelsims.RData")

labels <- c("(a)","(b)","(c)","(d)","(e)","(f)")
for (a in 1:length(ExBotSelsims)) {
  par(mfrow=c(3,2))
  par(mar = c(4.0, 4.0, 2.0, 2.0))
  myfile <- paste0("TajimasD_Extremebottleneck",myselection[[a]][1],"_",myselection[[a]][2],".pdf")
  for (b in 1:length(mytimes)) {
    mytitle <- paste0((labels[b]))
    plot(density(ExBotSelsims[[a]][[b]]$taji_d,na.rm=TRUE, adjust=1.5), main= mytitle, xlab = "Tajima's D value", xlim = c(-3,3), ylim = c(0.0,3.5))
    lines(density(extremesim$taji_d, na.rm=TRUE, adjust=1.5), lty= 3, lwd= 2, col = "green")
    lines(density(sim$taji_d, na.rm=TRUE, adjust=1.5), lty= 3, lwd= 2, col = "red")
  }
  dev.print(device=pdf,file=myfile)
}  


### a less extreme bottleneck situation - loss of 50% 5kya, regain Ne at time=0
MBotSelsims <- vector(mode="list", length=length(myselection))

for (a in 1:length(myselection)) {
  for (b in 1:length(mytimes)) {
    print(paste0("MBotSelsims",a,"-",b," AA: ",myselection[[a]][1]," Aa: ",myselection[[a]][2]," time: ",mytimes[b]))
    SelMod2 <- mildmodel + feat_selection(
      strength_AA = myselection[[a]][1], 
      strength_Aa = myselection[[a]][2], 
      strength_aa = 0,
      strength_A = NULL,
      population = "all",
      time = mytimes[b], 
      start = TRUE,
      start_frequency = 5e-04, 
      Ne = 10000,
      position = 0.5,
      force_keep = TRUE,
      locus_group = "all")
    check_model(SelMod2)
    MBotSelsims[[a]][[b]] <- simulate(SelMod2) 
  }
}

save(MBotSelsims, file = "MBotSelsims.RData")

labels <- c("(a)","(b)","(c)","(d)","(e)","(f)")
for (a in 1:length(MBotSelsims)) {
  par(mfrow=c(3,2))
  par(mar = c(4.0, 4.0, 2.0, 2.0))
  myfile <- paste0("TajimasD_Mildbottleneck",myselection[[a]][1],"_",myselection[[a]][2],".pdf")
  for (b in 1:length(mytimes)) {
    mytitle <- paste0((labels[b]))
    plot(density(MBotSelsims[[a]][[b]]$taji_d,na.rm=TRUE, adjust=1.5), main= mytitle, xlab = "Tajima's D value", xlim = c(-3,3), ylim = c(0.0,2.0))
    lines(density(mildsim$taji_d, na.rm=TRUE, adjust=1.5), lty= 3, lwd= 2, col = "blue")
    lines(density(sim$taji_d, na.rm=TRUE, adjust=1.5), lty= 3, lwd= 2, col = "red")
  }
  dev.print(device=pdf,file=myfile)
}  


####################################### STATS #########################################3

### EXTREME ###

load("ExBotSelsims.RData")

nvector <- extremesim$taji_d
percentile <- quantile(nvector,0.01)

FPRs1 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector = ExBotSelsims[[1]][[a]][["taji_d"]]
  FPRs1[a] = length(vector[vector <= percentile])/length(vector)
}

FPRs2 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector2 = ExBotSelsims[[2]][[a]][["taji_d"]]
  FPRs2[a] = length(vector2[vector2 <= percentile])/length(vector2)
}

save(FPRs1, file = "TajExBotpwr1.RData")
save(FPRs2, file = "TajExBotpwr2.RData")

### MILD ###

load("MBotSelsims.RData")

nvector <- mildsim$taji_d
percentile <- quantile(nvector,0.01)

FPRs3 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector3 = MBotSelsims[[1]][[a]][["taji_d"]]
  FPRs3[a] = length(vector3[vector3 <= percentile])/length(vector3)
}

FPRs4 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector4 = MBotSelsims[[2]][[a]][["taji_d"]]
  FPRs4[a] = length(vector4[vector4 <= percentile])/length(vector4)
}

save(FPRs3, file = "TajMBotpwr1.RData")
save(FPRs4, file = "TajMBotpwr2.RData")




