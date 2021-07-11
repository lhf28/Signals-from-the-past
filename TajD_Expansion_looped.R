########################################################### SETUP ##################################################################
setwd("~/Desktop/DissertationR")

library(coala)
activate_msms(jar = "~/Desktop/DissertationR/msms/lib/msms.jar", priority = 500, download = FALSE)

list_simulators()
library(BSDA)

####################################################### NEUTRAL MODELS ############################################################ 
load("Taj_model.RData")
load("Taj_neutral_sim_modified.RData")

extrememodel2 <- model + feat_size_change(1.0, time = 0.01) + feat_size_change(10, time = 0)
check_model(extrememodel2)
extremesim2 <- simulate(extrememodel2)
save(extremesim2, file = "updated_extremesim2.RData")

mildmodel2 <- model + feat_size_change(1.0, time = 0.01)+ feat_size_change(5, time = 0)
check_model(mildmodel2)
mildsim2 <- simulate(mildmodel2)
save(mildsim2, file = "updated_mildsim2.RData")

par(mfrow=c(1,1))
plot(density(sim$taji_d,na.rm=TRUE, adjust=1.5), lwd= 2, col = "black", main = "Distributions of Tajima's D under population expansion scenarios", xlab = "Tajima's D value", xlim = c(-1.5,1.5), ylim = c(0.0,2.3))
lines(density(extremesim2$taji_d, na.rm=TRUE, adjust=1.5), lty= 3, lwd= 2, col = "green")
lines(density(mildsim2$taji_d, na.rm=TRUE, adjust=1.5), lty= 3, lwd= 2, col = "blue")
colours=c("black","green","blue")
legend("topleft", 
       legend = c("Neutral","Extreme expansion","Mild expansion"), 
       fill = colours
)

### STATS ###
neutral <- sim$taji_d
extreme2 <- extremesim2$taji_d
mild2 <- mildsim2$taji_d
qqnorm(neutral)
qqnorm(extreme2)
qqnorm(mild2)
sigma.y <- sd(neutral)
sigma.x1 <- sd(extreme2)
sigma.x2 <- sd(mild2)
z.test(extreme2,neutral,"two.sided",sigma.x=sigma.x1,sigma.y = sigma.y)
z.test(mild2,neutral,"two.sided",sigma.x=sigma.x2,sigma.y = sigma.y)


####################################################### SELECTION MODELS ############################################################

### more extreme population expansion - Ne increases from 10,000 to 100,000 individuals from 10kya ###

mytimes <- c(0.005,0.01,0.02,0.05,0.1,0.2)

myselection <- vector(mode="list",length=2)
myselection[[1]] <- c(800,400) 
myselection[[2]] <- c(1200,600)
Selsims <- vector(mode="list", length=length(myselection))

for (a in 1:length(myselection)) {
  for (b in 1:length(mytimes)) {
    print(paste0("Selsims",a,"-",b," AA: ",myselection[[a]][1]," Aa: ",myselection[[a]][2]," time: ",mytimes[b]))
    SelMod1 <- extrememodel2 + feat_selection(
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
    Selsims[[a]][[b]] <- simulate(SelMod1) 
  }
}

save(Selsims, file = "ExtremeExpansionSims.RData")

labels <- c("(a)","(b)","(c)","(d)","(e)","(f)")
for (a in 1:length(Selsims)) {
  par(mfrow=c(3,2))
  par(mar = c(4.0, 4.0, 2.0, 2.0))
  myfile <- paste0("TajimasD_ExtremeExpansion_updated",myselection[[a]][1],"_",myselection[[a]][2],".pdf")
  for (b in 1:length(mytimes)) {
    mytitle <- paste0((labels[b]))
    plot(density(Selsims[[a]][[b]]$taji_d,na.rm=TRUE, adjust=1.5), main= mytitle, xlab = "Tajima's D value", xlim = c(-3,3), ylim = c(0.0,2.3))
    lines(density(extremesim2$taji_d, na.rm = TRUE, adjust = 1.5),lty= 3, lwd= 2, col = "green")
    lines(density(sim$taji_d, na.rm = TRUE, adjust = 1.5),lty= 3, lwd= 2, col = "red")
  }
  dev.print(device=pdf,file=myfile)
}  


### a less extreme expansion event - Ne increases from 10,000 to 50,000 individuals from 10kya ###

for (a in 1:length(myselection)) {
  for (b in 1:length(mytimes)) {
    print(paste0("Selsims",a,"-",b," AA: ",myselection[[a]][1]," Aa: ",myselection[[a]][2]," time: ",mytimes[b]))
    SelMod2 <- mildmodel2 + feat_selection(
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
    Selsims[[a]][[b]] <- simulate(SelMod2) 
  }
}

save(Selsims, file = "MildExpansionSims.RData")

labels <- c("(a)","(b)","(c)","(d)","(e)","(f)")
for (a in 1:length(Selsims)) {
  par(mfrow=c(3,2))
  par(mar = c(4.0, 4.0, 2.0, 2.0))
  myfile <- paste0("TajimasD_MildExpansion",myselection[[a]][1],"_",myselection[[a]][2],".pdf")
  for (b in 1:length(mytimes)) {
    mytitle <- paste0((labels[b]))
    plot(density(Selsims[[a]][[b]]$taji_d,na.rm=TRUE, adjust=1.5), main= mytitle, xlab = "Tajima's D value", xlim = c(-3,3), ylim = c(0.0,2.0))
    lines(density(mildsim2$taji_d, na.rm = TRUE, adjust = 1.5),lty= 3, lwd= 2, col = "blue")
    lines(density(sim$taji_d, na.rm = TRUE, adjust = 1.5),lty= 3, lwd= 2, col = "red")
  }
  dev.print(device=pdf,file=myfile)
}  


########################################### STATS ###########################################
nvector <- extremesim2$taji_d
percentile <- quantile(nvector,0.01)

### EXTREME ###

load("ExtremeExpansionSims.RData")

FPRs1 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector = Selsims[[1]][[a]][["taji_d"]]
  FPRs1[a] = length(vector[vector <= percentile])/length(vector)
}

FPRs2 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector2 = Selsims[[2]][[a]][["taji_d"]]
  FPRs2[a] = length(vector2[vector2 <= percentile])/length(vector2)
}

save(FPRs1, file = "TajDExExpanpwr1.RData")
save(FPRs2, file = "TajDExExpanpwr2.RData")

### MILD ###

nvector <- mildsim2$taji_d
percentile <- quantile(nvector,0.01)

load("MildExpansionSims.RData")

FPRs3 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector3 = Selsims[[1]][[a]][["taji_d"]]
  FPRs3[a] = length(vector3[vector3 <= percentile])/length(vector3)
}

FPRs4 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector4 = Selsims[[2]][[a]][["taji_d"]]
  FPRs4[a] = length(vector4[vector4 <= percentile])/length(vector4)
}

save(FPRs3, file = "TajDMExpanpwr1.RData")
save(FPRs4, file = "TajDMExpanpwr2.RData")
