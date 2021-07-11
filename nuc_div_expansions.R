########################################################### SETUP ##################################################################
setwd("~/Desktop/DissertationR")

library(coala)
activate_msms(jar = "~/Desktop/DissertationR/msms/lib/msms.jar", priority = 500, download = FALSE)

list_simulators()
library(BSDA)

####################################################### NEUTRAL MODELS ############################################################ 
load("nuc_div_model.RData")
load("nuc_div_sim.RData")

extrememodel2 <- model + feat_size_change(1.0, time = 0.01) + feat_size_change(10, time = 0)
check_model(extrememodel2)
extremesim2 <- simulate(extrememodel2)
save(extremesim2, file = "nuc_div_extremesim2.RData")

mildmodel2 <- model + feat_size_change(1.0, time = 0.01) + feat_size_change(5, time = 0)
check_model(mildmodel2)
mildsim2 <- simulate(mildmodel2)
save(mildsim2, file = "nuc_div_mildsim2.RData")

par(mfrow=c(1,1))
plot(density(sim$nuc_div,na.rm=TRUE, adjust=1.5), lwd= 2, col = "black", main = "Distributions of nucleotide diversity (π) under population expansion scenarios", xlab = "Nucleotide diversity (π) ", xlim = c(-5,200), ylim = c(0.0,0.05))
lines(density(extremesim2$nuc_div, na.rm=TRUE, adjust=1.5), lty= 3, lwd= 2, col = "green")
lines(density(mildsim2$nuc_div, na.rm=TRUE, adjust=1.5), lty= 3, lwd= 2, col = "blue")
colours=c("black","green","blue")
legend("topright", 
       legend = c("Neutral","Extreme expansion","Mild expansion"), 
       fill = colours
)

### STATS ###
neutral <- sim$nuc_div
extreme2 <- extremesim2$nuc_div
mild2 <- mildsim2$nuc_div
qqnorm(neutral)
qqnorm(extreme2)
qqnorm(mild2)
sigma.y <- sd(neutral)
sigma.x1 <- sd(extreme2)
sigma.x2 <- sd(mild2)
z.test(extreme2,neutral,"two.sided",sigma.x=sigma.x1,sigma.y = sigma.y)
z.test(mild2,neutral,"two.sided",sigma.x=sigma.x2,sigma.y = sigma.y)


####################################################### SELECTION MODELS ############################################################

### more extreme population expansion - Ne increases from 10,000 to 100,000 individuals by 10kya ###

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

save(Selsims, file = "nuc_div_ExtremeExpansionSims.RData")

labels <- c("(a)","(b)","(c)","(d)","(e)","(f)")
for (a in 1:length(Selsims)) {
  par(mfrow=c(3,2))
  par(mar = c(4.0, 4.0, 2.0, 2.0))
  myfile <- paste0("nuc_div_ExtremeExpansion",myselection[[a]][1],"_",myselection[[a]][2],".pdf")
  for (b in 1:length(mytimes)) {
    mytitle <- paste0((labels[b]))
    plot(density(Selsims[[a]][[b]]$nuc_div,na.rm=TRUE, adjust=1.5), main= mytitle, xlab = "Nucleotide diversity score", xlim = c(-5,650), ylim = c(0.0,0.04))
    lines(density(extremesim2$nuc_div, na.rm = TRUE, adjust = 1.5),lty= 3, lwd= 2, col = "green")
    lines(density(sim$nuc_div, na.rm = TRUE, adjust = 1.5),lty= 3, lwd= 2, col = "red")
  }
  dev.print(device=pdf,file=myfile)
}  


### a less extreme expansion event - Ne increases from 10,000 to 50,000 individuals by 10kya ###

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

save(Selsims, file = "nuc_div_MildExpansionSims.RData")

labels <- c("(a)","(b)","(c)","(d)","(e)","(f)")
for (a in 1:length(Selsims)) {
  par(mfrow=c(3,2))
  par(mar = c(4.0, 4.0, 2.0, 2.0))
  myfile <- paste0("nuc_div_MildExpansion",myselection[[a]][1],"_",myselection[[a]][2],".pdf")
  for (b in 1:length(mytimes)) {
    mytitle <- paste0((labels[b]))
    plot(density(Selsims[[a]][[b]]$nuc_div,na.rm=TRUE, adjust=1.5), main= mytitle, xlab = "Nucleotide diversity score", xlim = c(-5,300), ylim = c(0.0,0.04))
    lines(density(mildsim2$nuc_div, na.rm = TRUE, adjust = 1.5),lty= 3, lwd= 2, col = "blue")
    lines(density(sim$nuc_div, na.rm = TRUE, adjust = 1.5),lty= 3, lwd= 2, col = "red")
  }
  dev.print(device=pdf,file=myfile)
}  



################################# PLOTS ###################################
### EXTREME ###
nvector <- extremesim2$nuc_div
percentile <- quantile(nvector,0.01)

load("nuc_div_ExtremeExpansionSims.RData")

FPRs1 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector = Selsims[[1]][[a]][["nuc_div"]]
  FPRs1[a] = length(vector[vector <= percentile])/length(vector)
}

FPRs2 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector2 = Selsims[[2]][[a]][["nuc_div"]]
  FPRs2[a] = length(vector2[vector2 <= percentile])/length(vector2)
}

save(FPRs1, file = "nuc_divExExpanpwr1.RData")
save(FPRs2, file = "nuc_divExExpanpwr2.RData")

### MILD ###
nvector <- mildsim2$nuc_div
percentile <- quantile(nvector,0.01)

load("nuc_div_MildExpansionSims.RData")

FPRs3 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector3 = Selsims[[1]][[a]][["nuc_div"]]
  FPRs3[a] = length(vector3[vector3 <= percentile])/length(vector3)
}

FPRs4 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector4 = Selsims[[2]][[a]][["nuc_div"]]
  FPRs4[a] = length(vector4[vector4 <= percentile])/length(vector4)
}

save(FPRs3, file = "nuc_divMExpanpwr1.RData")
save(FPRs4, file = "nuc_divMExpanpwr2.RData")

