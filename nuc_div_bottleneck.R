########################################################### SETUP ##################################################################
setwd("~/Desktop/DissertationR")

library(coala)
activate_msms(jar = "~/Desktop/DissertationR/msms/lib/msms.jar", priority = 500, download = FALSE)

list_simulators()
library(BSDA)

####################################################### NEUTRAL MODELS ############################################################ 
load("nuc_div_model.RData")
load("nuc_div_sim.RData")

extrememodel <- model + feat_size_change(0.1, time = 0.005) + feat_size_change(1.0, time = 0)
check_model(extrememodel)
extremesim <- simulate(extrememodel)
save(extremesim, file = "extremesim_nuc_div.RData")

mildmodel <- model + feat_size_change(0.5, time = 0.005) + feat_size_change(1.0, time = 0)
check_model(mildmodel)
mildsim <- simulate(mildmodel)
save(mildsim, file = "mildsim_nuc_div.RData")

par(mfrow=c(1,1))
plot(density(sim$nuc_div,na.rm=TRUE, adjust=1.5), lwd= 2, col = "black", main = "Distributions of nucleotide diversity (π) under population bottleneck scenarios", xlab = "Nucleotide diversity (π)", xlim = c(-5,160), ylim = c(0.0,0.12))
lines(density(extremesim$nuc_div, na.rm=TRUE, adjust=1.5), lty= 3, lwd= 2, col = "green")
lines(density(mildsim$nuc_div, na.rm=TRUE, adjust=1.5), lty= 3, lwd= 2, col = "blue")
colours=c("black","green","blue")
legend("topright", 
       legend = c("Neutral","Extreme bottleneck","Mild bottleneck"), 
       fill = colours
)


### STATS ###
neutral <- sim$nuc_div
extreme <- extremesim$nuc_div
mild <- mildsim$nuc_div
qqnorm(neutral)
qqnorm(extreme)
qqnorm(mild)
sigma.y <- sd(neutral)
sigma.x1 <- sd(extreme)
sigma.x2 <- sd(mild)
z.test(extreme,neutral,"two.sided",sigma.x=sigma.x1,sigma.y = sigma.y)
z.test(mild,neutral,"two.sided",sigma.x=sigma.x2,sigma.y = sigma.y)


####################################################### SELECTION MODELS ############################################################

### extreme bottleneck situation

mytimes <- c(0.005,0.01,0.02,0.05,0.1,0.2)

myselection <- vector(mode="list",length=2)
myselection[[1]] <- c(800,400) 
myselection[[2]] <- c(1200,600)
ExBotSelsims_nuc_div <- vector(mode="list", length=length(myselection))

for (a in 1:length(myselection)) {
  for (b in 1:length(mytimes)) {
    print(paste0("ExBotSelsims_nuc_div",a,"-",b," AA: ",myselection[[a]][1]," Aa: ",myselection[[a]][2]," time: ",mytimes[b]))
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
    ExBotSelsims_nuc_div[[a]][[b]] <- simulate(SelMod1) 
  }
}

save(ExBotSelsims_nuc_div, file = "ExBotSelsims_nuc_div.RData")

labels <- c("(a)","(b)","(c)","(d)","(e)","(f)")
for (a in 1:length(ExBotSelsims_nuc_div)) {
  par(mfrow=c(3,2))
  par(mar = c(4.0, 4.0, 2.0, 2.0))
  myfile <- paste0("nuc_div_Extremebottleneck",myselection[[a]][1],"_",myselection[[a]][2],".pdf")
  for (b in 1:length(mytimes)) {
    mytitle <- paste0((labels[b]))
    plot(density(ExBotSelsims_nuc_div[[a]][[b]]$nuc_div,na.rm=TRUE, adjust=1.5), main= mytitle, xlab = "Nucleotide diversity score", xlim = c(-5,160), ylim = c(0.0,0.5))
    lines(density(extremesim$nuc_div, na.rm=TRUE, adjust=1.5), lty= 3, lwd= 2, col = "green")
    lines(density(sim$nuc_div, na.rm=TRUE, adjust=1.5), lty= 3, lwd= 2, col = "red")
  }
  dev.print(device=pdf,file=myfile)
}  


### a less extreme bottleneck situation - loss of 50% 5kya, regain Ne at time=0
MBotSelsims_nuc_div <- vector(mode="list", length=length(myselection))

for (a in 1:length(myselection)) {
  for (b in 1:length(mytimes)) {
    print(paste0("MBotSelsims_nuc_div",a,"-",b," AA: ",myselection[[a]][1]," Aa: ",myselection[[a]][2]," time: ",mytimes[b]))
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
    MBotSelsims_nuc_div[[a]][[b]] <- simulate(SelMod2) 
  }
}

save(MBotSelsims_nuc_div, file = "MBotSelsims_nuc_div.RData")

labels <- c("(a)","(b)","(c)","(d)","(e)","(f)")
for (a in 1:length(MBotSelsims_nuc_div)) {
  par(mfrow=c(3,2))
  par(mar = c(4.0, 4.0, 2.0, 2.0))
  myfile <- paste0("nuc_div_Mildbottleneck",myselection[[a]][1],"_",myselection[[a]][2],".pdf")
  for (b in 1:length(mytimes)) {
    mytitle <- paste0((labels[b]))
    plot(density(MBotSelsims_nuc_div[[a]][[b]]$nuc_div,na.rm=TRUE, adjust=1.5), main= mytitle, xlab = "Nucleotide diversity score", xlim = c(-5,160), ylim = c(0.0,0.15))
    lines(density(mildsim$nuc_div, na.rm=TRUE, adjust=1.5), lty= 3, lwd= 2, col = "blue")
    lines(density(sim$nuc_div, na.rm=TRUE, adjust=1.5), lty= 3, lwd= 2, col = "red")
  }
  dev.print(device=pdf,file=myfile)
}  


################################# PLOTS ###################################
### EXTREME ###
nvector <- extremesim$nuc_div
percentile <- quantile(nvector,0.01)

FPRs1 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector = ExBotSelsims_nuc_div[[1]][[a]][["nuc_div"]]
  FPRs1[a] = length(vector[vector <= percentile])/length(vector)
}

FPRs2 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector2 = ExBotSelsims_nuc_div[[2]][[a]][["nuc_div"]]
  FPRs2[a] = length(vector2[vector2 <= percentile])/length(vector2)
}

save(FPRs1, file = "nuc_divExBotpwr1.RData")
save(FPRs2, file = "nuc_divExBotpwr2.RData")

### MILD ###
nvector <- mildsim$nuc_div
percentile <- quantile(nvector,0.01)

FPRs3 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector3 = MBotSelsims_nuc_div[[1]][[a]][["nuc_div"]]
  FPRs3[a] = length(vector3[vector3 <= percentile])/length(vector3)
}

FPRs4 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector4 = MBotSelsims_nuc_div[[2]][[a]][["nuc_div"]]
  FPRs4[a] = length(vector4[vector4 <= percentile])/length(vector4)
}

save(FPRs3, file = "nuc_divMBotpwr1.RData")
save(FPRs4, file = "nuc_divMBotpwr2.RData")


