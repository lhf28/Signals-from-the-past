########################################################### SETUP ##############################################################
setwd("~/Desktop/DissertationR")

library(coala)
activate_msms(jar = "~/Desktop/DissertationR/msms/lib/msms.jar", priority = 500, download = FALSE)

list_simulators()

######################################################## NEUTRAL MODEL #########################################################
model <- coal_model(
  sample_size = 300,
  loci_number = 3000,
  loci_length = 200000,
  ploidy = 1)
model <-
  model + feat_mutation(
    rate = 100,
    model = "IFS",
    base_frequencies = NA,
    tstv_ratio = NA,
    gtr_rates = NA,
    fixed_number = FALSE,
    locus_group = "all"
  ) + feat_recombination(80) + sumstat_nucleotide_div("nuc_div")
check_model(model)
sim <- simulate(model)
summary(sim$nuc_div)

save(model, file = "nuc_div_model.RData")
save(sim, file = "nuc_div_sim.RData")

par(mfrow=c(1,1))
plot(density(sim$nuc_div,na.rm=TRUE, adjust=1.5), main= "Neutral distribution of Nucleotide Diversity", xlab = "Nucleotide diversity (π)")

########################################################### SELECTION ##############################################################
mytimes <- c(0.005,0.01,0.02,0.05,0.1,0.2)

myselection <- vector(mode="list",length=4)
myselection[[1]] <- c(200,100)
myselection[[2]] <- c(800,400) 
myselection[[3]] <- c(1200,600) 
myselection[[4]] <- c(4000,2000) 

Selsims <- vector(mode="list", length=length(myselection))

for (a in 1:length(myselection)) {
  for (b in 1:length(mytimes)) {
  print(paste0("Selsims",a,"-",b," AA: ",myselection[[a]][1]," Aa: ",myselection[[a]][2]," time: ",mytimes[b]))
  SelMod <- model + feat_selection(
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
  check_model(SelMod)
  Selsims[[a]][[b]] <- simulate(SelMod) 
  }
}

save(Selsims, file = "nuc_div_Selsims.RData")

################################################## PLOTS ################################################################
labels <- c("(a)","(b)","(c)","(d)","(e)","(f)")
for (a in 1:length(Selsims)) {
  par(mfrow=c(3,2))
  par(mar = c(4.0, 4.0, 2.0, 2.0))
  myfile <- paste0("Nuc_div_",myselection[[a]][1],"_",myselection[[a]][2],".pdf")
  for (b in 1:length(mytimes)) {
    mytitle <- paste0((labels[b]))
    plot(density(Selsims[[a]][[b]]$nuc_div,na.rm=TRUE, adjust=1.5), main= mytitle, xlab = "Nucleotide diversity (π)", xlim = c(-30,160), ylim = c(0.0,0.05))
    lines(density(sim$nuc_div, na.rm=TRUE, adjust=1.5), lty= 3, lwd= 2, col = "red")
  }
  dev.print(device=pdf,file=myfile)
}  


################################################# STATS #################################################################
nvector <- sim$nuc_div
percentile <- quantile(nvector,0.01)

FPRs5 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector = Selsims[[1]][[a]][["nuc_div"]]
  FPRs5[a] = length(vector[vector <= percentile])/length(vector)
}

FPRs6 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector2 = Selsims[[2]][[a]][["nuc_div"]]
  FPRs6[a] = length(vector2[vector2 <= percentile])/length(vector2)
}

FPRs7 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector3 = Selsims[[3]][[a]][["nuc_div"]]
  FPRs7[a] = length(vector3[vector3 <= percentile])/length(vector3)
}

FPRs8 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector4 = Selsims[[4]][[a]][["nuc_div"]]
  FPRs8[a] = length(vector4[vector4 <= percentile])/length(vector4)
}

save(FPRs6, file = "nuc_divneutralpwr1.RData")
save(FPRs7, file = "nuc_divneutralpwr2.RData")

