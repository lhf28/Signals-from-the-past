########################################################### SETUP ##############################################################
setwd("~/Desktop/DissertationR")

library(coala)
activate_msms(jar = "~/Desktop/DissertationR/msms/lib/msms.jar", priority = 500, download = FALSE)

list_simulators()
library(BSDA)

######################################################## NEUTRAL MODEL #########################################################
# 5 populations l-r African, European, E.Asian, Mexican Americans, neanderthals
OoAmodel <- coal_model(
  sample_size = c(100,100,100,100,100),
  loci_number = 1000,
  loci_length = 200000,
  ploidy = 1)
# M-A and E.Asians split 22kya, European and E.Asian split 23kya, Eurasians and (w).Africans split 51kya + neanderthals merge with pop.1 500kya prior to sweep
OoAmodel <-
  OoAmodel + feat_pop_merge(time = 0.022,
                            pop_source = 4,
                            pop_target = 3) + feat_pop_merge(time = 0.023,
                                                             pop_source = 3,
                                                             pop_target = 2) + feat_pop_merge(time = 0.051,
                                                                                              pop_source = 2,
                                                                                              pop_target = 1) + feat_pop_merge (time = 0.5,
                                                                                                                                pop_source = 5,
                                                                                                                                pop_target = 1)
### migration between europeans and east asians
OoAmodel <- OoAmodel + feat_migration(
  rate = 1.244,
  pop_from = 2,
  pop_to = 3,
  symmetric = FALSE,
  time = 0,
  locus_group = "all"
) + feat_migration(
  rate = 1.244,
  pop_from = 3,
  pop_to = 2,
  symmetric = FALSE,
  time = 0,
  locus_group = "all"
) + feat_migration(
  rate = 0,
  pop_from = 2,
  pop_to = 3,
  symmetric = FALSE,
  time = 0.02299,
  locus_group = "all"
) + feat_migration(
  rate = 0,
  pop_from = 3,
  pop_to = 2,
  symmetric = FALSE,
  time = 0.02299,
  locus_group = "all"
) 
### migration between africans and east asians
OoAmodel <- OoAmodel + feat_migration(
  rate = 0.312,
  pop_from = 1,
  pop_to = 3,
  symmetric = FALSE,
  time = 0,
  locus_group = "all"
) + feat_migration(
  rate = 0.312,
  pop_from = 3,
  pop_to = 1,
  symmetric = FALSE,
  time = 0,
  locus_group = "all"
) + feat_migration(
  rate = 0,
  pop_from = 1,
  pop_to = 3,
  symmetric = FALSE,
  time = 0.02299,
  locus_group = "all"
) + feat_migration(
  rate = 0,
  pop_from = 3,
  pop_to = 1,
  symmetric = FALSE,
  time = 0.02299,
  locus_group = "all"
) 

### migration between africans and europeans
OoAmodel <- OoAmodel + feat_migration(
  rate = 1,
  pop_from = 2,
  pop_to = 1,
  symmetric = FALSE,
  time = 0,
  locus_group = "all"
) + feat_migration(
  rate = 1,
  pop_from = 1,
  pop_to = 2,
  symmetric = FALSE,
  time = 0,
  locus_group = "all"
) + feat_migration(
  rate = 0,
  pop_from = 2,
  pop_to = 1,
  symmetric = FALSE,
  time = 0.02299,
  locus_group = "all"
) + feat_migration(
  rate = 0,
  pop_from = 1,
  pop_to = 2,
  symmetric = FALSE,
  time = 0.02299,
  locus_group = "all"
)

### migration between africans and eurasians 
OoAmodel + feat_migration(
  rate = 6,
  pop_from = 2,
  pop_to = 1,
  symmetric = FALSE,
  time = 0.023001,
  locus_group = "all"
) + feat_migration(
  rate = 6,
  pop_from = 1,
  pop_to = 2,
  symmetric = FALSE,
  time = 0.023001,
  locus_group = "all"
) + feat_migration(
  rate = 0,
  pop_from = 2,
  pop_to = 1,
  symmetric = FALSE,
  time = 0.050999,
  locus_group = "all"
) + feat_migration(
  rate = 0,
  pop_from = 1,
  pop_to = 2,
  symmetric = FALSE,
  time = 0.050999,
  locus_group = "all"
)

### POPULATION GROWTH RATE ### is this backwards?
OoAmodel <-
  OoAmodel + feat_growth(
    rate = 0.0038,
    population = 2,
    time = 0,
    locus_group = "all"
  ) + feat_growth(
    rate = 0,
    population = 2,
    time = 0.02299,
    locus_group = "all"
  ) + feat_growth(
    rate = 0.0048,
    population = 3,
    time = 0,
    locus_group = "all"
  ) + feat_growth(
    rate = 0,
    population = 3,
    time = 0.02299,
    locus_group = "all"
  ) + feat_growth(
    rate = 0.005,
    population = 4,
    time = 0,
    locus_group = "all"
  ) + feat_growth(
    rate = 0,
    population = 4,
    time = 0.02199,
    locus_group = "all"
  )


### USING N ESTIMATES FROM GRAVEL AND GUTENKUNST
OoAmodel <-
  OoAmodel + feat_size_change(1.4474,
                              population = 1,
                              time = 0,
                              locus_group = "all") + feat_size_change(0.08,
                                                                      population = 4,
                                                                      time = 0.022,
                                                                      locus_group = "all") + feat_size_change(0.0554,
                                                                                                              population = 3,
                                                                                                              time = 0.023,
                                                                                                              locus_group = "all") + feat_size_change(0.1032,
                                                                                                                                                      population = 2,
                                                                                                                                                      time = 0.023,
                                                                                                                                                      locus_group = "all") + feat_size_change(0.1861,
                                                                                                                                                                                              population = 2,
                                                                                                                                                                                              time = 0.051,
                                                                                                                                                                                              locus_group = "all") + feat_size_change(0.73,
                                                                                                                                                                                                                                      population = 1,
                                                                                                                                                                                                                                      time = 0.148,
                                                                                                                                                                                                                                      locus_group = "all")

TajOoAmodel <-
  OoAmodel + feat_mutation(
    rate = 100,
    model = "IFS",
    base_frequencies = NA,
    tstv_ratio = NA,
    gtr_rates = NA,
    fixed_number = FALSE,
    locus_group = "all"
  ) + feat_recombination(80) + sumstat_tajimas_d("taji_d")

nuc_divOoAmodel <-
  OoAmodel + feat_mutation(
    rate = 100,
    model = "IFS",
    base_frequencies = NA,
    tstv_ratio = NA,
    gtr_rates = NA,
    fixed_number = FALSE,
    locus_group = "all"
  ) + feat_recombination(80) + sumstat_nucleotide_div("nuc_div")

check_model(TajOoAmodel)
save(TajOoAmodel, file = "TajOoAmodel.RData")
TajOoAsim <- simulate(TajOoAmodel)
check_model(nuc_divOoAmodel)
save(nuc_divOoAmodel, file = "nuc_divOoAmodel.RData")
nuc_divOoAsim <- simulate(nuc_divOoAmodel)

save(TajOoAsim, file = "TajOoAsim.RData")
save(nuc_divOoAsim, file = "nuc_divOoAsim.RData")

load("TajneutralNEW.RData")
plot(density(TajOoAsim$taji_d), xlab = "Tajima's D value", main = "Distribution of Tajima's D under the OoA model", col = "purple", lty = 3, lwd= 2, xlim = c(-2.0,1.5))
lines(density(sim$taji_d,na.rm=TRUE, adjust=1.5), lwd= 2, col = "black")
load("nuc_divneutralNEW.RData")
plot(density(nuc_divOoAsim$nuc_div), xlab = "Nucleotide diversity (π)", main = "Distribution of nucleotide diversity under the OoA model", col = "purple", lty = 3, lwd= 2, xlim = c(0,150))
lines(density(sim2$nuc_div,na.rm=TRUE, adjust=1.5), lwd= 2, col = "black")

neutralTaj <- sim$taji_d
neutralnuc_div <- sim2$nuc_div
OoATaj <- TajOoAsim$taji_d
OoAnuc_div <- nuc_divOoAsim$nuc_div
qqnorm(neutralTaj)
qqnorm(neutralnuc_div)
qqnorm(OoATaj)
qqnorm(OoAnuc_div)
sigma.y <- sd(neutralTaj)
sigma.x <- sd(OoATaj)
z.test(OoATaj,neutralTaj,"two.sided",sigma.x=sigma.x,sigma.y = sigma.y)
sigma.y <- sd(neutralnuc_div)
sigma.x <- sd(OoAnuc_div)
z.test(OoAnuc_div,neutralnuc_div,"two.sided",sigma.x=sigma.x,sigma.y = sigma.y)



### add selection
mytimes <- c(0.005,0.01,0.02,0.05,0.1,0.2)

myselection <- vector(mode="list",length=2)
myselection[[1]] <- c(800,400) 
myselection[[2]] <- c(1200,600) 

Selsims1 <- vector(mode="list", length=length(myselection))

for (a in 1:length(myselection)) {
  for (b in 1:length(mytimes)) {
    print(paste0("TajOoASelsims",a,"-",b," AA: ",myselection[[a]][1]," Aa: ",myselection[[a]][2]," time: ",mytimes[b]))
    SelMod1 <- TajOoAmodel + feat_selection(
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
    Selsims1[[a]][[b]] <- simulate(SelMod1) 
  }
}

save(Selsims1, file = "Taj_Selsims_OoA.RData")

labels <- c("(a)","(b)","(c)","(d)","(e)","(f)")
for (a in 1:length(Selsims1)) {
  par(mfrow=c(3,2))
  par(mar = c(4.0, 4.0, 2.0, 2.0))
  myfile <- paste0("TajimasD_OoA",myselection[[a]][1],"_",myselection[[a]][2],".pdf")
  for (b in 1:length(mytimes)) {
    mytitle <- paste0((labels[b]))
    plot(density(Selsims1[[a]][[b]]$taji_d,na.rm=TRUE, adjust=1.5), main= mytitle, xlab = "Tajima's D value", xlim = c(-3,3), ylim = c(0.0,3.0))
    lines(density(TajOoAsim$taji_d, na.rm=TRUE, adjust=1.5), lty= 3, lwd= 2, col = "red")
  }
  dev.print(device=pdf,file=myfile)
}



Selsims2 <- vector(mode="list", length=length(myselection))

for (a in 1:length(myselection)) {
  for (b in 1:length(mytimes)) {
    print(paste0("nuc_div_Selsims2",a,"-",b," AA: ",myselection[[a]][1]," Aa: ",myselection[[a]][2]," time: ",mytimes[b]))
    SelMod2 <- nuc_divOoAmodel + feat_selection(
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
    Selsims2[[a]][[b]] <- simulate(SelMod2) 
  }
}

save(Selsims2, file = "nuc_div_Selsims_OoA.RData")

for (a in 1:length(Selsims2)) {
  par(mfrow=c(3,2))
  par(mar = c(4.0, 4.0, 2.0, 2.0))
  myfile <- paste0("nuc_div_OoA",myselection[[a]][1],"_",myselection[[a]][2],".pdf")
  for (b in 1:length(mytimes)) {
    mytitle <- paste0((labels[b]))
    plot(density(Selsims2[[a]][[b]]$nuc_div,na.rm=TRUE, adjust=1.5), main= mytitle, xlab = "Nucleotide diversity (π)", ylim = c(0.0,0.040))
    lines(density(nuc_divOoAsim$nuc_div, na.rm=TRUE, adjust=1.5), lty= 3, lwd= 2, col = "red")
  }
  dev.print(device=pdf,file=myfile)
}


TajDnvector <- TajOoAsim$taji_d
Tajpercentile <- quantile(TajDnvector,0.01)
Nuc_divnvector <- nuc_divOoAsim$nuc_div
Nuc_divpercentile <- quantile(Nuc_divnvector,0.01)

FPRs1 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector1 = Selsims1[[1]][[a]][["taji_d"]]
  FPRs1[a] = length(vector1[vector1 <= Tajpercentile])/length(vector1)
}

FPRs2 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector2 = Selsims1[[2]][[a]][["taji_d"]]
  FPRs2[a] = length(vector2[vector2 <= Tajpercentile])/length(vector2)
}

FPRs3 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector3 = Selsims2[[1]][[a]][["nuc_div"]]
  FPRs3[a] = length(vector3[vector3 <= Nuc_divpercentile])/length(vector3)
}

FPRs4 <- vector(mode="list", length=length(mytimes))

for (a in 1:length(mytimes)) {
  vector4 = Selsims2[[2]][[a]][["nuc_div"]]
  FPRs4[a] = length(vector4[vector4 <= Nuc_divpercentile])/length(vector4)
}

save(FPRs1, file = "OoATajpwr1.RData")
save(FPRs2, file = "OoATajpwr2.RData")
save(FPRs3, file = "OoAnuc_divpwr1.RData")
save(FPRs4, file = "OoAnuc_divpwr2.RData")
