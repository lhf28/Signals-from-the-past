########################################################### SETUP ###################################################################
setwd("~/Desktop/DissertationR")

library(coala)
activate_msms(jar = "~/Desktop/DissertationR/msms/lib/msms.jar", priority = 500, download = FALSE)

list_simulators()
library(BSDA)
####################################################### SAMPLE FROM POPULATIONS #####################################################

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
  ) + feat_recombination(80)

TajAfrica <- TajOoAmodel + sumstat_tajimas_d("taji_d", population = 1)
simTajAfrica <- simulate(TajAfrica)
save(simTajAfrica, file = "simTajAfrica.RData")
TajEurope <- TajOoAmodel + sumstat_tajimas_d("taji_d", population = 2)
simTajEurope <- simulate(TajEurope)
save(simTajEurope, file = "simTajEurope.RData")
TajAsia <- TajOoAmodel + sumstat_tajimas_d("taji_d", population = 3)
simTajAsia <- simulate(TajAsia)
save(simTajAsia, file = "simTajAsia.RData")
TajMexAm <- TajOoAmodel + sumstat_tajimas_d("taji_d", population = 4)
simTajMexAm <- simulate(TajMexAm)
save(simTajMexAm, file = "simTajMexAm.RData")

load("TajOoAsim.RData")
load("Taj_neutral_sim_modified.RData")
plot(density(TajOoAsim$taji_d), xlab = "Tajima's D value", main = "Distributions of Tajima's D in pooled and individual samples", lty = 2, lwd= 2, xlim = c(-3.5,1.5), ylim = c(0.0,3.2))
lines(density(sim$taji_d), lty = 1, lwd= 2, col = "black")
lines(density(simTajAfrica$taji_d), lty = 3, lwd= 2, col = "red")
lines(density(simTajEurope$taji_d), lty = 3, lwd= 2, col = "blue")
lines(density(simTajAsia$taji_d), lty = 3, lwd= 2, col = "green")
lines(density(simTajMexAm$taji_d), lty = 3, lwd= 2, col = "orange")
legend("topright", c("Neutral", "Pooled", "African", "European", "Asian", "Mexican American"), lty = c(1,2,3,3,3,3), col = c("black", "black", "red", "blue", "green", "orange"), cex = 0.9, lwd= 2)
mean(TajOoAsim$taji_d)
mean(simTajAfrica$taji_d)
mean(simTajEurope$taji_d)
mean(simTajAsia$taji_d)
mean(simTajMexAm$taji_d)
### STATS ###
neutral <- sim$taji_d
Pooled <- TajOoAsim$taji_d
Africa <- simTajAfrica$taji_d
Europe <- simTajEurope$taji_d
Asia <- simTajAsia$taji_d
MexAm <- simTajMexAm$taji_d
qqnorm(neutral)
qqnorm(Pooled)
qqnorm(Africa)
qqnorm(Europe)
qqnorm(Asia)
qqnorm(MexAm)
sigma.y <- sd(neutral)
sigma.x1 <- sd(Africa)
sigma.x2 <- sd(Europe)
sigma.x3 <- sd(Asia)
sigma.x4 <- sd(MexAm)
sigma.x5 <- sd(Pooled)
z.test(Pooled,neutral,"two.sided",sigma.x=sigma.x5,sigma.y = sigma.y)
z.test(Africa,neutral,"two.sided",sigma.x=sigma.x1,sigma.y = sigma.y)
z.test(Europe,neutral,"two.sided",sigma.x=sigma.x2,sigma.y = sigma.y)
z.test(Asia,neutral,"two.sided",sigma.x=sigma.x3,sigma.y = sigma.y)
z.test(MexAm,neutral,"two.sided",sigma.x=sigma.x4,sigma.y = sigma.y)

nuc_divOoAmodel <-
  OoAmodel + feat_mutation(
    rate = 100,
    model = "IFS",
    base_frequencies = NA,
    tstv_ratio = NA,
    gtr_rates = NA,
    fixed_number = FALSE,
    locus_group = "all"
  ) + feat_recombination(80)

nuc_divAfrica <- nuc_divOoAmodel + sumstat_nucleotide_div("nuc_div", population = 1)
simnuc_divAfrica <- simulate(nuc_divAfrica)
save(simnuc_divAfrica, file = "simnuc_divAfrica.RData")
nuc_divEurope <- nuc_divOoAmodel + sumstat_nucleotide_div("nuc_div", population = 2)
simnuc_divEurope <- simulate(nuc_divEurope)
save(simnuc_divEurope, file = "simnuc_divEurope.RData")
nuc_divAsia <- nuc_divOoAmodel + sumstat_nucleotide_div("nuc_div", population = 3)
simnuc_divAsia <- simulate(nuc_divAsia)
save(simnuc_divAsia, file = "simnuc_divAsia.RData")
nuc_divMexAm <- nuc_divOoAmodel + sumstat_nucleotide_div("nuc_div", population = 4)
simnuc_divMexAm <- simulate(nuc_divMexAm)
save(simnuc_divMexAm, file = "simnuc_divMexAm.RData")

load("nuc_divOoAsim.RData")
load("nuc_div_sim.RData")
plot(density(nuc_divOoAsim$nuc_div), xlab = "Nucleotide diversity (π)", main = "Distributions of Nucleotide diversity (π) in pooled and individual samples", lty = 2, lwd= 2, xlim = c(-50,250), ylim = c(0.0,0.035))
lines(density(sim$nuc_div), lty = 1, lwd= 2, col = "black")
lines(density(simnuc_divAfrica$nuc_div), lty = 3, lwd= 2, col = "red")
lines(density(simnuc_divEurope$nuc_div), lty = 3, lwd= 2, col = "blue")
lines(density(simnuc_divAsia$nuc_div), lty = 3, lwd= 2, col = "green")
lines(density(simnuc_divMexAm$nuc_div), lty = 3, lwd= 2, col = "orange")
legend("topright", c("Neutral", "Pooled", "African", "European", "Asian", "Mexican American"), lty = c(1,2,3,3,3,3), col = c("black", "black", "red", "blue", "green", "orange"), cex = 0.9, lwd= 2)
mean(nuc_divOoAsim$nuc_div)
mean(simnuc_divAfrica$nuc_div)
mean(simnuc_divEurope$nuc_div)
mean(simnuc_divAsia$nuc_div)
mean(simnuc_divMexAm$nuc_div)
### STATS ###
neutral <- sim$nuc_div
Pooled <- nuc_divOoAsim$nuc_div
Africa <- simnuc_divAfrica$nuc_div
Europe <- simnuc_divEurope$nuc_div
Asia <- simnuc_divAsia$nuc_div
MexAm <- simnuc_divMexAm$nuc_div
qqnorm(neutral)
qqnorm(Pooled)
qqnorm(Africa)
qqnorm(Europe)
qqnorm(Asia)
qqnorm(MexAm)
sigma.y <- sd(neutral)
sigma.x1 <- sd(Africa)
sigma.x2 <- sd(Europe)
sigma.x3 <- sd(Asia)
sigma.x4 <- sd(MexAm)
sigma.x5 <- sd(Pooled)
z.test(Pooled,neutral,"two.sided",sigma.x=sigma.x5,sigma.y = sigma.y)
z.test(Africa,neutral,"two.sided",sigma.x=sigma.x1,sigma.y = sigma.y)
z.test(Europe,neutral,"two.sided",sigma.x=sigma.x2,sigma.y = sigma.y)
z.test(Asia,neutral,"two.sided",sigma.x=sigma.x3,sigma.y = sigma.y)
z.test(MexAm,neutral,"two.sided",sigma.x=sigma.x4,sigma.y = sigma.y)

####################################################### SELECTION #####################################################################################

mytimes <- c(0.005,0.01,0.02,0.05,0.1,0.2)

myselection <- vector(mode="list",length=2)
myselection[[1]] <- c(800,400) 
myselection[[2]] <- c(1200,600) 

Selsims1 <- vector(mode="list", length=length(myselection))

for (a in 1:length(myselection)) {
  for (b in 1:length(mytimes)) {
    print(paste0("selsimTajAfrica",a,"-",b," AA: ",myselection[[a]][1]," Aa: ",myselection[[a]][2]," time: ",mytimes[b]))
    SelMod1 <- TajAfrica + feat_selection(
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
save(Selsims1, file = "Selsims1OoA.RData")

Selsims2 <- vector(mode="list", length=length(myselection))

for (a in 1:length(myselection)) {
  for (b in 1:length(mytimes)) {
    print(paste0("selsimTajEurope",a,"-",b," AA: ",myselection[[a]][1]," Aa: ",myselection[[a]][2]," time: ",mytimes[b]))
    SelMod2 <- TajEurope + feat_selection(
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
save(Selsims2, file = "Selsims2OoA.RData")

Selsims3 <- vector(mode="list", length=length(myselection))

for (a in 1:length(myselection)) {
  for (b in 1:length(mytimes)) {
    print(paste0("selsimTajAsia",a,"-",b," AA: ",myselection[[a]][1]," Aa: ",myselection[[a]][2]," time: ",mytimes[b]))
    SelMod3 <- TajAsia + feat_selection(
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
    check_model(SelMod3)
    Selsims3[[a]][[b]] <- simulate(SelMod3) 
  }
}
save(Selsims3, file = "Selsims3OoA.RData")

Selsims4 <- vector(mode="list", length=length(myselection))

for (a in 1:length(myselection)) {
  for (b in 1:length(mytimes)) {
    print(paste0("selsimTajMexAm",a,"-",b," AA: ",myselection[[a]][1]," Aa: ",myselection[[a]][2]," time: ",mytimes[b]))
    SelMod4 <- TajMexAm + feat_selection(
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
    check_model(SelMod4)
    Selsims4[[a]][[b]] <- simulate(SelMod4) 
  }
}
save(Selsims4, file = "Selsims4OoA.RData")

Selsims5 <- vector(mode="list", length=length(myselection))

for (a in 1:length(myselection)) {
  for (b in 1:length(mytimes)) {
    print(paste0("selsimnuc_divAfrica",a,"-",b," AA: ",myselection[[a]][1]," Aa: ",myselection[[a]][2]," time: ",mytimes[b]))
    SelMod5 <- nuc_divAfrica + feat_selection(
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
    check_model(SelMod5)
    Selsims5[[a]][[b]] <- simulate(SelMod5) 
  }
}
save(Selsims5, file = "Selsims5OoA.RData")

Selsims6 <- vector(mode="list", length=length(myselection))

for (a in 1:length(myselection)) {
  for (b in 1:length(mytimes)) {
    print(paste0("selsimnuc_divEurope",a,"-",b," AA: ",myselection[[a]][1]," Aa: ",myselection[[a]][2]," time: ",mytimes[b]))
    SelMod6 <- nuc_divEurope + feat_selection(
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
    check_model(SelMod6)
    Selsims6[[a]][[b]] <- simulate(SelMod6) 
  }
}
save(Selsims6, file = "Selsims6OoA.RData")

Selsims7 <- vector(mode="list", length=length(myselection))

for (a in 1:length(myselection)) {
  for (b in 1:length(mytimes)) {
    print(paste0("selsimnuc_divAsia",a,"-",b," AA: ",myselection[[a]][1]," Aa: ",myselection[[a]][2]," time: ",mytimes[b]))
    SelMod7 <- nuc_divAsia + feat_selection(
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
    check_model(SelMod7)
    Selsims7[[a]][[b]] <- simulate(SelMod7) 
  }
}
save(Selsims7, file = "Selsims7OoA.RData")

Selsims8 <- vector(mode="list", length=length(myselection))

for (a in 1:length(myselection)) {
  for (b in 1:length(mytimes)) {
    print(paste0("selsimnuc_divMexAm",a,"-",b," AA: ",myselection[[a]][1]," Aa: ",myselection[[a]][2]," time: ",mytimes[b]))
    SelMod8 <- nuc_divMexAm + feat_selection(
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
    check_model(SelMod8)
    Selsims8[[a]][[b]] <- simulate(SelMod8) 
  }
}
save(Selsims8, file = "Selsims8OoA.RData")
