rm(list=ls())
here::i_am("1a_Scripts/3_Simulate_Long_Time.R")
source("./global_options.R")
source("./1a_Scripts/0_SEIR.R")

N <- 83; N1 <- 15; pop.size <- 2e4; seed <- 60 # Kansas example has 83 counties with 15 treated; average pop size is 23539
T0 <- 5*7; T1 <- 20*7; burnin <- 2*7 # T0+burnin = 7 weeks, T1 = 20 weeks
sim.param <- expand.grid(trans_prob.base2=1.15/inf_days, trans_prob.ratio=1.1,
                         eff.multi=c(0.9, 0.95, 1)) %>%
  mutate(trans_prob.base1=trans_prob.base2*trans_prob.ratio)

nsim <- 100 # Note: this simulation takes a long time run because of the long T1, so I split it up to run 100 simulations at a time.
sim.out <- data.frame()
for (j in 1:nrow(sim.param)) {
  print(j)
  set.seed(j, kind = "L'Ecuyer-CMRG") # seeds were set accordingly as j, 100+j, 200+j, etc.
  out <- foreach(s = 1:nsim,
                 .combine = "rbind",
                 .errorhandling = "remove") %dopar%
    {
      SEIR_sim(pop.size=pop.size, N=N, N1=N1, seed1=seed, seed2=seed,
               T0=T0, T1=T1, burnin=burnin, inf_days=inf_days, delta=delta,
               trans_prob.base1=sim.param$trans_prob.base1[j],
               trans_prob.base2=sim.param$trans_prob.base2[j],
               eff.multi1=sim.param$eff.multi[j], parallel.id=s)
    }
  sim.out <- rbind(sim.out, out)
}
sim.out$inf_days <- inf_days; sim.out$delta <- delta

filename <- "./4_Output/SEIR_long_time.rds"
if (file.exists(filename)) {
  saveRDS(rbind(readRDS(filename), sim.out), filename)
} else {
  saveRDS(sim.out, filename)
}
