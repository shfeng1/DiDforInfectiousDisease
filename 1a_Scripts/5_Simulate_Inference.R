rm(list=ls())
here::i_am("1a_Scripts/4_Simulate_Inference.R")
source("./global_options.R")
source("./1a_Scripts/0_SIR.R")

N <- 50; N1 <- 5
sim.param <- expand.grid(trans_prob.base2=1.15/inf_mean, trans_prob.ratio=1.1, eff.multi=c(0.9, 0.95, 1, 1.05, 1.1)) %>%
  mutate(trans_prob.base1=trans_prob.base2*trans_prob.ratio)

nsim <- 5000
sim.out <- data.frame()
for (j in 1:nrow(sim.param)) {
  print(j)
  set.seed(j, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
  out <- foreach(s = 1:nsim,
                 .combine = "rbind",
                 .errorhandling = "remove") %dopar%
    {
      inference_sim(pop.size=pop.size, N=N, N1=N1, seed1=seed1, seed2=seed2,
                    T0=T0, T1=T1, burnin=burnin, inf_mean=inf_mean,
                    trans_prob.base1=sim.param$trans_prob.base1[j],
                    trans_prob.base2=sim.param$trans_prob.base2[j],
                    eff.multi1=sim.param$eff.multi[j], parallel.id=s)
    }
  sim.out <- rbind(sim.out, out)
}

saveRDS(sim.out, "./4_Output/SIR_N1=5.rds")