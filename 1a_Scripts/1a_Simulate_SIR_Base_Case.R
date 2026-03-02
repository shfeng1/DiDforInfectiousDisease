rm(list=ls())
here::i_am("1a_Scripts/1a_Simulate_SIR_Base_Case.R")
source("./global_options.R")
source("./1a_Scripts/0_SIR.R")

sim.param <- expand.grid(trans_prob.base2=1.15/inf_mean, trans_prob.ratio=c(1, 1.1),
                         eff.multi=c(0.8, 0.9, 0.95, 1, 1.05, 1.1, 1.2)) %>%
  mutate(trans_prob.base1=trans_prob.base2*trans_prob.ratio)

nsim <- 1000 # 1K run at a time. With a more powerful machine, this can be directly set to 5000.
sim.out <- data.frame()
for (j in 1:nrow(sim.param)) {
  print(j)
  # seeds were set as (j, 1000+j, 2000+j, 3000+j, 4000+j) # I split 5K simulations into 5 runs on my local machine
  set.seed(j, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
  out <- foreach(s = 1:nsim,
                 .combine = "rbind",
                 .errorhandling = "remove") %dopar%
    {
      SIR_sim(pop.size=pop.size, N=N, N1=N1, seed1=seed1, seed2=seed2,
              T0=T0, T1=T1, burnin=burnin, inf_mean=inf_mean,
              trans_prob.base1=sim.param$trans_prob.base1[j],
              trans_prob.base2=sim.param$trans_prob.base2[j],
              eff.multi1=sim.param$eff.multi[j], parallel.id=s)
    }
  sim.out <- rbind(sim.out, out)
}

# saveRDS(sim.out, "./4_Output/SIR_base_case.rds")
# saveRDS(rbind(sim.out, readRDS("./4_Output/SIR_base_case.rds")), "./4_Output/SIR_base_case.rds")
