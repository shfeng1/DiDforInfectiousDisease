rm(list=ls())
here::i_am("1a_Scripts/2b_Calculate_SEIR_True_RR.R")
source("./global_options.R")
source("./1a_Scripts/0_SEIR.R")

sim.param <- expand.grid(trans_prob.base2=1.15/inf_mean, trans_prob.ratio=1.1,
                         eff.multi1=c(0.8, 0.9, 0.95, 1, 1.05, 1.1, 1.2)) %>% 
  mutate(trans_prob.base1=trans_prob.base2*trans_prob.ratio)

nsim <- 1000
sim.out <- data.frame()
for (j in 1:nrow(sim.param)) {
  print(j)
  set.seed(j, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
  out <- foreach(s = 1:nsim, 
                 .combine = "rbind",
                 .errorhandling = "remove") %dopar% 
    {
      SEIR_true_eff(trans_prob.base1 = sim.param$trans_prob.base1[j],
                    trans_prob.base2 = sim.param$trans_prob.base2[j],
                    eff.multi1 = sim.param$eff.multi1[j])
    }
  sim.out <- rbind(sim.out, out)
}
saveRDS(sim.out, "./4_Output/SEIR_RR.rds")

