rm(list=ls())
here::i_am("1a_Scripts/3a_Misspecify_GI.R")
source("./global_options.R")
source("./1a_Scripts/0_Misspecify.R")

# translate Geometric(1/inf_mean) to a Discrete Gamma:
mean_true <- std_true <- 1 / -log(1 - 1/inf_mean)
var_true <- std_true^2

# fix baseline transmission rates for misspecification simulations
trans_prob.base2 <- 0.115; trans_prob.base1 <- trans_prob.base2*1.1

nsim <- 1000
sim.param <- rbind(expand.grid(mean_spe=c(0.8, 1, 1.2), var_spe=1, eff.multi=c(0.8, 1, 1.2)),
                   expand.grid(mean_spe=1, var_spe=c(0.8, 1.2), eff.multi=c(0.8, 1, 1.2)))

sim.out <- data.frame()
for (j in 1:nrow(sim.param)) {
  print(j)
  set.seed(j, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
  out <- foreach(s = 1:nsim,
                 .combine = "rbind",
                 .errorhandling = "remove") %dopar%
    {
      sim_misspecify_GI(mean_true = mean_true, var_true = var_true,
                        eff.multi1 = sim.param$eff.multi[j], 
                        mean_spe = sim.param$mean_spe[j], 
                        var_spe = sim.param$var_spe[j],
                        parallel.id = s)
    }
  sim.out <- rbind(sim.out, out)
}
saveRDS(sim.out, "./4_Output/misspecify_GI.rds")
