rm(list=ls())
here::i_am("1a_Scripts/2b_Calculate_SEIR_True_RR.R")
source("./global_options.R")
source("./1a_Scripts/0_SEIR.R")
output.file <- "./4_Output/SEIR_RR.rds"
nsim <- 1000 # 1K run at a time.
sim.out <- data.frame()

sim.param <- expand.grid(trans_prob.base2=1.15/inf_mean, trans_prob.ratio=1.1,
                         eff.multi=c(0.8, 0.9, 0.95, 1, 1.05, 1.1, 1.2)) %>% 
  mutate(trans_prob.base1=trans_prob.base2*trans_prob.ratio)

# split 5K simulations into 5 batches
for (sim_batch in seq(0, 4000, by=1000)) { # seeds set as (j, 1000+j, 2000+j, 3000+j, 4000+j)
  for (j in 1:nrow(sim.param)) {
    print(j)
    set.seed(sim_batch+j, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
    out <- foreach(s = 1:nsim, 
                   .combine = "rbind",
                   .errorhandling = "stop") %dopar% 
      {
        SEIR_true_eff(trans_prob.base1 = sim.param$trans_prob.base1[j],
                      trans_prob.base2 = sim.param$trans_prob.base2[j],
                      eff.multi1 = sim.param$eff.multi[j])
      }
    sim.out <- rbind(sim.out, out)
  }
  if (file.exists(output.file)) sim.out <- rbind(readRDS(output.file), sim.out)
  saveRDS(sim.out, output.file)
}
