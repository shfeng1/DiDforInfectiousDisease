rm(list=ls())
here::i_am("1a_Scripts/2a_Calculate_SIR_True_RR.R")
source("./global_options.R")
source("./1a_Scripts/0_SIR.R")

sim.param <- expand.grid(trans_prob.base2=1.15/inf_mean, trans_prob.ratio=c(1, 1.1),
                         eff.multi1=c(0.8, 0.9, 0.95, 1, 1.05, 1.1, 1.2)) %>% 
  mutate(trans_prob.base1=trans_prob.base2*trans_prob.ratio)

nsim <- 1000
sim.out <- data.frame()
for (j in 1:nrow(sim.param)) {
  print(j)
  set.seed(j, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
  out <- foreach(s = 1:nsim, 
                 .combine = "rbind",
                 .errorhandling = "stop") %dopar% 
    {
      # Skip only that understood stochastic outcome; all other errors stop.
      tryCatch(
        SIR_true_eff(trans_prob.base1 = sim.param$trans_prob.base1[j],
                     trans_prob.base2 = sim.param$trans_prob.base2[j],
                     eff.multi1 = sim.param$eff.multi1[j]),
        epidemic_extinct = function(e) NULL,
        error = function(e) handle_simulation_error(e, j, s)
      )
    }
  sim.out <- rbind(sim.out, out)
}

if (nrow(sim.out) == 0L) stop("No SIR true-effect simulations completed successfully.")
output.file <- "./4_Output/SIR_RR.rds"
dir.create(dirname(output.file), recursive=TRUE, showWarnings=FALSE)
if (file.exists(output.file)) sim.out <- rbind(readRDS(output.file), sim.out)
saveRDS(sim.out, output.file)

