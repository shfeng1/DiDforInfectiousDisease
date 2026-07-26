rm(list=ls())
here::i_am("1a_Scripts/4_Simulate_Long_Time.R")
source("./global_options.R")
source("./1a_Scripts/0_SEIR.R")

N <- 83; N1 <- 13; pop.size <- 2e4; seed <- 60
T0 <- 5*7; T1 <- 20*7; burnin <- 2*7 # T0+burnin = 7 weeks, T1 = 20 weeks
sim.param <- expand.grid(trans_prob.base2=1.15/inf_mean, trans_prob.ratio=1.1,
                         eff.multi=c(0.9, 0.95, 1)) %>%
  mutate(trans_prob.base1=trans_prob.base2*trans_prob.ratio)

nsim <- 100 # takes about 1.7 hours to run 50
sim.out <- data.frame()
for (j in 1:nrow(sim.param)) {
  print(j)
  set.seed(j, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
  out <- foreach(s = 1:nsim,
                 .combine = "rbind",
                 .errorhandling = "stop") %dopar%
    {
      tryCatch(
        SEIR_sim(pop.size=pop.size, N=N, N1=N1, seed1=seed, seed2=seed,
                 T0=T0, T1=T1, burnin=burnin, inf_mean=inf_mean, delta=delta,
                 trans_prob.base1=sim.param$trans_prob.base1[j],
                 trans_prob.base2=sim.param$trans_prob.base2[j],
                 eff.multi1=sim.param$eff.multi[j], parallel.id=s),
        
        # Allow epidemic extinction to pass
        epidemic_extinct = function(e) NULL,
        error = function(e) handle_simulation_error(e, j, s))
    }
  sim.out <- rbind(sim.out, out)
}
if (nrow(sim.out) == 0L) stop("No long-time SEIR simulations completed successfully.")
sim.out$inf_mean <- inf_mean; sim.out$delta <- delta

output.file <- "./4_Output/SEIR_long_time.rds"
dir.create(dirname(output.file), recursive=TRUE, showWarnings=FALSE)
if (file.exists(output.file)) sim.out <- rbind(readRDS(output.file), sim.out)
saveRDS(sim.out, output.file)
