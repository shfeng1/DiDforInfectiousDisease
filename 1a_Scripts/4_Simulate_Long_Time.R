rm(list=ls())
here::i_am("1a_Scripts/4_Simulate_Long_Time.R")
source("./global_options.R")
source("./1a_Scripts/0_SEIR.R")
output.file <- "./4_Output/SEIR_long_time.rds"
nsim <- 100 # takes about 3 hours to run 100
sim.out <- data.frame()

N <- 83; N1 <- 15; pop.size <- 2e4; seed <- 60 # Kansas example has 83 counties with 15 treated; average pop size is 23539
T0 <- 5*7; T1 <- 20*7; burnin <- 2*7 # T0+burnin = 7 weeks, T1 = 20 weeks
sim.param <- expand.grid(trans_prob.base2=1.15/inf_mean, trans_prob.ratio=1.1,
                         eff.multi=c(0.9, 0.95, 1)) %>%
  mutate(trans_prob.base1=trans_prob.base2*trans_prob.ratio)

# split 5K simulations into 50 batches
for (sim_batch in seq(100, 5000, by=100)) { # seeds set as (100+j, 200+j, ..., 5000+j)
  for (j in 1:nrow(sim.param)) {
    print(j)
    set.seed(sim_batch+j, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
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
          
          # Allow the known writer error
          error = function(e) {
            msg <- conditionMessage(e)
            
            allowed_write_error <-
              grepl("unable to open file for writing", msg, fixed = TRUE) &&
              (
                grepl("Stale NFS file handle", msg, fixed = TRUE) ||
                  grepl("Operation timed out", msg, fixed = TRUE)
              )
            
            if (allowed_write_error) {
              return(NULL)
            }
            
            # Any other error stops the parallel batch.
            stop(e)
          })
      }
    sim.out <- rbind(sim.out, out)
  }
  if (file.exists(output.file)) sim.out <- rbind(readRDS(output.file), sim.out)
  saveRDS(sim.out, output.file)
}
