rm(list=ls())
here::i_am("1a_Scripts/2a_Calculate_SIR_True_RR.R")
source("./global_options.R")
source("./1a_Scripts/0_SIR.R")
output.file <- "./4_Output/SIR_RR.rds"
nsim <- 1000 # 1K run at a time
sim.out <- data.frame()

sim.param <- expand.grid(trans_prob.base2=1.15/inf_mean, trans_prob.ratio=c(1, 1.1),
                         eff.multi=c(0.8, 0.9, 0.95, 1, 1.05, 1.1, 1.2)) %>% 
  mutate(trans_prob.base1=trans_prob.base2*trans_prob.ratio)

# split 5K simulations into 5 batches
for (sim_batch in seq(0, 4000, by=1000)) { # seeds set as (j, 1000+j, 2000+j, 3000+j, 4000+j)
  sim.out <- data.frame()
  for (j in 1:nrow(sim.param)) {
    print(j)
    set.seed(sim_batch+j, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
    out <- foreach(s = 1:nsim, 
                   .combine = "rbind",
                   .errorhandling = "stop") %dopar% 
      {
        tryCatch(
          SIR_true_eff(trans_prob.base1 = sim.param$trans_prob.base1[j],
                       trans_prob.base2 = sim.param$trans_prob.base2[j],
                       eff.multi1 = sim.param$eff.multi[j]),
          
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
    rownames(out) <- NULL
    sim.out <- rbind(sim.out, out)
  }
  if (file.exists(output.file)) sim.out <- rbind(readRDS(output.file), sim.out)
  saveRDS(sim.out, output.file)
}
