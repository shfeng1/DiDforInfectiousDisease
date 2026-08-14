rm(list=ls())
here::i_am("1a_Scripts/3b_Misspecify_SEIR_to_SIR.R")
source("./global_options.R")
source("./1a_Scripts/0_Misspecify.R")
output.file <- "./4_Output/misspecify_SEIR_to_SIR.rds"
nsim <- 1000 # 1K run at a time
sim.out <- data.frame()

trans_prob.base2 <- 0.115; trans_prob.base1 <- trans_prob.base2*1.1
sim.param <- expand.grid(eff.multi=c(0.8, 1, 1.1, 1.2))

# split 5K simulations into 5 batches
for (sim_batch in seq(0, 4000, by=1000)) { # seeds set as (j, 1000+j, 2000+j, 3000+j, 4000+j)
  for (j in 1:nrow(sim.param)) {
    print(j)
    set.seed(sim_batch+j, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
    out <- foreach(s = 1:nsim,
                   .combine = "rbind",
                   .errorhandling = "stop") %dopar%
      {
        tryCatch(
          sim_misspecify_SEIR(eff.multi1 = sim.param$eff.multi[j], parallel.id = s),
          
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