rm(list=ls())
here::i_am("1a_Scripts/3b_Misspecify_SEIR_to_SIR.R")
source("./global_options.R")
source("./1a_Scripts/0_Misspecify.R")

# fix baseline transmission rates for misspecification simulations
trans_prob.base2 <- 0.115; trans_prob.base1 <- trans_prob.base2*1.1

nsim <- 1000
sim.param <- expand.grid(eff.multi=c(0.8, 1, 1.1, 1.2))

sim.out <- data.frame()
for (j in 1:nrow(sim.param)) {
  print(j)
  set.seed(1000+j, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
  out <- foreach(s = 1:nsim,
                 .combine = "rbind",
                 .errorhandling = "remove") %dopar%
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

# saveRDS(sim.out, "./4_Output/misspecify_SEIR_to_SIR_inf_days.rds")
saveRDS(rbind(sim.out, readRDS("./4_Output/misspecify_SEIR_to_SIR_inf_days.rds")), "./4_Output/misspecify_SEIR_to_SIR_inf_days.rds")
