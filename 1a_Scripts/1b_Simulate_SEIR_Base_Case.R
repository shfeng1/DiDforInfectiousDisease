rm(list=ls())
here::i_am("1a_Scripts/1b_Simulate_SEIR_Base_Case.R")
source("./global_options.R")
source("./1a_Scripts/0_SEIR.R")

smearing <- FALSE
smearing.method <- "local"
smearing.reps <- 500L

sim.param <- expand.grid(trans_prob.base2=1.15/inf_mean, trans_prob.ratio=1.1,
                         eff.multi=c(0.8, 1, 1.1, 1.2)) %>% # only need to simulate for the extreme ends
  mutate(trans_prob.base1=trans_prob.base2*trans_prob.ratio)

nsim <- 1000
sim.out <- data.frame()
for (j in 1:nrow(sim.param)) {
  print(j)
  # seeds were set as (j, 1000+j, 2000+j, 3000+j, 4000+j)
  set.seed(4000+j, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
  out <- foreach(s = 1:nsim,
                 .combine = "rbind",
                 .errorhandling = "stop") %dopar%
    {
      tryCatch(
        SEIR_sim(pop.size=pop.size, N=N, N1=N1, seed1=seed1, seed2=seed2,
                 T0=T0, T1=T1, burnin=burnin, inf_mean=inf_mean, delta=delta,
                 trans_prob.base1=sim.param$trans_prob.base1[j],
                 trans_prob.base2=sim.param$trans_prob.base2[j],
                 eff.multi1=sim.param$eff.multi[j], parallel.id=s,
                smearing=smearing, smearing_reps=smearing.reps,
                smearing_method=smearing.method),
        
        # Allow epidemic extinction to pass
        epidemic_extinct = function(e) NULL,
        # Skip only the two explicitly allowed transient NFS writer errors.
        # Every other failure now reports the parameter row, replicate, and
        # exact my_RStata stage instead of collapsing to "cannot open the connection".
        error = function(e) handle_simulation_error(e, j, s))
    }
  sim.out <- rbind(sim.out, out)
}

if (nrow(sim.out) == 0L) stop("No SEIR base-case simulations completed successfully.")
output.file <- "./4_Output/SEIR_base_case.rds"
dir.create(dirname(output.file), recursive=TRUE, showWarnings=FALSE)
if (file.exists(output.file)) sim.out <- rbind(readRDS(output.file), sim.out)
saveRDS(sim.out, output.file)
