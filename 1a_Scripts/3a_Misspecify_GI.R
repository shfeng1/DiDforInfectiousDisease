rm(list=ls())
here::i_am("1a_Scripts/3a_Misspecify_GI.R")
source("./global_options.R")
source("./1a_Scripts/0_Misspecify.R")
output.file <- "./4_Output/misspecify_GI.rds"

# translate Geometric(1/inf_mean) to a Discrete Gamma:
mean_true <- std_true <- 1 / -log(1 - 1/inf_mean)
var_true <- std_true^2

# fix baseline transmission rates for misspecification simulations
trans_prob.base2 <- 0.115; trans_prob.base1 <- trans_prob.base2*1.1

# each effect size is repeated 5 times in total, so 
#   1K simulation over this grid gives 5K total simulations per effect size
sim.param <- rbind(expand.grid(mean_spe=c(0.8, 1, 1.2), var_spe=1, eff.multi=c(0.8, 1, 1.1, 1.2)),
                   expand.grid(mean_spe=1, var_spe=c(0.8, 1.2), eff.multi=c(0.8, 1, 1.1, 1.2)))

models <- c("Rt_exposure", "beta_exposure")
sim.out <- list()
for (j in 1:nrow(sim.param)) {
  print(j)
  rng.streams <- make_rng_streams(j, nsim)
  chunks <- split(seq_len(nsim), ceiling(seq_len(nsim)/simulation.chunk.size))
  out <- foreach(chunk=chunks, .combine="rbind", .errorhandling="stop") %dopar% {
    run_simulation_chunk(
      chunk, rng.streams,
      function(s) {
        sim_misspecify_GI(mean_true = mean_true, var_true = var_true, eff.multi1 = sim.param$eff.multi[j],
                          mean_spe = sim.param$mean_spe[j], var_spe = sim.param$var_spe[j], parallel.id = s,
                          calculate_p=FALSE, return_data=TRUE)
      }, models, paste0("misspecify_GI_", j, "_", chunk[1L])
    )
  }
  sim.out[[j]] <- out
}
sim.out <- rbindlist(sim.out, use.names=TRUE, fill=TRUE) %>% data.frame()
rownames(sim.out) <- NULL
saveRDS(sim.out, output.file)
