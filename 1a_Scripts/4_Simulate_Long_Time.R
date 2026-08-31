rm(list=ls())
here::i_am("1a_Scripts/4_Simulate_Long_Time.R")
source("./global_options.R")
source("./1a_Scripts/0_SEIR.R")
output.file <- "./4_Output/SEIR_long_time_T1.rds"
nsim <- 500

models <- c("growth", "Rt_exposure", "beta_exposure")
sim.param <- expand.grid(trans_prob.base2=1.15/inf_mean, trans_prob.ratio=1.1,
  eff.multi=c(0.9, 0.95, 1), T1.weeks=c(3, 5, 10, 15, 20)) %>%
  mutate(trans_prob.base1=trans_prob.base2*trans_prob.ratio)

# Split 1K simulations into 2 batches of 500.
for (sim_batch in seq(0, 500, by=500)) {
  sim.out <- list()
  for (j in 1:nrow(sim.param)) {
    print(j)
    T1 <- sim.param$T1.weeks[j]*agg
    rng.streams <- make_rng_streams(sim_batch+j, nsim)
    chunks <- split(seq_len(nsim), ceiling(seq_len(nsim)/simulation.chunk.size))
    out <- foreach(chunk=chunks, .combine="rbind", .errorhandling="stop") %dopar% {
      run_pure_simulation_chunk(
        chunk, rng.streams,
        function(s) {
          SEIR_sim(pop.size=pop.size, N=N, N1=N1, seed1=seed1, seed2=seed2,
                   T0=T0, T1=T1, burnin=burnin, inf_mean=inf_mean, delta=delta,
                   trans_prob.base1=sim.param$trans_prob.base1[j],
                   trans_prob.base2=sim.param$trans_prob.base2[j],
                   eff.multi1=sim.param$eff.multi[j], parallel.id=s,
                   calculate_p=FALSE)
        }
      )
    }
    sim.out[[j]] <- out
  }
  sim.out <- rbindlist(sim.out, use.names=TRUE, fill=TRUE) %>%
    data.frame() %>%
    filter(model %in% models)
  rownames(sim.out) <- NULL
  if (file.exists(output.file)) sim.out <- rbind(readRDS(output.file), sim.out)
  saveRDS(sim.out, output.file)
}
