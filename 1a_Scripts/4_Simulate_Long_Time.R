rm(list=ls())
here::i_am("1a_Scripts/4_Simulate_Long_Time.R")
source("./global_options.R")
source("./1a_Scripts/0_SIR.R")
output.file <- "./4_Output/SIR_long_time_T1.rds"

# consider a much larger population to account for the prolonged time periods
pop.size <- 1e6; seed1 <- seed2 <- 1000

sim.param <- expand.grid(trans_prob.base2=1.01/inf_mean, trans_prob.ratio=1,
  eff.multi=1, T1.weeks=c(5, 10, 15, 20)) %>%
  mutate(trans_prob.base1=trans_prob.base2*trans_prob.ratio)

# split 3K simulations into 3 batches
for (sim_batch in seq(0, 2000, by=1000)) {
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
          SIR_sim(pop.size=pop.size, N=N, N1=N1, seed1=seed1, seed2=seed2,
                   T0=T0, T1=T1, burnin=burnin, inf_mean=inf_mean,
                   trans_prob.base1=sim.param$trans_prob.base1[j],
                   trans_prob.base2=sim.param$trans_prob.base2[j],
                   eff.multi1=sim.param$eff.multi[j], parallel.id=s,
                   calculate_p=FALSE)
        }
      )
    }
    sim.out[[j]] <- out
  }
  sim.out <- rbindlist(sim.out, use.names=TRUE, fill=TRUE) %>% data.frame()
  rownames(sim.out) <- NULL
  if (file.exists(output.file)) sim.out <- rbind(readRDS(output.file), sim.out)
  saveRDS(sim.out, output.file)
}
