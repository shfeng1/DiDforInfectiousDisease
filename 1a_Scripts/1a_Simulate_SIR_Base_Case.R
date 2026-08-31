rm(list=ls())
here::i_am("1a_Scripts/1a_Simulate_SIR_Base_Case.R")
source("./global_options.R")
source("./1a_Scripts/0_SIR.R")
output.file <- "./4_Output/SIR_base_case.rds"

sim.param <- expand.grid(trans_prob.base2=1.15/inf_mean, trans_prob.ratio=c(1, 1.1),
                         eff.multi=c(0.8, 0.9, 0.95, 1, 1.05, 1.1, 1.2)) %>%
  mutate(trans_prob.base1=trans_prob.base2*trans_prob.ratio)

# split 5K simulations into 5 batches
for (sim_batch in seq(0, 4000, by=1000)) { # seeds set as (j, 1000+j, 2000+j, 3000+j, 4000+j)
  sim.out <- list()
  for (j in 1:nrow(sim.param)) {
    print(j)
    rng.streams <- make_rng_streams(sim_batch+j, nsim)
    chunks <- split(seq_len(nsim), ceiling(seq_len(nsim)/simulation.chunk.size))
    out <- foreach(chunk=chunks, .combine="rbind", .errorhandling="stop") %dopar% {
      run_simulation_chunk(
        chunk, rng.streams,
        function(s) {
          SIR_sim(pop.size=pop.size, N=N, N1=N1, seed1=seed1, seed2=seed2,
                  T0=T0, T1=T1, burnin=burnin, inf_mean=inf_mean,
                  trans_prob.base1=sim.param$trans_prob.base1[j],
                  trans_prob.base2=sim.param$trans_prob.base2[j],
                  eff.multi1=sim.param$eff.multi[j], parallel.id=s,
                  calculate_p=FALSE, return_data=TRUE)
        }, model.list, paste0("SIR_", sim_batch, "_", j, "_", chunk[1L])
      )
    }
    sim.out[[j]] <- out
  }
  sim.out <- rbindlist(sim.out, use.names=TRUE, fill=TRUE) %>% data.frame()
  rownames(sim.out) <- NULL
  if (file.exists(output.file)) sim.out <- rbind(readRDS(output.file), sim.out)
  saveRDS(sim.out, output.file)
}
