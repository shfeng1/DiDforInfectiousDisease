rm(list=ls())
here::i_am("1a_Scripts/3b_Misspecify_SEIR_to_SIR.R")
source("./global_options.R")
source("./1a_Scripts/0_Misspecify.R")
output.file <- "./4_Output/misspecify_SEIR_to_SIR.rds"

trans_prob.base2 <- 0.115; trans_prob.base1 <- trans_prob.base2*1.1
sim.param <- expand.grid(eff.multi=c(0.8, 1, 1.1, 1.2))

# split 5K simulations into 5 batches
models <- c("Rt_exposure", "beta_exposure")
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
          sim_misspecify_SEIR(eff.multi1=sim.param$eff.multi[j], parallel.id=s,
                              calculate_p=FALSE, return_data=TRUE)
        }, models, paste0("misspecify_SEIR_", sim_batch, "_", j, "_", chunk[1L])
      )
    }
    sim.out[[j]] <- out
  }
  sim.out <- rbindlist(sim.out, use.names=TRUE, fill=TRUE) %>% data.frame()
  rownames(sim.out) <- NULL
  if (file.exists(output.file)) sim.out <- rbind(readRDS(output.file), sim.out)
  saveRDS(sim.out, output.file)
}
