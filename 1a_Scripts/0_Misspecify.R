source("./1a_Scripts/0_SIR.R")
source("./1a_Scripts/0_SEIR.R")

# mean_spe, var_spe are the multiplicative factor on the misspecification
# they range from (0.8, 1.2)
sim_misspecify_GI <- function(mean_true, var_true, eff.multi1, mean_spe, var_spe, parallel.id=0,
                              end_buffer=NULL, calculate_p=TRUE,
                              return_data=FALSE) {
  if (is.null(end_buffer)) end_buffer <- get("end_buffer", envir=.GlobalEnv)
  parallel.id <- paste0("misspecify_GI_", parallel.id)
  trt.IDs <- 1:N1
  total.days <- T0+T1+burnin+end_buffer
  out.sim <- lapply(1:N, function(ind) { 
    if (ind %in% trt.IDs) {
      run_SIR_varying(pop.size=pop.size, seeds=seed1, time_steps=total.days, inf_mean=mean_true, inf_var=var_true,
                      trans_prob=c(rep(trans_prob.base1, T0+burnin),
                                   rep(trans_prob.base1*eff.multi1, T1+end_buffer)))
    } else {
      run_SIR_varying(pop.size=pop.size, seeds=seed2, time_steps=total.days, inf_mean=mean_true, inf_var=var_true,
                      trans_prob=rep(trans_prob.base2, total.days))
    }})
  out.df <- rbindlist(out.sim) %>%
    mutate(unit=rep(1:N, each=total.days),
           initial_seed=ifelse(unit %in% trt.IDs, seed1, seed2))
  
  # Misspecification: (mean_spe = var_spe = 1) corresponds to the base case of NO MISPECIFICATION
  inf_mean_spe <- mean_true * mean_spe
  inf_var_spe  <- var_true * var_spe
  # Convert (inf_mean_spe, inf_var_spe) into a Gamma(shape, scale):
  shape_spe <- inf_mean_spe^2 / inf_var_spe
  scale_spe <- inf_var_spe / inf_mean_spe
  
  # Estimate prevalence using Gamma. Prepend the initial seed at t=0 before
  # evaluating the incidence-history convolution, consistent with all other
  # simulation-based prevalence/Rt/beta calculations.
  survival <- 1-pgamma(0:total.days, shape=shape_spe, scale=scale_spe)
  out.df$prevalence_gamma <- NA
  for (i in 1:N) {
    idx <- which(out.df$unit == i)
    inc_i <- out.df$inc[idx]
    seed_i <- unique(out.df$initial_seed[idx])
    if (length(seed_i) != 1L || !is.finite(seed_i) || seed_i < 0) {
      stop("Each simulated unit must have one finite, nonnegative initial seed.")
    }
    inc_history <- c(seed_i, inc_i)
    for (t in seq_len(total.days)) {
      lags <- 0:t
      out.df$prevalence_gamma[idx[t]] <-
        sum(survival[lags+1L] * inc_history[(t+1L)-lags])
    }
  }
  
  df.agg <- out.df %>% 
    group_by(unit) %>%
    arrange(t) %>%
    mutate(week = ceiling(t/agg),
           S_lag = ifelse(t==1, pop.size-initial_seed, lag(S, 1)),
           prevalence_lag = lag(prevalence_gamma, 1),
           inc_lag = lag(inc, 1),
           I_lag = ifelse(t==1, initial_seed, lag(I, 1))) %>%
    filter(t >= 3) %>%
    group_by(unit, week) %>%
    summarise(inc = sum(inc),
              growth = sum(inc) / sum(inc_lag),
              S_frac = sum(S_lag) / (pop.size*agg),
              R_true = sum(inc) / sum(I_lag),
              Rt_exposure = sum(inc) / sum(prevalence_lag),
              initial_seed = first(initial_seed)) %>%
    group_by(unit) %>%
    mutate(beta_true = R_true / S_frac,
           beta_exposure = Rt_exposure / S_frac,
           trt.unit = (unit %in% trt.IDs),
           trt.time = (week > (T0+burnin)/agg),
           trt_post = (trt.unit & trt.time),
           time_to_trt = ifelse(!trt_post, 0, week - (T0+burnin)/agg))
  
  data.in <- df.agg %>% data.frame() %>%
    mutate(unit = relevel(factor(unit), ref = N),
           week = week - burnin/agg) %>%
    filter(week > 0, week <= max(week) - end_buffer/agg)
  ############################################################################################################################
  # Fit regression models
  Rt_exposure.out <- run_Rt(data.in, out.df, dgp="SIR", inf_mean=inf_mean_spe, inf_var=inf_var_spe,
                            parallel.id=parallel.id, calculate_p=calculate_p)
  beta.out <- run_beta(data.in, out.df, dgp="SIR", inf_mean=inf_mean_spe, inf_var=inf_var_spe,
                       parallel.id=parallel.id, calculate_p=calculate_p)
  Y.untrt.true <- run_true(out.df, trans_prob.base1, dgp="SIR")
  ################################################################################################################
  # summarize outputs
  out <- rbind(Rt_exposure.out, beta.out) %>%
    mutate(eff.multi=eff.multi1, seed=seed1, mean_spe=mean_spe, var_spe=var_spe, 
           inf_mean=mean_true, S_frac.mean=mean(data.in$S_frac[data.in$trt.time]), S_frac.min=min(data.in$S_frac),
           Y.trt=mean(data.in$inc[data.in$trt_post]), Y.untrt.true=Y.untrt.true)
  if (return_data) return(list(result=out, bootstrap_data=data.in))
  out
}

sim_misspecify_SEIR <- function(eff.multi1, parallel.id=0,
                                end_buffer=NULL, calculate_p=TRUE,
                                return_data=FALSE) {
  if (is.null(end_buffer)) end_buffer <- get("end_buffer", envir=.GlobalEnv)
  parallel.id <- paste0("misspecify_SEIR", parallel.id)
  out.df <- gen_SEIR(trans_prob.base1, trans_prob.base2, eff.multi1, inf_mean,
                     delta=delta, N=N, N1=N1, pop.size=pop.size,
                     seed1=seed1, seed2=seed2, T0=T0, T1=T1,
                     burnin=burnin, end_buffer=end_buffer)
  data.in <- process_data(out.df, inf_mean, agg, dgp="SIR",
                          discard_start=burnin, discard_end=end_buffer) # misspecifying as SIR
  ################################################################################################################################
  Rt_exposure.out <- run_Rt(data.in, out.df, dgp="SIR", inf_mean,
                            parallel.id=parallel.id, calculate_p=calculate_p)
  beta.out <- run_beta(data.in, out.df, dgp="SIR", inf_mean,
                       parallel.id=parallel.id, calculate_p=calculate_p)
  Y.untrt.true <- run_true(out.df, trans_prob.base1, dgp="SEIR") # true under SEIR
  ############################################################################################################################
  # summarize outputs
  out <- rbind(Rt_exposure.out, beta.out) %>%
    mutate(dgp_true="SEIR", dgp_spe="SIR", N=N, N1=N1, trans_prob.base1=trans_prob.base1, trans_prob.base2=trans_prob.base2, 
           pop.size=pop.size, seed=seed1, eff.multi=eff.multi1, burnin=burnin, T0=T0, T1=T1, 
           S_frac.mean=mean(data.in$S_frac[data.in$trt.time]), S_frac.min=min(data.in$S_frac),
           Y.trt=mean(data.in$inc[data.in$trt_post]), Y.untrt.true=Y.untrt.true)
  if (return_data) return(list(result=out, bootstrap_data=data.in))
  out
}
