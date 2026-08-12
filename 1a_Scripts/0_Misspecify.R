source("./1a_Scripts/0_SIR.R")
source("./1a_Scripts/0_SEIR.R")

# mean_spe, var_spe are the multiplicative factor on the misspecification
# they range from (0.8, 1.2)
sim_misspecify_GI <- function(mean_true, var_true, eff.multi1, mean_spe, var_spe, parallel.id=0) {
  parallel.id <- paste0("misspecify_GI_", parallel.id)
  trt.IDs <- 1:N1
  out.sim <- lapply(1:N, function(ind) { 
    if (ind %in% trt.IDs) {
      run_SIR_varying(pop.size=pop.size, seeds=seed1, time_steps=(T0+T1+burnin*3), inf_mean=mean_true, inf_var=var_true,
                      trans_prob = c(rep(trans_prob.base1, (T0+burnin)), rep(trans_prob.base1*eff.multi1, (T1+burnin*2))))
    } else {
      run_SIR_varying(pop.size=pop.size, seeds=seed2, time_steps=(T0+T1+burnin*3), inf_mean=mean_true, inf_var=var_true,
                      trans_prob = rep(trans_prob.base2, (T0+T1+burnin*3)))
    }})
  out.df <- rbindlist(out.sim) %>%
    mutate(unit=rep(1:N, each=(T0+T1+burnin*3)),
           initial_seed=ifelse(unit %in% trt.IDs, seed1, seed2))
  
  # Misspecification: (mean_spe = var_spe = 1) corresponds to the base case of NO MISPECIFICATION
  inf_mean_spe <- mean_true * mean_spe
  inf_var_spe  <- var_true * var_spe
  inf_std_spe <- sqrt(inf_var_spe)
  
  # Convert (inf_mean_spe, inf_var_spe) into a Gamma(shape, scale):
  shape_spe <- inf_mean_spe^2 / inf_var_spe
  scale_spe <- inf_var_spe / inf_mean_spe
  
  # Estimate prevalence using Gamma. Prepend the initial seed at t=0 before
  # evaluating the incidence-history convolution, consistent with all other
  # simulation-based prevalence/Rt/beta calculations.
  out.df$prevalence_gamma <- NA_real_
  for (i in 1:N) {
    idx <- which(out.df$unit == i)
    inc_i <- out.df$inc[idx]
    seed_i <- unique(out.df$initial_seed[idx])
    if (length(seed_i) != 1L || !is.finite(seed_i) || seed_i < 0) {
      stop("Each simulated unit must have one finite, nonnegative initial seed.")
    }
    inc_history <- c(seed_i, inc_i)
    last_t <- min(T0 + T1 + burnin, length(inc_i))
    for (t in seq_len(last_t)) {
      lags <- 0:t
      surv <- 1 - pgamma(lags, shape = shape_spe, scale = scale_spe)
      out.df$prevalence_gamma[idx[t]] <- sum(surv * inc_history[(t + 1L) - lags])
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
              R_est = sum(inc) / sum(prevalence_lag),
              initial_seed = first(initial_seed)) %>%
    group_by(unit) %>%
    mutate(beta_true = R_true / S_frac,
           beta_est = R_est / S_frac,
           trt.unit = (unit %in% trt.IDs),
           trt.time = (week > (T0+burnin)/agg),
           trt_post = (trt.unit & trt.time),
           time_to_trt = ifelse(!trt_post, 0, week - (T0+burnin)/agg))
  
  df.agg <- df.agg %>%
    group_by(unit) %>%
    mutate(inc_for_Rt = inc + ifelse(week == min(week), initial_seed, 0)) %>%
    ungroup()

  keep_wt <- list()
  for(i in unique(df.agg$unit)){
    vec <- df.agg$inc_for_Rt[df.agg$unit==i]
    temp_wt <- wallinga_teunis(vec, method="parametric_si",
                               config=list(
                                 t_start=2:((T0+T1+burnin*3)/agg-1),
                                 t_end=3:((T0+T1+burnin*3)/agg),
                                 method="parametric_si", 
                                 mean_si=inf_mean_spe/agg, std_si=inf_std_spe/agg,
                                 n_sim=0))$R %>%
      mutate(week=t_end, unit=i, type="processed", R_wt=`Mean(R)`)
    keep_wt <- rbindlist(list(keep_wt, temp_wt))
  }
  keep <- keep_wt %>% dplyr::select(unit, week, R_wt)
  df <- df.agg %>% left_join(keep, c("unit"="unit", "week"="week")) %>%
    mutate(R_wt = R_wt / inf_mean_spe, beta_wt = R_wt / S_frac)
  
  # chop off burnin in the beginning and 2*burnin periods in the end
  data.in <- df %>% data.frame() %>%
    mutate(unit = relevel(factor(unit), ref = N),
           week = week - burnin/agg) %>%
    filter(week > 0, week <= max(week) - burnin*2/agg)
  ############################################################################################################################
  # Fit regression models
  Rt_wt.out <- run_Rt(data.in, out.df, type="wt", dgp="SIR", inf_mean=inf_mean_spe, inf_var=inf_var_spe, parallel.id=parallel.id)
  Rt_est.out <- run_Rt(data.in, out.df, type="est", dgp="SIR", inf_mean=inf_mean_spe, inf_var=inf_var_spe, parallel.id=parallel.id)
  beta.out <- run_beta(data.in, out.df, dgp="SIR", inf_mean=inf_mean_spe, inf_var=inf_var_spe, parallel.id=parallel.id)
  Y.untrt.true <- run_true(out.df, trans_prob.base1, dgp="SIR")
  ################################################################################################################
  # summarize outputs
  out <- rbind(Rt_wt.out, Rt_est.out, beta.out) %>% 
    mutate(eff.multi=eff.multi1, seed=seed1, mean_spe=mean_spe, var_spe=var_spe, 
           inf_mean=mean_true, S_frac.mean=mean(data.in$S_frac[data.in$trt.time]), S_frac.min=min(data.in$S_frac),
           Y.trt=mean(data.in$inc[data.in$trt_post]), Y.untrt.true=Y.untrt.true)
  out
}

sim_misspecify_SEIR <- function(eff.multi1, parallel.id=0) {
  parallel.id <- paste0("misspecify_SEIR", parallel.id)
  out.df <- gen_SEIR(trans_prob.base1, trans_prob.base2, eff.multi1, inf_mean)
  data.in <- process_data(out.df, inf_mean, agg, dgp="SIR") # misspecifying as SIR
  ################################################################################################################################
  Rt_wt.out <- run_Rt(data.in, out.df, type="wt", dgp="SIR", inf_mean, parallel.id=parallel.id)
  Rt_est.out <- run_Rt(data.in, out.df, type="est", dgp="SIR", inf_mean, parallel.id=parallel.id)
  beta.out <- run_beta(data.in, out.df, dgp="SIR", inf_mean, parallel.id=parallel.id)
  Y.untrt.true <- run_true(out.df, trans_prob.base1, dgp="SEIR") # true under SEIR
  ############################################################################################################################
  # summarize outputs
  out <- rbind(Rt_wt.out, Rt_est.out, beta.out) %>% 
    mutate(dgp_true="SEIR", dgp_spe="SIR", N=N, N1=N1, trans_prob.base1=trans_prob.base1, trans_prob.base2=trans_prob.base2, 
           pop.size=pop.size, seed=seed1, eff.multi=eff.multi1, burnin=burnin, T0=T0, T1=T1, 
           S_frac.mean=mean(data.in$S_frac[data.in$trt.time]), S_frac.min=min(data.in$S_frac),
           Y.trt=mean(data.in$inc[data.in$trt_post]), Y.untrt.true=Y.untrt.true)
  out
}
