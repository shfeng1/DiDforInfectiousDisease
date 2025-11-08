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
  out.df <- rbindlist(out.sim) %>% mutate(unit=rep(1:N, each=(T0+T1+burnin*3)))
  
  # Misspecification: (mean_spe = var_spe = 1) corresponds to the base case of NO MISPECIFICATION
  inf_mean_spe <- mean_true * mean_spe
  inf_var_spe  <- var_true * var_spe
  inf_std_spe <- sqrt(inf_var_spe)
  
  # Convert (inf_mean_spe, inf_var_spe) into a Gamma(shape, scale):
  shape_spe <- inf_mean_spe^2 / inf_var_spe
  scale_spe <- inf_var_spe / inf_mean_spe
  
  # Estimate prevalence using Gamma
  out.df$prevalence_gamma <- NA
  for (i in 1:N) {
    for (t in 2:(T0 + T1 + burnin)) {
      out.df$prevalence_gamma[out.df$unit==i & out.df$t==t] <- 
        sum(sapply(0:(t-1), function(j) {
          surv_j <- 1 - pgamma(j, shape = shape_spe, scale = scale_spe)
          surv_j * out.df$inc[out.df$unit==i & out.df$t==t-j]
        }))
    }
  }
  
  df.agg <- out.df %>% 
    group_by(unit) %>%
    arrange(t) %>%
    mutate(week = ceiling(t/agg),
           S_lag = ifelse(t==1, (pop.size-seed1), lag(S, 1)),
           prevalence_lag = lag(prevalence_gamma, 1),
           inc_lag = lag(inc, 1),
           I_lag = ifelse(t==1, seed1, lag(I, 1))) %>%
    filter(t >= 3) %>%
    group_by(unit, week) %>%
    summarise(inc = sum(inc),
              growth = sum(inc) / sum(inc_lag),
              S_frac = sum(S_lag) / (pop.size*agg),
              R_true = sum(inc) / sum(I_lag),
              R_est = sum(inc) / sum(prevalence_lag)) %>%
    group_by(unit) %>%
    mutate(beta_true = R_true / S_frac,
           beta_est = R_est / S_frac,
           trt.unit = (unit %in% trt.IDs),
           trt.time = (week > (T0+burnin)/agg),
           trt_post = (trt.unit & trt.time),
           time_to_trt = ifelse(!trt_post, 0, week - (T0+burnin)/agg))
  
  keep_wt <- list()
  for(i in unique(df.agg$unit)){
    vec <- df.agg$inc[df.agg$unit==i]
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
    mutate(R_wt = R_wt / (inf_mean * mean_spe), beta_wt = R_wt / S_frac)
  
  # chop off burnin in the beginning and 2*burnin periods in the end
  data.in <- df %>% data.frame() %>%
    mutate(unit = relevel(factor(unit), ref = N),
           week = week - burnin/agg) %>%
    filter(week > 0, week <= max(week) - burnin*2/agg)
  ############################################################################################################################
  # Fit regression models
  Rt_wt.out <- run_Rt(data.in, out.df, type="wt", dgp="SIR", inf_mean_spe, inf_var_spe, parallel.id=parallel.id)
  Rt_est.out <- run_Rt(data.in, out.df, type="est", dgp="SIR", inf_mean_spe, inf_var_spe, parallel.id=parallel.id)
  beta.out <- run_beta(data.in, out.df, dgp="SIR", inf_mean_spe, inf_var_spe, parallel.id=parallel.id)
  Y.untrt.true <- run_true(out.df, trans_prob.base1, dgp="SIR")
  ################################################################################################################
  # summarize outputs
  out <- rbind(Rt_wt.out, Rt_est.out, beta.out) %>% 
    mutate(eff.multi=eff.multi1, seed=seed1, mean_spe=mean_spe, var_spe=var_spe, 
           inf_mean=inf_mean, S_frac.mean=mean(data.in$S_frac[data.in$trt.time]), S_frac.min=min(data.in$S_frac),
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
