# Calculate the expected number of active infection on day t by looking at past incidences
compute_prevalence <- function(inf_mean, ID, inc, time, Ttot) {
  prevalence <- rep(NA, length(inc))
  for (i in unique(ID)) {
    for (t in 2:Ttot) {
      prevalence[ID==i & time==t] <- sum(sapply(0:(t-1), function(j) {
        # sum up (inc[t] + w*inc[t-1] + w^2*inc[t-2] + ...)
        (1 - 1/inf_mean)^j * inc[ID==i & time==(t-j)]
      }), na.rm = T)
    }
  }
  return(prevalence)
}

# Calculate the expected number of individuals currently in the Exposed state on day t
compute_infected <- function(delta, ID, inc, time, Ttot) {
  infected_est <- rep(NA, length(inc))
  for (i in unique(ID)) {
    for (t in 2:Ttot) {
      inc_today <- inc[ID==i & time==t]
      inc_tomorrow <- inc[ID==i & time==(t+1)]
      inc_today <- ifelse(is.na(inc_today), 0, inc_today)
      inc_tomorrow <- ifelse(is.na(inc_tomorrow), 0, inc_tomorrow)
      
      infected_est[ID==i & time==t] <- delta * inc_tomorrow - delta * (1 - 1/delta) * inc_today
    }
  }
  return(infected_est)
}

# Estimate cohort-based Rt according to Wallinga Teunis
compute_Rt_wt <- function(inf_mean, ID, inc, agg=7) {
  keep_wt <- list()
  for(i in unique(ID)){
    vec <- inc[ID==i]
    temp_wt <- wallinga_teunis(vec, method="parametric_si",
                               config=list(
                                 t_start=2:((T0+T1+burnin*3)/agg-1),
                                 t_end=3:((T0+T1+burnin*3)/agg),
                                 method="parametric_si", 
                                 mean_si=inf_mean/agg, std_si=inf_mean/agg,
                                 n_sim=0))$R %>%
      mutate(week=t_end, unit=i, type="processed", R_wt=`Mean(R)`)
    keep_wt <- rbindlist(list(keep_wt, temp_wt))
  }
  return(keep_wt %>% dplyr::select(unit, week, R_wt))
}

# Estimate the true value for cohort-based Rt
compute_Rt_cohort <- function(incidence, ID, inf_mean, Ttot) {
  w <- (1 - 1/inf_mean)^(1:Ttot)
  w <- w / sum(w)  # normalize
  Rt_cohort <- rep(NA, length(incidence))
  
  for (i in unique(ID)) {
    for (t in 1:(Ttot - 1)) {
      future_cases <- incidence[ID==i][(t + 1):Ttot]
      current_case <- incidence[ID==i][t]
      weights <- w[1:length(future_cases)]
      
      if (current_case == 0) {
        Rt_cohort[ID==i][t] <- NA
      } else {
        Rt_cohort[ID==i][t] <- sum(future_cases * weights) / current_case
      }
    }
  }

  return(Rt_cohort)
}

# 1. Simulate data according to SIR or SEIR
# 2. Aggregate data to the desired level (default: weekly)
# 3. Estimate R_t and beta_t
gen_SIR <- function(trans_prob.base1, trans_prob.base2, eff.multi1, inf_mean) {
  out.sim <- lapply(1:N, function(ind) { 
    if (ind %in% (1:N1)) { # generating treated units
      run_SIR_varying(pop.size=pop.size, seeds=seed1, time_steps=(T0+T1+burnin*3), inf_mean=inf_mean,
                      trans_prob = c(rep(trans_prob.base1, (T0+burnin)), rep(trans_prob.base1*eff.multi1, (T1+burnin*2))))
    } else { # generating comparison units
      run_SIR_varying(pop.size=pop.size, seeds=seed2, time_steps=(T0+T1+burnin*3), inf_mean=inf_mean,
                      trans_prob = c(rep(trans_prob.base2, (T0+burnin)), rep(trans_prob.base2, (T1+burnin*2))))
    }})
  out.df <- rbindlist(out.sim) %>% mutate(unit=rep(1:N, each=(T0+T1+burnin*3)))
  out.df$prevalence <- compute_prevalence(inf_mean=inf_mean, ID=out.df$unit, inc=out.df$inc, time=out.df$t, Ttot=T0+T1+burnin)
  out.df$R_cohort <- compute_Rt_cohort(out.df$inc, out.df$unit, inf_mean, (T0+T1+burnin*3))
  
  return(out.df)
}

gen_SEIR <- function(trans_prob.base1, trans_prob.base2, eff.multi1, inf_mean) {
  out.sim <- lapply(1:N, function(ind) { 
    if (ind %in% (1:N1)) {
      run_SEIR_varying(pop.size=pop.size, I0=seed1, time_steps=(T0+T1+burnin*3), inf_mean=inf_mean, delta=delta,
                       trans_prob = c(rep(trans_prob.base1, (T0+burnin)), rep(trans_prob.base1*eff.multi1, (T1+burnin*2))))
    } else {
      run_SEIR_varying(pop.size=pop.size, I0=seed2, time_steps=(T0+T1+burnin*3), inf_mean=inf_mean, delta=delta,
                       trans_prob = c(rep(trans_prob.base2, (T0+burnin)), rep(trans_prob.base2, (T1+burnin*2))))
    }})
  out.df <- rbindlist(out.sim) %>% mutate(unit=rep(1:N, each=(T0+T1+burnin*3)))
  out.df$prevalence <- compute_prevalence(inf_mean=inf_mean, ID=out.df$unit, inc=out.df$inc, time=out.df$t, Ttot=T0+T1+burnin)
  out.df$infected_est <- compute_infected(delta=delta, ID=out.df$unit, inc=out.df$inc, time=out.df$t, Ttot=T0+T1+burnin)
  out.df$R_cohort <- compute_Rt_cohort(out.df$inc, out.df$unit, inf_mean, (T0+T1+burnin*3))
  
  return(out.df)
}

# dgp is either "SIR" or "SEIR"
process_data <- function(out.df, inf_mean, agg=7, dgp="SIR") {
  df.agg <- out.df %>% 
    group_by(unit) %>%
    arrange(t) %>%
    mutate(week = ceiling(t/agg),
           S_lag = ifelse(t==1, (pop.size-seed1), lag(S, 1)),
           prevalence_lag = lag(prevalence, 1),
           inc_lag = lag(inc, 1),
           I_lag = ifelse(t==1, seed1, lag(I, 1))) %>%
    filter(t >= 2) %>%
    group_by(unit, week) %>%
    summarise(inc = sum(inc),
              growth = sum(inc) / sum(inc_lag),
              S_frac = sum(S_lag) / (pop.size*agg),
              R_true = sum(inc) / sum(I_lag),
              # if SIR, the numerator is the total number of observed incidence
              # if SEIR, the numerator is the estimated total number of infected
              R_est = ifelse(dgp=="SIR", sum(inc) / sum(prevalence_lag),
                             sum(infected_est) / sum(prevalence_lag)),
              R_cohort = mean(R_cohort, na.rm = TRUE) / inf_mean) %>%
    group_by(unit) %>%
    mutate(beta_est = R_est / S_frac,
           trt.unit = (unit %in% (1:N1)),
           trt.time = (week > (T0+burnin)/agg),
           trt_post = (trt.unit & trt.time))
  
  wt.df <- compute_Rt_wt(inf_mean=inf_mean, ID=df.agg$unit, inc=df.agg$inc)
  data.clean <- df.agg %>% left_join(wt.df, c("unit"="unit", "week"="week")) %>%
    mutate(R_wt = R_wt / inf_mean, beta_wt = R_wt / S_frac)
  
  # chop off burnin in the beginning and 2*burnin periods in the end
  df <- data.clean %>% data.frame() %>%
    mutate(unit = relevel(factor(unit), ref = max(data.clean$unit)),
           week = week - burnin/agg) %>%
    filter(week > 0, 
           week <= max(week) - burnin*2/agg)
  
  return(df)
}
