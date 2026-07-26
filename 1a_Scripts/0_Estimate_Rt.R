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

# Estimate cohort-based Rt according to Wallinga-Teunis.
compute_Rt_wt <- function(inf_mean, ID, inc, agg = 7, incubation_mean = 0) {
  # Under the SEIR DGP, transmission generations include both the latent and
  # infectious intervals.
  serial_interval_mean <- inf_mean + incubation_mean
  n_periods <- (T0 + T1 + burnin * 3) / agg

  if (length(n_periods) != 1L || !is.finite(n_periods) || n_periods < 3 ||
      abs(n_periods - round(n_periods)) > sqrt(.Machine$double.eps)) {
    stop("The Wallinga-Teunis analysis period must contain an integer number of at least three aggregation periods.")
  }

  n_periods <- as.integer(round(n_periods))
  t_start <- 2:(n_periods - 1L)
  t_end <- 3:n_periods
  units <- unique(ID)
  keep_wt <- vector("list", length(units))

  for (k in seq_along(units)) {
    i <- units[k]
    vec <- as.numeric(inc[ID == i])
    missing_result <- data.table(unit = i, week = t_end, R_wt = NA_real_)

    # Rt is undefined for a unit with no usable incidence.  This is a
    # unit-level estimator limitation, not extinction of the entire simulated
    # dataset; retain the replicate and return NA only for this unit's W-T Rt.
    if (length(vec) < n_periods || any(!is.finite(vec)) || any(vec < 0) ||
        sum(vec, na.rm = TRUE) <= 0) {
      keep_wt[[k]] <- missing_result
      next
    }

    sparse_or_terminally_extinct <-
      sum(vec > 0, na.rm = TRUE) < 2L ||
      !any(is.finite(tail(vec, 2L)) & tail(vec, 2L) > 0)

    temp_wt <- tryCatch(
      EpiEstim::wallinga_teunis(
        vec,
        method = "parametric_si",
        config = list(
          t_start = t_start,
          t_end = t_end,
          method = "parametric_si",
          mean_si = serial_interval_mean / agg,
          std_si = serial_interval_mean / agg,
          n_sim = 0
        )
      )$R,
      error = function(e) {
        # EpiEstim can emit this non-specific base-R error when a sparse or
        # terminally extinct unit does not support W-T estimation.  Handle only
        # that understood unit-level case; every other error remains fatal.
        if (sparse_or_terminally_extinct &&
            identical(conditionMessage(e), "missing value where TRUE/FALSE needed")) {
          return(NULL)
        }
        stop(e)
      }
    )

    if (is.null(temp_wt)) {
      keep_wt[[k]] <- missing_result
    } else {
      keep_wt[[k]] <- temp_wt %>%
        mutate(week = t_end, unit = i, type = "processed", R_wt = `Mean(R)`) %>%
        dplyr::select(unit, week, R_wt)
    }
  }

  rbindlist(keep_wt, use.names = TRUE, fill = TRUE)
}

# Estimate the true value for cohort-based Rt
compute_Rt_cohort <- function(incidence, ID, inf_mean, Ttot) {
  w <- (1 - 1/inf_mean)^(1:Ttot)
  w <- w / sum(w)
  Rt_cohort <- rep(NA_real_, length(incidence))

  for (i in unique(ID)) {
    idx <- which(ID == i)
    n_i <- length(idx)
    if (n_i < 2L) next

    for (t in seq_len(n_i - 1L)) {
      current_case <- incidence[idx[t]]
      future_cases <- incidence[idx[(t + 1L):n_i]]
      weights <- w[seq_along(future_cases)]

      if (length(current_case) != 1L || !is.finite(current_case) || current_case <= 0) {
        Rt_cohort[idx[t]] <- NA_real_
      } else {
        Rt_cohort[idx[t]] <- sum(future_cases * weights, na.rm = TRUE) / current_case
      }
    }
  }

  Rt_cohort
}

# simulation error handler
sim_error <- function(class, message) {
  stop(structure(
    list(message = message, call = NULL),
    class = c(class, "error", "condition")
  ))
}

# 1. Simulate data according to SIR or SEIR
# 2. Aggregate data to the desired level (default: weekly)
# 3. Estimate R_t and beta_t
gen_SIR <- function(trans_prob.base1, trans_prob.base2, eff.multi1, inf_mean, calculate_prevalence = T, equal_pop = T) {
  if (equal_pop) {
    pop.size1 <- pop.size2 <- pop.size
  }
  
  out.sim <- lapply(1:N, function(ind) { 
    if (ind %in% (1:N1)) { # generating treated units
      run_SIR_varying(pop.size=pop.size1, seeds=seed1, time_steps=(T0+T1+burnin*3), inf_mean=inf_mean,
                      trans_prob = c(rep(trans_prob.base1, (T0+burnin)), rep(trans_prob.base1*eff.multi1, (T1+burnin*2))))
    } else { # generating comparison units
      run_SIR_varying(pop.size=pop.size2, seeds=seed2, time_steps=(T0+T1+burnin*3), inf_mean=inf_mean,
                      trans_prob = c(rep(trans_prob.base2, (T0+burnin)), rep(trans_prob.base2, (T1+burnin*2))))
    }})
  out.df <- rbindlist(out.sim) %>% mutate(unit=rep(1:N, each=(T0+T1+burnin*3)))
  if (calculate_prevalence) {
    out.df$prevalence <- compute_prevalence(inf_mean=inf_mean, ID=out.df$unit, inc=out.df$inc, time=out.df$t, Ttot=T0+T1+burnin)
    out.df$R_cohort <- compute_Rt_cohort(out.df$inc, out.df$unit, inf_mean, (T0+T1+burnin*3)) 
  }
  
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
process_data <- function(out.df, inf_mean, agg=7, dgp="SIR", incubation_mean=0,
                         discard_start=burnin, discard_end=burnin*2) {
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
  
  # Classify only replicate-level extinction.  A single unit with no cases
  # must not discard the entire simulation; unit-level W-T non-estimability is
  # handled inside compute_Rt_wt() by returning NA for that unit.
  analysis.window <- df.agg %>%
    ungroup() %>%
    mutate(analysis_week = week - discard_start / agg)

  last.analysis.week <- max(analysis.window$analysis_week, na.rm = TRUE) -
    discard_end / agg

  analysis.window <- analysis.window %>%
    filter(analysis_week > 0, analysis_week <= last.analysis.week)

  if (nrow(analysis.window) == 0L) {
    sim_error(
      "epidemic_extinct",
      "No retained analysis weeks were available after burn-in removal."
    )
  }

  activity.by.group <- analysis.window %>%
    mutate(simulation_group = ifelse(trt.unit, "treated", "comparison")) %>%
    group_by(simulation_group) %>%
    summarise(
      active = any(is.finite(inc) & inc > 0),
      .groups = "drop"
    )

  inactive.groups <- activity.by.group$simulation_group[!activity.by.group$active]
  if (length(inactive.groups) > 0L) {
    sim_error(
      "epidemic_extinct",
      paste0(
        "No usable incidence was generated in the retained analysis period for group(s): ",
        paste(inactive.groups, collapse = ", ")
      )
    )
  }

  # NECESSARY CHANGE [E16-CALL]: receive the incubation mean explicitly.
  # R uses lexical (not caller/dynamic) scoping, so looking for `delta` inside
  # this helper can silently pick up an unrelated global value instead of the
  # SEIR parameter supplied to the simulation wrapper.
  wt.incubation <- if (dgp == "SEIR") incubation_mean else 0
  if (dgp == "SEIR" &&
      (length(wt.incubation) != 1 || is.na(wt.incubation) ||
       !is.finite(wt.incubation) || wt.incubation <= 0)) {
    stop("SEIR Wallinga-Teunis estimation requires a positive incubation_mean.")
  }

  wt.df <- compute_Rt_wt(
    inf_mean = inf_mean,
    ID = df.agg$unit,
    inc = df.agg$inc,
    agg = agg,
    incubation_mean = wt.incubation
  )
  data.clean <- df.agg %>% left_join(wt.df, c("unit" = "unit", "week" = "week")) %>%
    mutate(R_wt = R_wt / inf_mean, beta_wt = R_wt / S_frac)

  # chop off burnin in the beginning and 2*burnin periods in the end
  df <- data.clean %>% data.frame() %>%
    mutate(unit = relevel(factor(unit), ref = max(data.clean$unit)),
           week = week - discard_start/agg) %>%
    filter(week > 0,
           week <= max(week) - discard_end/agg)
  
  return(df)
}
