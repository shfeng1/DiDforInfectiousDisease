# Calculate the expected number of active infection on day t by looking at past incidences
compute_prevalence <- function(inf_mean, ID, inc, time, Ttot) {
  prevalence <- rep(NA, length(inc))
  ID.key <- if (inherits(ID, "haven_labelled")) {
    as.character(unclass(ID))
  } else {
    as.character(ID)
  }
  for (i in unique(ID.key)) {
    for (t in 2:Ttot) {
      prevalence[ID.key==i & time==t] <- sum(sapply(0:(t-1), function(j) {
        # sum up (inc[t] + w*inc[t-1] + w^2*inc[t-2] + ...)
        (1 - 1/inf_mean)^j * inc[ID.key==i & time==(t-j)]
      }), na.rm = T)
    }
  }
  return(prevalence)
}

# Calculate the expected number of individuals currently in the Exposed state on day t
compute_infected <- function(delta, ID, inc, time, Ttot) {
  infected_est <- rep(NA, length(inc))
  ID.key <- if (inherits(ID, "haven_labelled")) {
    as.character(unclass(ID))
  } else {
    as.character(ID)
  }
  for (i in unique(ID.key)) {
    for (t in 2:Ttot) {
      inc_today <- inc[ID.key==i & time==t]
      inc_tomorrow <- inc[ID.key==i & time==(t+1)]
      inc_today <- ifelse(is.na(inc_today), 0, inc_today)
      inc_tomorrow <- ifelse(is.na(inc_tomorrow), 0, inc_tomorrow)
      
      infected_est[ID.key==i & time==t] <- delta * inc_tomorrow - delta * (1 - 1/delta) * inc_today
    }
  }
  return(infected_est)
}

compute_simulation_prevalence <- function(inf_mean, ID, inc, time, Ttot) {
  prevalence <- rep(NA, length(inc))
  ID.key <- as.character(ID)
  decay <- 1-1/inf_mean
  for (i in unique(ID.key)) {
    idx <- which(ID.key==i)
    idx <- idx[order(time[idx])]
    unit.time <- time[idx]
    first <- match(1L, unit.time)
    first.inc <- if (is.na(first)) NA else inc[idx[first]]
    active <- if (is.finite(first.inc)) first.inc else 0
    if (Ttot>=2L) {
      for (t in 2:as.integer(Ttot)) {
        pos <- match(t, unit.time)
        current <- if (is.na(pos)) NA else inc[idx[pos]]
        active <- (if (is.finite(current)) current else 0)+decay*active
        if (!is.na(pos)) prevalence[idx[pos]] <- active
      }
    }
  }
  prevalence
}

compute_simulation_infected <- function(delta, ID, inc, time, Ttot) {
  infected_est <- rep(NA, length(inc))
  ID.key <- as.character(ID)
  for (i in unique(ID.key)) {
    idx <- which(ID.key==i)
    idx <- idx[order(time[idx])]
    unit.time <- time[idx]
    for (t in 2:as.integer(Ttot)) {
      pos <- match(t, unit.time)
      if (is.na(pos)) next
      tomorrow <- match(t+1L, unit.time)
      inc.today <- inc[idx[pos]]
      inc.tomorrow <- if (is.na(tomorrow)) NA else inc[idx[tomorrow]]
      if (!is.finite(inc.today)) inc.today <- 0
      if (!is.finite(inc.tomorrow)) inc.tomorrow <- 0
      infected_est[idx[pos]] <- delta*inc.tomorrow-
        delta*(1-1/delta)*inc.today
    }
  }
  infected_est
}

# Build the full incidence history used by prevalence, Rt, and beta estimators.
# The simulation functions return only observations after the initial state, so
# prepend each unit's initial infectious seed as incidence at t=0.
build_seeded_incidence_history <- function(out.df) {
  required <- c("unit", "t", "inc", "I", "initial_seed")
  missing <- setdiff(required, names(out.df))
  if (length(missing)) {
    stop("Seeded incidence history requires column(s): ", paste(missing, collapse=", "), ".")
  }

  seed.table <- out.df %>%
    transmute(.unit_key=as.character(unit), initial_seed=as.numeric(initial_seed)) %>%
    group_by(.unit_key) %>%
    summarise(initial_seed={
      values <- unique(initial_seed)
      if (length(values)==1L) values else NA
    }, .groups="drop")

  if (anyNA(seed.table$initial_seed) || any(!is.finite(seed.table$initial_seed)) ||
      any(seed.table$initial_seed < 0)) {
    stop("Each simulated unit must have one finite, nonnegative initial_seed value.")
  }

  observed <- out.df %>%
    transmute(.unit_key=as.character(unit), original_t=as.integer(t),
              seeded_t=as.integer(t)+1L, inc=as.numeric(inc), I=as.numeric(I),
              .seed_row=FALSE)
  seed.rows <- seed.table %>%
    transmute(.unit_key, original_t=0L, seeded_t=1L,
              inc=initial_seed, I=initial_seed, .seed_row=TRUE)

  bind_rows(seed.rows, observed) %>% arrange(.unit_key, seeded_t)
}

map_seeded_estimate_to_output <- function(out.df, seeded.history, value.column) {
  lookup <- seeded.history %>%
    filter(!.seed_row) %>%
    transmute(.key=paste(.unit_key, original_t, sep="::"),
              .value=.data[[value.column]])
  out.key <- paste(as.character(out.df$unit), as.integer(out.df$t), sep="::")
  lookup$.value[match(out.key, lookup$.key)]
}

compute_prevalence_seeded <- function(out.df, inf_mean, Ttot) {
  history <- build_seeded_incidence_history(out.df)
  history$prevalence <- compute_simulation_prevalence(
    inf_mean=inf_mean, ID=history$.unit_key, inc=history$inc,
    time=history$seeded_t, Ttot=as.integer(Ttot)+1L)
  map_seeded_estimate_to_output(out.df, history, "prevalence")
}

compute_infected_seeded <- function(out.df, delta, Ttot) {
  history <- build_seeded_incidence_history(out.df)
  history$infected_est <- compute_simulation_infected(
    delta=delta, ID=history$.unit_key, inc=history$inc,
    time=history$seeded_t, Ttot=as.integer(Ttot))
  map_seeded_estimate_to_output(out.df, history, "infected_est")
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
gen_SIR <- function(trans_prob.base1, trans_prob.base2, eff.multi1, inf_mean,
                    calculate_prevalence=TRUE, equal_pop=TRUE,
                    N=NULL, N1=NULL, pop.size=NULL, seed1=NULL, seed2=NULL,
                    T0=NULL, T1=NULL, burnin=NULL, end_buffer=NULL) {
  if (is.null(N)) N <- get("N", envir=.GlobalEnv)
  if (is.null(N1)) N1 <- get("N1", envir=.GlobalEnv)
  if (is.null(pop.size)) pop.size <- get("pop.size", envir=.GlobalEnv)
  if (is.null(seed1)) seed1 <- get("seed1", envir=.GlobalEnv)
  if (is.null(seed2)) seed2 <- get("seed2", envir=.GlobalEnv)
  if (is.null(T0)) T0 <- get("T0", envir=.GlobalEnv)
  if (is.null(T1)) T1 <- get("T1", envir=.GlobalEnv)
  if (is.null(burnin)) burnin <- get("burnin", envir=.GlobalEnv)
  if (is.null(end_buffer)) end_buffer <- get("end_buffer", envir=.GlobalEnv)
  total.days <- T0+T1+burnin+end_buffer
  if (equal_pop) pop.size1 <- pop.size2 <- pop.size
  out.sim <- lapply(1:N, function(ind) {
    if (ind %in% (1:N1)) {
      run_SIR_varying(pop.size=pop.size1, seeds=seed1,
                      time_steps=total.days, inf_mean=inf_mean,
                      trans_prob=c(rep(trans_prob.base1, T0+burnin),
                                   rep(trans_prob.base1*eff.multi1, T1+end_buffer)))
    } else {
      run_SIR_varying(pop.size=pop.size2, seeds=seed2,
                      time_steps=total.days, inf_mean=inf_mean,
                      trans_prob=c(rep(trans_prob.base2, T0+burnin),
                                   rep(trans_prob.base2, T1+end_buffer)))
    }
  })
  out.df <- rbindlist(out.sim) %>%
    mutate(unit=rep(1:N, each=total.days),
           initial_seed=ifelse(unit %in% (1:N1), seed1, seed2))
  if (calculate_prevalence) {
    out.df$prevalence <- compute_prevalence_seeded(out.df, inf_mean, total.days)
  }
  out.df
}

gen_SEIR <- function(trans_prob.base1, trans_prob.base2, eff.multi1, inf_mean,
                     delta=NULL, N=NULL, N1=NULL, pop.size=NULL,
                     seed1=NULL, seed2=NULL, T0=NULL, T1=NULL,
                     burnin=NULL, end_buffer=NULL) {
  if (is.null(delta)) delta <- get("delta", envir=.GlobalEnv)
  if (is.null(N)) N <- get("N", envir=.GlobalEnv)
  if (is.null(N1)) N1 <- get("N1", envir=.GlobalEnv)
  if (is.null(pop.size)) pop.size <- get("pop.size", envir=.GlobalEnv)
  if (is.null(seed1)) seed1 <- get("seed1", envir=.GlobalEnv)
  if (is.null(seed2)) seed2 <- get("seed2", envir=.GlobalEnv)
  if (is.null(T0)) T0 <- get("T0", envir=.GlobalEnv)
  if (is.null(T1)) T1 <- get("T1", envir=.GlobalEnv)
  if (is.null(burnin)) burnin <- get("burnin", envir=.GlobalEnv)
  if (is.null(end_buffer)) end_buffer <- get("end_buffer", envir=.GlobalEnv)
  total.days <- T0+T1+burnin+end_buffer
  out.sim <- lapply(1:N, function(ind) {
    if (ind %in% (1:N1)) {
      run_SEIR_varying(pop.size=pop.size, I0=seed1,
                       time_steps=total.days, inf_mean=inf_mean, delta=delta,
                       trans_prob=c(rep(trans_prob.base1, T0+burnin),
                                    rep(trans_prob.base1*eff.multi1, T1+end_buffer)))
    } else {
      run_SEIR_varying(pop.size=pop.size, I0=seed2,
                       time_steps=total.days, inf_mean=inf_mean, delta=delta,
                       trans_prob=c(rep(trans_prob.base2, T0+burnin),
                                    rep(trans_prob.base2, T1+end_buffer)))
    }
  })
  out.df <- rbindlist(out.sim) %>%
    mutate(unit=rep(1:N, each=total.days),
           initial_seed=ifelse(unit %in% (1:N1), seed1, seed2))
  out.df$prevalence <- compute_prevalence_seeded(out.df, inf_mean, total.days)
  out.df$infected_est <- compute_infected_seeded(out.df, delta, total.days)
  out.df
}

# dgp is either "SIR" or "SEIR"
process_data <- function(out.df, inf_mean, agg=7, dgp="SIR", incubation_mean=0,
                         discard_start=NULL, discard_end=NULL) {
  if (is.null(discard_start)) discard_start <- get("burnin", envir=.GlobalEnv)
  if (is.null(discard_end)) discard_end <- get("end_buffer", envir=.GlobalEnv)
  if (!("initial_seed" %in% names(out.df))) {
    stop("Simulated data passed to process_data() must contain initial_seed.")
  }
  df.agg <- out.df %>%
    group_by(unit) %>% arrange(t) %>%
    mutate(week=ceiling(t/agg),
           S_lag=ifelse(t==1, pop.size-initial_seed, lag(S,1)),
           prevalence_lag=lag(prevalence,1), inc_lag=lag(inc,1),
           I_lag=ifelse(t==1, initial_seed, lag(I,1))) %>%
    filter(t>=2) %>% group_by(unit, week) %>%
    summarise(inc=sum(inc), growth=sum(inc)/sum(inc_lag),
              S_frac=sum(S_lag)/(pop.size*agg),
              R_true=ifelse(dgp=="SIR", sum(inc)/sum(I_lag),
                            sum(infected)/sum(I_lag)),
              Rt_exposure=ifelse(dgp=="SIR", sum(inc)/sum(prevalence_lag),
                                 sum(infected_est)/sum(prevalence_lag)),
              .groups="drop_last") %>%
    mutate(beta_exposure=Rt_exposure/S_frac,
           trt.unit=unit %in% (1:N1),
           trt.time=week > (T0+burnin)/agg,
           trt_post=trt.unit & trt.time)

  analysis.window <- df.agg %>% ungroup() %>%
    mutate(analysis_week=week-discard_start/agg)
  last.analysis.week <- max(analysis.window$analysis_week, na.rm=TRUE)-discard_end/agg
  analysis.window <- analysis.window %>%
    filter(analysis_week>0, analysis_week<=last.analysis.week)
  if (nrow(analysis.window)==0L) {
    sim_error("epidemic_extinct", "No retained analysis weeks were available after burn-in removal.")
  }
  activity.by.group <- analysis.window %>%
    mutate(simulation_group=ifelse(trt.unit, "treated", "comparison")) %>%
    group_by(simulation_group) %>%
    summarise(active=any(is.finite(inc) & inc>0), .groups="drop")
  inactive.groups <- activity.by.group$simulation_group[!activity.by.group$active]
  if (length(inactive.groups)>0L) {
    sim_error("epidemic_extinct",
              paste0("No usable incidence was generated in the retained analysis period for group(s): ",
                     paste(inactive.groups, collapse=", ")))
  }

  df.agg %>% data.frame() %>%
    mutate(unit=relevel(factor(unit), ref=max(df.agg$unit)),
           week=week-discard_start/agg) %>%
    filter(week>0, week<=max(week)-discard_end/agg)
}
