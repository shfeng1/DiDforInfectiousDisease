# The helper functions below run each estimator in the following steps:
# 1. Fit the model.
# 2. Recover the untreated incidence trajectory.
# 3. Apply the manuscript bias correction to the fitted untreated trajectory.
# 4. Calculate AME as observed treated incidence minus fitted untreated incidence.
#
# For simulation-study calls (coef=NULL), the wild-bootstrap p-value is still
# calculated because the simulation summaries use it for power/type-I error.
# For Kansas AME conversion (coef supplied), no p-value is calculated or returned.

source("./1a_Scripts/0_RStata.R")
source("./1a_Scripts/0_Bias_Correction.R")

get_unit_population <- function(ind, unit_population=NULL) {
  if (is.null(unit_population)) pop.size else unname(unit_population[as.character(ind)])
}

aggregate_simulated_incidence <- function(sim.data, trt.IDs=NULL,
                                          incidence_aggregation="sum",
                                          unit_population=NULL,
                                          incidence_scale=NULL,
                                          week_offset=0) {
  if (!("unit" %in% names(sim.data))) {
    sim.data$unit <- rep(trt.IDs, each=nrow(sim.data)/length(trt.IDs))
  }
  sim.data <- sim.data %>%
    mutate(week=ceiling(t/agg)+week_offset) %>%
    group_by(unit, week) %>%
    summarise(inc=if (incidence_aggregation=="mean") mean(inc) else sum(inc), .groups="drop")
  if (!is.null(incidence_scale)) {
    sim.data$inc <- sim.data$inc / unname(unit_population[as.character(sim.data$unit)]) * incidence_scale
  }
  sim.data
}

adjust_fitted_untreated <- function(untreated, variance, approach="approach2") {
  if (approach=="approach1") untreated/(1+variance/2) else untreated/exp(variance/2)
}

prepare_ame_simulation_window <- function(data_deagg, ind, untreated_rate, dgp,
                                          simulate_from_trt=FALSE,
                                          treated_rate=NULL) {
  unit.data <- data_deagg[as.character(data_deagg$unit)==as.character(ind), , drop=FALSE]
  unit.data <- unit.data[order(unit.data$t), , drop=FALSE]
  if (simulate_from_trt) {
    path.data <- unit.data[unit.data$trt.time, , drop=FALSE]
    state <- unit.data[unit.data$t==min(path.data$t)-1, , drop=FALSE]
  } else {
    path.data <- unit.data
    state <- unit.data[1, , drop=FALSE]
  }
  list(
    unit=as.character(ind),
    time_steps=nrow(path.data),
    I0=as.numeric(state$I),
    E0=if (dgp=="SEIR") as.numeric(state$E) else NULL,
    R0=as.numeric(state$R),
    S0=if (simulate_from_trt) as.numeric(state$S) else NULL,
    treated_path=if (is.null(treated_rate)) NULL else as.numeric(path.data[[treated_rate]]),
    untreated_path=as.numeric(path.data[[untreated_rate]]),
    week_offset=if (simulate_from_trt) min(path.data$week)-1 else 0
  )
}

simulate_untreated_ame_trajectory <- function(spec, dgp, pop.ind, inf_mean,
                                               delta=NULL, inf_var=NULL) {
  out <- if (dgp=="SIR") {
    run_SIR_varying(pop.size=pop.ind, seeds=spec$I0, recovered=spec$R0,
                    S0=spec$S0, trans_prob=spec$untreated_path,
                    time_steps=spec$time_steps, inf_mean=inf_mean, inf_var=inf_var)
  } else {
    run_SEIR_varying(pop.size=pop.ind, I0=spec$I0, E0=spec$E0,
                     recovered=spec$R0, S0=spec$S0,
                     trans_prob=spec$untreated_path,
                     time_steps=spec$time_steps, inf_mean=inf_mean, delta=delta)
  }
  out$unit <- spec$unit
  out
}

# Paired treated/untreated mechanistic trajectories. This is used only when
# difference=TRUE. Both paths use the same U(0,1) draw at each time step.
simulate_paired_ame_trajectories <- function(spec, dgp, pop.ind, inf_mean,
                                             delta=NULL, inf_var=NULL) {
  Ttot <- spec$time_steps+1
  beta1 <- c(NA, spec$treated_path)
  beta0 <- c(NA, spec$untreated_path)
  S1 <- S0 <- I1 <- I0 <- R1 <- R0 <- numeric(Ttot)
  inc1 <- inc0 <- numeric(Ttot)
  if (is.null(spec$S0)) {
    start.S <- pop.ind - spec$I0 - spec$R0
    if (dgp == "SEIR") start.S <- start.S - spec$E0
  } else {
    start.S <- spec$S0
  }
  S1[1] <- S0[1] <- start.S
  I1[1] <- I0[1] <- spec$I0
  R1[1] <- R0[1] <- spec$R0
  inc1[1] <- inc0[1] <- spec$I0

  if (dgp=="SEIR") {
    E1 <- E0 <- numeric(Ttot)
    E1[1] <- E0[1] <- spec$E0
  } else if (!is.null(inf_var)) {
    shape <- inf_mean^2/inf_var
    scale <- inf_var/inf_mean
    F.gamma <- pgamma(seq_len(Ttot), shape=shape, scale=scale)
    pmf <- F.gamma-c(0, head(F.gamma, -1))
  }

  for (tt in 2:Ttot) {
    lambda1 <- beta1[tt]*I1[tt-1]*S1[tt-1]/pop.ind
    lambda0 <- beta0[tt]*I0[tt-1]*S0[tt-1]/pop.ind
    u <- runif(1)
    new1 <- qpois(u, lambda1)
    new0 <- qpois(u, lambda0)
    S1[tt] <- S1[tt-1]-new1
    S0[tt] <- S0[tt-1]-new0

    if (dgp=="SEIR") {
      inc1[tt] <- E1[tt-1]/delta
      inc0[tt] <- E0[tt-1]/delta
      E1[tt] <- (1-1/delta)*E1[tt-1]+new1
      E0[tt] <- (1-1/delta)*E0[tt-1]+new0
      I1[tt] <- (1-1/inf_mean)*I1[tt-1]+inc1[tt]
      I0[tt] <- (1-1/inf_mean)*I0[tt-1]+inc0[tt]
      R1[tt] <- R1[tt-1]+I1[tt-1]/inf_mean
      R0[tt] <- R0[tt-1]+I0[tt-1]/inf_mean
    } else if (is.null(inf_var)) {
      inc1[tt] <- new1
      inc0[tt] <- new0
      I1[tt] <- (1-1/inf_mean)*I1[tt-1]+new1
      I0[tt] <- (1-1/inf_mean)*I0[tt-1]+new0
      R1[tt] <- R1[tt-1]+I1[tt-1]/inf_mean
      R0[tt] <- R0[tt-1]+I0[tt-1]/inf_mean
    } else {
      inc1[tt] <- new1
      inc0[tt] <- new0
      lag.index <- seq_len(tt-1)
      removal1 <- sum(pmf[lag.index]*inc1[tt-lag.index])
      removal0 <- sum(pmf[lag.index]*inc0[tt-lag.index])
      I1[tt] <- I1[tt-1]+new1-removal1
      I0[tt] <- I0[tt-1]+new0-removal0
      R1[tt] <- R1[tt-1]+removal1
      R0[tt] <- R0[tt-1]+removal0
    }
  }

  list(
    treated=data.table(t=seq_len(spec$time_steps), inc=inc1[-1], unit=spec$unit),
    untreated=data.table(t=seq_len(spec$time_steps), inc=inc0[-1], unit=spec$unit)
  )
}

# Optional null-preserving empirical construction.
# difference=TRUE uses exactly:
#   Y.untrt = Y.obs + (Y.model.untrt - Y.model.trt)
construct_observed_difference_counterfactual <- function(observed, modeled_treated, modeled_untreated) {
  observed + (modeled_untreated-modeled_treated)
}

run_inc <- function(data.in, parallel.id) {
  inc.fit <- lm(inc ~ -1 + factor(week) + factor(unit) + factor(trt_post), data=data.in)
  command <- "set matsize 1000
    glm inc i.unit i.week i.trt_post, family(gaussian) link(identity)
    boottest 1.trt_post, cluster(unit) reps(10000) quietly
    gen p = r(p) in 1
    keep p
    keep if _n==1"
  stata.out <- my_RStata(src=command, data.in=data.in, data.out=TRUE, stata.echo=FALSE, id=parallel.id)
  out <- data.frame(model="inc", effect=tail(coef(inc.fit),1), p=stata.out$p,
                    AME=tail(coef(inc.fit),1), AME.adj1=NA, AME.adj2=NA)
  rownames(out) <- NULL
  out
}

run_loginc <- function(data.in, parallel.id) {
  loginc.fit <- glm(inc ~ -1 + factor(week) + factor(unit) + factor(trt_post), family=poisson, data=data.in)
  command <- "set matsize 1000
    glm inc i.unit i.week 1.trt_post, family(poisson) link(log)
    boottest 1.trt_post, cluster(unit) reps(10000) quietly
    gen p = r(p) in 1
    keep p
    keep if _n==1"
  stata.out <- my_RStata(src=command, data.in=data.in, data.out=TRUE, stata.echo=FALSE, id=parallel.id)
  data.untrt <- data.in %>% filter(trt.unit) %>% mutate(trt_post=FALSE)
  data.untrt$loginc_fit <- predict(loginc.fit, newdata=data.untrt, type="response")
  out <- data.frame(model="loginc", effect=exp(tail(coef(loginc.fit),1)), p=stata.out$p,
                    AME=mean(data.in$inc[data.in$trt_post])-mean(data.untrt$loginc_fit[data.untrt$trt.time]),
                    AME.adj1=NA, AME.adj2=NA)
  rownames(out) <- NULL
  out
}

run_growth <- function(data.in, parallel.id=0, trt.IDs=1:N1, coef=NULL) {
  growth.fit <- glm(growth ~ -1 + factor(week) + factor(unit) + factor(trt_post), family=poisson, data=data.in)
  if (is.null(coef)) {
    command <- "set matsize 1000
      glm growth i.unit i.week 1.trt_post, family(poisson) link(log)
      boottest 1.trt_post, cluster(unit) reps(10000) quietly
      gen p = r(p) in 1
      keep p
      keep if _n==1"
    stata.out <- my_RStata(src=command, data.in=data.in, data.out=TRUE, stata.echo=FALSE, id=parallel.id)
  }
  data.untrt <- data.in %>% filter(trt.unit) %>% mutate(trt_post=FALSE)
  data.untrt$growth_fit <- if (is.null(coef)) predict(growth.fit, newdata=data.untrt, type="response") else data.untrt$growth/exp(coef)
  growth.df <- data.untrt %>%
    filter(week>=T0/agg) %>%
    mutate(Y.untrt.growth=ifelse(trt.time, NA, inc)) %>%
    dplyr::select(unit, week, inc, trt.time, Y.untrt.growth, growth_fit)
  for (unit in unique(growth.df$unit)) {
    for (time in (T0/agg+1):((T0+T1)/agg)) {
      growth.df$Y.untrt.growth[growth.df$unit==unit & growth.df$week==time] <-
        growth.df$Y.untrt.growth[growth.df$unit==unit & growth.df$week==time-1] *
        growth.df$growth_fit[growth.df$unit==unit & growth.df$week==time]
    }
  }
  data.untrt$Y.untrt.growth <- NA_real_
  data.untrt$Y.untrt.growth[data.untrt$week>=T0/agg] <- growth.df$Y.untrt.growth
  data.untrt$Y.untrt.growth.adj1 <- data.untrt$Y.untrt.growth.adj2 <- NA_real_
  for (i in trt.IDs) {
    for (time in (T0/agg+1):((T0+T1)/agg)) {
      idx <- data.untrt$unit==i & data.untrt$week==time
      growth.var <- get_var(growth.fit, time, T0=T0/agg, id=i)
      data.untrt$Y.untrt.growth.adj1[idx] <- adjust_fitted_untreated(data.untrt$Y.untrt.growth[idx], growth.var, "approach1")
      data.untrt$Y.untrt.growth.adj2[idx] <- adjust_fitted_untreated(data.untrt$Y.untrt.growth[idx], growth.var, "approach2")
    }
  }
  post <- data.untrt$trt.time
  ame <- mean(data.untrt$inc[post] - data.untrt$Y.untrt.growth[post])
  ame1 <- mean(data.untrt$inc[post] - data.untrt$Y.untrt.growth.adj1[post])
  ame2 <- mean(data.untrt$inc[post] - data.untrt$Y.untrt.growth.adj2[post])
  if (is.null(coef)) {
    data.frame(model="growth", effect=exp(tail(stats::coef(growth.fit),1)), p=stata.out$p,
               AME=ame, AME.adj1=ame1, AME.adj2=ame2)
  } else {
    data.frame(model="growth", effect=exp(coef), AME=ame, AME.adj1=ame1, AME.adj2=ame2)
  }
}

run_Rt <- function(data.in, out.df, type="est", dgp, inf_mean, delta=NULL,
                   inf_var=NULL, trt.IDs=1:N1, coef=NULL, parallel.id=0,
                   unit_population=NULL, incidence_scale=NULL,
                   incidence_aggregation="sum", simulate_from_trt=FALSE,
                   difference=FALSE) {
  if (type=="wt") data.in$Rt <- data.in$R_wt else data.in$Rt <- data.in$R_est
  Rt.fit <- glm(Rt ~ -1 + factor(week) + factor(unit) + factor(trt_post), family=poisson, data=data.in)
  if (is.null(coef)) {
    command <- "set matsize 1000
      glm Rt i.unit i.week 1.trt_post, family(poisson) link(log)
      boottest 1.trt_post, cluster(unit) reps(10000) quietly
      gen p = r(p) in 1
      keep p
      keep if _n==1"
    stata.out <- my_RStata(src=command, data.in=data.in, data.out=TRUE, stata.echo=FALSE, id=parallel.id)
  }
  data.untrt <- data.in %>% filter(trt.unit) %>% mutate(trt_post=FALSE)
  data.untrt$Rt_fit <- if (is.null(coef)) predict(Rt.fit, newdata=data.untrt, type="response") else data.untrt$Rt/exp(coef)
  data.deagg <- out.df %>%
    filter(unit %in% trt.IDs, t<=T0+T1+burnin) %>%
    mutate(week=ceiling(t/agg), unit=factor(unit)) %>%
    merge(data.untrt %>% dplyr::select(unit, week, Rt, Rt_fit) %>% mutate(week=week+burnin/agg),
          by=c("unit","week"), all.x=TRUE) %>%
    mutate(trt.time=t>(T0+burnin), Rt_fit=ifelse(trt.time, Rt_fit, Rt),
           .trans.trt=Rt/S_frac, .trans.untrt=Rt_fit/S_frac) %>%
    filter(t>burnin)

  if (difference) {
    paired <- lapply(trt.IDs, function(ind) {
      pop.ind <- get_unit_population(ind, unit_population)
      spec <- prepare_ame_simulation_window(data.deagg, ind, ".trans.untrt", dgp, simulate_from_trt, ".trans.trt")
      paths <- simulate_paired_ame_trajectories(spec, dgp, pop.ind, inf_mean, delta, inf_var)
      list(
        trt=aggregate_simulated_incidence(paths$treated, incidence_aggregation=incidence_aggregation,
                                          unit_population=unit_population, incidence_scale=incidence_scale,
                                          week_offset=spec$week_offset),
        untrt=aggregate_simulated_incidence(paths$untreated, incidence_aggregation=incidence_aggregation,
                                            unit_population=unit_population, incidence_scale=incidence_scale,
                                            week_offset=spec$week_offset)
      )
    })
    Rt.trt <- rbindlist(lapply(paired, `[[`, "trt"))
    Rt.untrt <- rbindlist(lapply(paired, `[[`, "untrt"))
    data.untrt <- data.untrt %>%
      mutate(.unit_key=as.character(unit)) %>%
      left_join(Rt.trt %>% transmute(.unit_key=as.character(unit), week, .model.trt=inc), by=c(".unit_key","week")) %>%
      left_join(Rt.untrt %>% transmute(.unit_key=as.character(unit), week, .model.untrt=inc), by=c(".unit_key","week")) %>%
      dplyr::select(-.unit_key)
    post <- data.untrt$trt.time
    data.untrt$Y.untrt.Rt <- data.untrt$inc
    data.untrt$Y.untrt.Rt[post] <- construct_observed_difference_counterfactual(
      data.untrt$inc[post], data.untrt$.model.trt[post], data.untrt$.model.untrt[post])
  } else {
    specs <- lapply(trt.IDs, function(ind) prepare_ame_simulation_window(data.deagg, ind, ".trans.untrt", dgp, simulate_from_trt))
    Rt.untrt <- rbindlist(lapply(seq_along(trt.IDs), function(k) {
      simulate_untreated_ame_trajectory(specs[[k]], dgp, get_unit_population(trt.IDs[k], unit_population), inf_mean, delta, inf_var)
    }))
    week.offset <- if (simulate_from_trt) min(data.untrt$week[data.untrt$trt.time])-1 else 0
    Rt.untrt <- aggregate_simulated_incidence(Rt.untrt, trt.IDs, incidence_aggregation,
                                              unit_population, incidence_scale, week.offset)
    data.untrt <- data.untrt %>%
      mutate(.unit_key=as.character(unit)) %>%
      left_join(Rt.untrt %>% transmute(.unit_key=as.character(unit), week, Y.untrt.Rt=inc), by=c(".unit_key","week")) %>%
      dplyr::select(-.unit_key)
    post <- data.untrt$trt.time
  }

  data.untrt$Y.untrt.Rt.adj1 <- data.untrt$Y.untrt.Rt.adj2 <- NA_real_
  for (i in trt.IDs) {
    for (time in (T0/agg+1):((T0+T1)/agg)) {
      idx <- as.character(data.untrt$unit)==as.character(i) & data.untrt$week==time
      Rt.var <- get_var(Rt.fit, time, T0=T0/agg, id=i)
      data.untrt$Y.untrt.Rt.adj1[idx] <- adjust_fitted_untreated(data.untrt$Y.untrt.Rt[idx], Rt.var, "approach1")
      data.untrt$Y.untrt.Rt.adj2[idx] <- adjust_fitted_untreated(data.untrt$Y.untrt.Rt[idx], Rt.var, "approach2")
    }
  }

  ame <- mean(data.untrt$inc[post] - data.untrt$Y.untrt.Rt[post])
  ame1 <- mean(data.untrt$inc[post] - data.untrt$Y.untrt.Rt.adj1[post])
  ame2 <- mean(data.untrt$inc[post] - data.untrt$Y.untrt.Rt.adj2[post])
  if (is.null(coef)) {
    data.frame(model=paste0("Rt_",type), effect=exp(tail(stats::coef(Rt.fit),1)), p=stata.out$p,
               AME=ame, AME.adj1=ame1, AME.adj2=ame2)
  } else {
    data.frame(model=paste0("Rt_",type), effect=exp(coef),
               AME=ame, AME.adj1=ame1, AME.adj2=ame2)
  }
}

run_beta <- function(data.in, out.df, dgp, inf_mean, delta=NULL, inf_var=NULL,
                     trt.IDs=1:N1, coef=NULL, parallel.id=0,
                     unit_population=NULL, incidence_scale=NULL,
                     incidence_aggregation="sum", simulate_from_trt=FALSE,
                     difference=FALSE) {
  beta.fit <- glm(beta_est ~ -1 + factor(week) + factor(unit) + factor(trt_post), family=poisson, data=data.in)
  if (is.null(coef)) {
    command <- "set matsize 600
      glm beta_est i.unit i.week 1.trt_post, family(poisson) link(log)
      boottest 1.trt_post, cluster(unit) reps(10000) quietly
      gen p = r(p) in 1
      keep p
      keep if _n==1"
    stata.out <- my_RStata(src=command, data.in=data.in, data.out=TRUE, stata.echo=FALSE, id=parallel.id)
  }
  data.untrt <- data.in %>% filter(trt.unit) %>% mutate(trt_post=FALSE)
  data.untrt$beta_fit <- if (is.null(coef)) predict(beta.fit, newdata=data.untrt, type="response") else data.untrt$beta_est/exp(coef)
  data.deagg <- out.df %>%
    filter(unit %in% trt.IDs, t<=T0+T1+burnin) %>%
    mutate(week=ceiling(t/agg), unit=factor(unit)) %>%
    merge(data.untrt %>% dplyr::select(unit, week, beta_est, beta_fit) %>% mutate(week=week+burnin/agg),
          by=c("unit","week"), all.x=TRUE) %>%
    mutate(trt.time=t>(T0+burnin), beta_fit=ifelse(trt.time, beta_fit, beta_est)) %>%
    filter(t>burnin)

  if (difference) {
    paired <- lapply(trt.IDs, function(ind) {
      pop.ind <- get_unit_population(ind, unit_population)
      spec <- prepare_ame_simulation_window(data.deagg, ind, "beta_fit", dgp, simulate_from_trt, "beta_est")
      paths <- simulate_paired_ame_trajectories(spec, dgp, pop.ind, inf_mean, delta, inf_var)
      list(
        trt=aggregate_simulated_incidence(paths$treated, incidence_aggregation=incidence_aggregation,
                                          unit_population=unit_population, incidence_scale=incidence_scale,
                                          week_offset=spec$week_offset),
        untrt=aggregate_simulated_incidence(paths$untreated, incidence_aggregation=incidence_aggregation,
                                            unit_population=unit_population, incidence_scale=incidence_scale,
                                            week_offset=spec$week_offset)
      )
    })
    beta.trt <- rbindlist(lapply(paired, `[[`, "trt"))
    beta.untrt <- rbindlist(lapply(paired, `[[`, "untrt"))
    data.untrt <- data.untrt %>%
      mutate(.unit_key=as.character(unit)) %>%
      left_join(beta.trt %>% transmute(.unit_key=as.character(unit), week, .model.trt=inc), by=c(".unit_key","week")) %>%
      left_join(beta.untrt %>% transmute(.unit_key=as.character(unit), week, .model.untrt=inc), by=c(".unit_key","week")) %>%
      dplyr::select(-.unit_key)
    post <- data.untrt$trt.time
    data.untrt$Y.untrt.beta <- data.untrt$inc
    data.untrt$Y.untrt.beta[post] <- construct_observed_difference_counterfactual(
      data.untrt$inc[post], data.untrt$.model.trt[post], data.untrt$.model.untrt[post])
  } else {
    specs <- lapply(trt.IDs, function(ind) prepare_ame_simulation_window(data.deagg, ind, "beta_fit", dgp, simulate_from_trt))
    beta.untrt <- rbindlist(lapply(seq_along(trt.IDs), function(k) {
      simulate_untreated_ame_trajectory(specs[[k]], dgp, get_unit_population(trt.IDs[k], unit_population), inf_mean, delta, inf_var)
    }))
    week.offset <- if (simulate_from_trt) min(data.untrt$week[data.untrt$trt.time])-1 else 0
    beta.untrt <- aggregate_simulated_incidence(beta.untrt, trt.IDs, incidence_aggregation,
                                                unit_population, incidence_scale, week.offset)
    data.untrt <- data.untrt %>%
      mutate(.unit_key=as.character(unit)) %>%
      left_join(beta.untrt %>% transmute(.unit_key=as.character(unit), week, Y.untrt.beta=inc), by=c(".unit_key","week")) %>%
      dplyr::select(-.unit_key)
    post <- data.untrt$trt.time
  }

  data.untrt$Y.untrt.beta.adj1 <- data.untrt$Y.untrt.beta.adj2 <- NA_real_
  for (i in trt.IDs) {
    for (time in (T0/agg+1):((T0+T1)/agg)) {
      idx <- as.character(data.untrt$unit)==as.character(i) & data.untrt$week==time
      beta.var <- get_var(beta.fit, time, T0=T0/agg, id=i)
      data.untrt$Y.untrt.beta.adj1[idx] <- adjust_fitted_untreated(data.untrt$Y.untrt.beta[idx], beta.var, "approach1")
      data.untrt$Y.untrt.beta.adj2[idx] <- adjust_fitted_untreated(data.untrt$Y.untrt.beta[idx], beta.var, "approach2")
    }
  }

  ame <- mean(data.untrt$inc[post] - data.untrt$Y.untrt.beta[post])
  ame1 <- mean(data.untrt$inc[post] - data.untrt$Y.untrt.beta.adj1[post])
  ame2 <- mean(data.untrt$inc[post] - data.untrt$Y.untrt.beta.adj2[post])
  if (is.null(coef)) {
    data.frame(model="beta", effect=exp(tail(stats::coef(beta.fit),1)), p=stata.out$p,
               AME=ame, AME.adj1=ame1, AME.adj2=ame2)
  } else {
    data.frame(model="beta", effect=exp(coef),
               AME=ame, AME.adj1=ame1, AME.adj2=ame2)
  }
}

# Calculate the true untreated incidence counterfactual
run_true <- function(out.df, trans_prob.base1, dgp, trt.IDs=1:N1) {
  true.untrt <- rbindlist(lapply(trt.IDs, function(ind) {
    I0 <- out.df$I[out.df$unit==ind & out.df$t==burnin+1]
    recovered <- out.df$R[out.df$unit==ind & out.df$t==burnin+1]
    if (dgp=="SIR") {
      run_SIR_varying(pop.size=pop.size, seeds=I0, recovered=recovered, time_steps=(T0+T1), 
                      inf_mean=inf_mean, trans_prob=rep(trans_prob.base1, (T0+T1)))
    } else if (dgp=="SEIR") {
      E0 <- out.df$E[out.df$unit==ind & out.df$t==burnin+1]
      run_SEIR_varying(pop.size=pop.size, I0=I0, E0=E0, recovered=recovered, time_steps=(T0+T1), 
                       inf_mean=inf_mean, delta=delta, trans_prob=rep(trans_prob.base1, (T0+T1)))
    }
  })) %>%
    mutate(unit=rep(trt.IDs, each=(T0+T1)), week = ceiling(t/agg)) %>%
    group_by(unit, week) %>%
    summarise(inc = sum(inc))
  Y.untrt.true <- mean(true.untrt$inc[true.untrt$week > (T0/agg)])
  return(Y.untrt.true)
}

# Helper functions for School Masking example
loginc_AME <- function(coef, subset=NULL) {
  ATT_gt.tmp <- data.frame(time_to_trt = time_to_trt, coef = coef)
  loginc.obs <- df_sunab %>% mutate(time_to_trt = time - group)
  
  if (!is.null(subset)) {
    loginc.obs <- loginc.obs %>% filter(time_to_trt %in% subset)
  } else {
    loginc.obs <- loginc.obs %>% filter(time_to_trt >= 0)
  }
    
  loginc.obs <- loginc.obs %>%
    merge(ATT_gt.tmp, by = "time_to_trt") %>%
    mutate(y.ctl = y/exp(coef), # fitted control potential outcome on case scale, cannot convert on the original log scale b/c of zeros
           diff = y - y.ctl) %>% # difference to get marginal effect for each unit i to ME
    group_by(time_to_trt) %>%
    summarise(diff = mean(diff)) # average over all units to get the average marginal effect
  sum(loginc.obs$diff)
}

growth_AME <- function(coef, subset=NULL) {
  ATT_gt.tmp <- data.frame(time_to_trt = time_to_trt, coef = coef)
  
  df.all <- df_sunab %>% mutate(time_to_trt = time - group)
  df <- df.all %>%
    filter(time_to_trt >= -1) %>%
    merge(ATT_gt.tmp, by = "time_to_trt", all.x = T) %>%
    mutate(coef = ifelse(time_to_trt==-1, 0, coef),
           growth.ctl = log(y) - coef,
           Pos.ctl = NA)
  
  for (unit in unique(df$ID)) {
    for (time in sort(unique(df$time[df$ID==unit]))) {
      if ((time+1) == unique(df$group[df$ID==unit])) { # last period before intervention
        inc.last <- df$PosPer1K[df$ID==unit & df$time==time]
        if (inc.last==0) {
          df$Pos.ctl[df$ID==unit & df$time==time] <- inc.last
        } else {
          df$Pos.ctl[df$ID==unit & df$time==time] <- mean(df.all$PosPer1K[df.all$ID==unit & df.all$time %in% (time-4):time])
        }
      } else { # recover the untreated trajectory for the treated group
        time.ind <- ifelse(nrow(df[df$ID==unit & df$time==time-1,])==0, 2, 1)
        df$Pos.ctl[df$ID==unit & df$time==time] <- df$Pos.ctl[df$ID==unit & df$time==time-time.ind] * 
          exp(df$growth.ctl[df$ID==unit & df$time==time])
      }
    }
  }
  
  if (!is.null(subset)) df <- df %>% filter(time_to_trt %in% subset)
  
  df <- df %>% filter(time >= group) %>%
    mutate(diff = PosPer1K - Pos.ctl) %>%
    group_by(ID) %>% summarise(diff = sum(diff)) # difference to get marginal effect for each unit
  mean(df$diff) # average over all units to get the AME
}
