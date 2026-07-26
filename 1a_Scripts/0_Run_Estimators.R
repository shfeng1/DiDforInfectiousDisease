# The helper functions below run each of the estimators in the following steps
# 1. Fit lm or glm model
# 2. get p-value from calling Stata using wild score cluster bootstrap
# 3. Reover AMEs
# 4. Bias correction if needed
# 5. Summarize and return results

source("./1a_Scripts/0_RStata.R")
source("./1a_Scripts/0_Bias_Correction.R")
source("./1a_Scripts/0_Smearing_Correction.R")

# Resolve the population used by a mechanistic AME simulation.
# Defaults preserve the simulation-study behavior, which uses the global scalar
# pop.size.  Kansas passes a named vector so each county is simulated at its
# actual population rather than at the unrelated simulation-study population.
get_unit_population <- function(ind, unit_population = NULL) {
  if (is.null(unit_population)) {
    return(pop.size)
  }
  if (is.null(names(unit_population))) {
    stop("unit_population must be a named numeric vector keyed by unit ID.")
  }

  value <- unname(unit_population[as.character(ind)])
  if (length(value) != 1 || is.na(value) || !is.finite(value) || value <= 0) {
    stop("Missing or invalid population for unit ", as.character(ind), ".")
  }
  value
}

# Convert simulated unit-level counts to the manuscript's requested incidence
# scale.  For Kansas, incidence_aggregation='mean' produces the weekly mean of
# daily infections (matching the observed weekly outcome) and incidence_scale
# converts county counts to cases per 100,000.  Defaults retain weekly sums in
# the simulation study.
aggregate_simulated_incidence <- function(sim.data, trt.IDs,
                                          incidence_aggregation = c("sum", "mean"),
                                          unit_population = NULL,
                                          incidence_scale = NULL) {
  incidence_aggregation <- match.arg(incidence_aggregation)

  out <- sim.data %>%
    mutate(unit = rep(trt.IDs, each = (T0 + T1)), week = ceiling(t / agg)) %>%
    group_by(unit, week) %>%
    summarise(
      inc = if (incidence_aggregation == "sum") sum(inc) else mean(inc),
      .groups = "drop"
    )

  if (!is.null(incidence_scale)) {
    if (is.null(unit_population)) {
      stop("incidence_scale requires unit_population.")
    }
    population <- unname(unit_population[as.character(out$unit)])
    if (anyNA(population) || any(!is.finite(population)) || any(population <= 0)) {
      stop("Cannot standardize simulated incidence because a unit population is missing or invalid.")
    }
    out$inc <- out$inc / population * incidence_scale
  }

  out
}

# Central helpers for the manuscript-aligned AME calculation.
#
# The fitted untreated incidence is bias-corrected FIRST.  The treated/observed
# outcome is never used in the bias correction itself:
#
#   Approach 1: Y.untrt.adj1 = Y.untrt / (1 + variance / 2)
#   Approach 2: Y.untrt.adj2 = Y.untrt / exp(variance / 2)
#
# Only after the untreated counterfactual has been corrected do we calculate:
#
#   AME      = Y.obs - Y.untrt
#   AME.adj1 = Y.obs - Y.untrt.adj1
#   AME.adj2 = Y.obs - Y.untrt.adj2
#
# Y.obs depends on the existing sim handler:
#   sim = FALSE: observed treated incidence from the input data;
#   sim = TRUE:  separately simulated treated incidence.
ame_observed_minus_untreated <- function(observed, untreated) {
  if (length(observed) != length(untreated) && length(observed) != 1L && length(untreated) != 1L) {
    stop("Observed and untreated outcomes must have compatible lengths.")
  }
  observed - untreated
}

adjust_fitted_untreated <- function(untreated, variance,
                                    approach = c("approach1", "approach2")) {
  approach <- match.arg(approach)

  if (any(!is.finite(untreated)) || any(!is.finite(variance))) {
    stop("Untreated outcomes and bias-correction variances must be finite.")
  }
  if (any(variance < -1e-10)) {
    stop("A negative cumulative variance was supplied for bias correction.")
  }

  # Permit only negligible negative values caused by floating-point arithmetic.
  variance <- pmax(variance, 0)
  correction_factor <- if (approach == "approach1") {
    1 + variance / 2
  } else {
    exp(variance / 2)
  }

  untreated / correction_factor
}

mean_selected <- function(x, selected, label) {
  if (length(x) != length(selected)) {
    stop(label, " and its selection indicator have different lengths.")
  }
  values <- x[which(selected %in% TRUE)]
  if (length(values) == 0L) {
    stop("No post-intervention values were available for ", label, ".")
  }
  if (anyNA(values) || any(!is.finite(values))) {
    stop("Missing or non-finite post-intervention values were found for ", label, ".")
  }
  mean(values)
}

# Every estimator returns the underlying outcomes as well as the AMEs.  Summary
# scripts therefore read Y.untrt/Y.untrt.adj1/Y.untrt.adj2 directly and never
# reverse-engineer a counterfactual from a stored AME.
build_ame_output <- function(model, effect, p, observed, untreated,
                             untreated_adj1 = untreated,
                             untreated_adj2 = untreated,
                             untreated_smear = NULL,
                             smear_factor = NA_real_,
                             smear_ratio = NA_real_,
                             post_index = rep(TRUE, length(observed))) {
  Y.obs <- mean_selected(observed, post_index, "Y.obs")
  Y.untrt <- mean_selected(untreated, post_index, "Y.untrt")
  Y.untrt.adj1 <- mean_selected(untreated_adj1, post_index, "Y.untrt.adj1")
  Y.untrt.adj2 <- mean_selected(untreated_adj2, post_index, "Y.untrt.adj2")

  if (is.null(untreated_smear)) {
    Y.untrt.smear <- NA_real_
    AME.smear <- NA_real_
  } else {
    Y.untrt.smear <- mean_selected(untreated_smear, post_index, "Y.untrt.smear")
    AME.smear <- ame_observed_minus_untreated(Y.obs, Y.untrt.smear)
  }

  data.frame(
    model = model,
    effect = effect,
    p = p,
    Y.obs = Y.obs,
    Y.untrt = Y.untrt,
    Y.untrt.adj1 = Y.untrt.adj1,
    Y.untrt.adj2 = Y.untrt.adj2,
    Y.untrt.smear = Y.untrt.smear,
    AME = ame_observed_minus_untreated(Y.obs, Y.untrt),
    AME.adj1 = ame_observed_minus_untreated(Y.obs, Y.untrt.adj1),
    AME.adj2 = ame_observed_minus_untreated(Y.obs, Y.untrt.adj2),
    AME.smear = AME.smear,
    smear.factor = smear_factor,
    smear.ratio = smear_ratio,
    row.names = NULL
  )
}

run_inc <- function(data.in, parallel.id) {
  inc.fit <- lm(inc ~ -1 + factor(week) + factor(unit) + factor(trt_post), data = data.in)

  command <- "set matsize 1000
    glm inc i.unit i.week i.trt_post, family(gaussian) link(identity)
    boottest 1.trt_post, cluster(unit) reps(10000) quietly
    gen p = r(p) in 1
    keep p
    keep if _n==1"
  stata.out <- my_RStata(src=command, data.in=data.in, data.out=TRUE,
                         stata.echo=FALSE, id=parallel.id)

  data.untrt <- data.in %>% filter(trt.unit) %>% mutate(trt_post = FALSE)
  data.untrt$Y.untrt <- predict(inc.fit, newdata = data.untrt, type = "response")

  # No log-normal transformation correction applies to the incidence model, so
  # both adjusted untreated columns equal the fitted untreated outcome.
  build_ame_output(
    model = "inc",
    effect = tail(coef(inc.fit), 1),
    p = stata.out$p,
    observed = data.untrt$inc,
    untreated = data.untrt$Y.untrt,
    post_index = data.untrt$trt.time
  )
}

run_loginc <- function(data.in, parallel.id) {
  loginc.fit <- glm(inc ~ -1 + factor(week) + factor(unit) + factor(trt_post),
                    family = poisson, data = data.in)

  command <- "set matsize 1000
    glm inc i.unit i.week 1.trt_post, family(poisson) link(log)
    boottest 1.trt_post, cluster(unit) reps(10000) quietly
    gen p = r(p) in 1
    keep p
    keep if _n==1"
  stata.out <- my_RStata(src=command, data.in=data.in, data.out=TRUE,
                         stata.echo=FALSE, id=parallel.id)

  data.untrt <- data.in %>% filter(trt.unit) %>% mutate(trt_post = FALSE)
  data.untrt$Y.untrt <- predict(loginc.fit, newdata = data.untrt, type = "response")

  # The manuscript's recursive bias correction is for log growth, log Rt, and
  # log beta.  Accordingly, log-incidence adjusted columns equal Y.untrt.
  build_ame_output(
    model = "loginc",
    effect = exp(tail(coef(loginc.fit), 1)),
    p = stata.out$p,
    observed = data.untrt$inc,
    untreated = data.untrt$Y.untrt,
    post_index = data.untrt$trt.time
  )
}

run_growth <- function(data.in, parallel.id=0, trt.IDs=1:N1, coef=NULL, sim=FALSE) {
  growth.fit <- glm(growth ~ -1 + factor(week) + factor(unit) + factor(trt_post),
                    family = poisson, data = data.in)
  command <- "set matsize 1000
    glm growth i.unit i.week 1.trt_post, family(poisson) link(log)
    boottest 1.trt_post, cluster(unit) reps(10000) quietly
    gen p = r(p) in 1
    keep p
    keep if _n==1"
  stata.out <- my_RStata(src=command, data.in=data.in, data.out=TRUE,
                         stata.echo=FALSE, id=parallel.id)

  data.untrt <- data.in %>% filter(trt.unit) %>% mutate(trt_post = FALSE)
  if (is.null(coef)) {
    data.untrt$growth_fit <- predict(growth.fit, newdata = data.untrt,
                                     type = "response")
  } else {
    data.untrt$growth_fit <- data.untrt$growth / exp(coef)
  }

  growth.df <- data.untrt %>%
    filter(trt.unit, week >= T0/agg) %>%
    mutate(
      trt_post = trt.unit & trt.time,
      Y.untrt.growth = ifelse(trt_post, NA, inc),
      Y.trt.growth = ifelse(trt_post, NA, inc)
    ) %>%
    dplyr::select(unit, week, S_frac, trt.unit, trt.time, trt_post,
                  inc, Y.untrt.growth, growth_fit, Y.trt.growth, growth)

  for (unit in unique(growth.df$unit)) {
    for (time in (T0/agg+1):((T0+T1)/agg)) {
      growth.df$Y.untrt.growth[growth.df$unit==unit & growth.df$week==time] <-
        growth.df$Y.untrt.growth[growth.df$unit==unit & growth.df$week==time-1] *
        growth.df$growth_fit[growth.df$unit==unit & growth.df$week==time]
      growth.df$Y.trt.growth[growth.df$unit==unit & growth.df$week==time] <-
        growth.df$Y.trt.growth[growth.df$unit==unit & growth.df$week==time-1] *
        growth.df$growth[growth.df$unit==unit & growth.df$week==time]
    }
  }

  data.untrt$Y.untrt.growth <- NA_real_
  data.untrt$Y.trt.growth <- NA_real_
  data.untrt$Y.untrt.growth[data.untrt$week >= T0/agg] <- growth.df$Y.untrt.growth
  data.untrt$Y.trt.growth[data.untrt$week >= T0/agg] <- growth.df$Y.trt.growth
  data.untrt$Y.untrt.growth.adj1 <- NA_real_
  data.untrt$Y.untrt.growth.adj2 <- NA_real_

  # Manuscript-aligned correction: adjust the fitted untreated trajectory itself.
  for (i in trt.IDs) {
    for (time in (T0/agg+1):((T0+T1)/agg)) {
      idx <- data.untrt$unit==i & data.untrt$week==time
      growth.var <- get_var(growth.fit, time, T0=T0/agg, id=i)
      data.untrt$Y.untrt.growth.adj1[idx] <- adjust_fitted_untreated(
        data.untrt$Y.untrt.growth[idx], growth.var, "approach1")
      data.untrt$Y.untrt.growth.adj2[idx] <- adjust_fitted_untreated(
        data.untrt$Y.untrt.growth[idx], growth.var, "approach2")
    }
  }

  observed <- if (sim) data.untrt$Y.trt.growth else data.untrt$inc
  build_ame_output(
    model = "growth",
    effect = exp(tail(coef(growth.fit), 1)),
    p = stata.out$p,
    observed = observed,
    untreated = data.untrt$Y.untrt.growth,
    untreated_adj1 = data.untrt$Y.untrt.growth.adj1,
    untreated_adj2 = data.untrt$Y.untrt.growth.adj2,
    post_index = data.untrt$trt.time
  )
}

# type is either "wt" for cohort-based or "est" for prevalence-based.
# dgp is either "SIR" or "SEIR".
run_Rt <- function(data.in, out.df, type="est", dgp, inf_mean, delta=NULL,
                   inf_var=NULL, trt.IDs=1:N1, coef=NULL, parallel.id=0,
                   sim=FALSE, unit_population=NULL, incidence_scale=NULL,
                   incidence_aggregation=c("sum", "mean"),
                   smearing=FALSE, smearing_reps=500L,
                   smearing_method=c("local", "cumulative"),
                   p_value=NULL) {
  incidence_aggregation <- match.arg(incidence_aggregation)
  smearing_method <- match.arg(smearing_method)
  if (type == "wt") {
    data.in$Rt <- data.in$R_wt
  } else if (type == "est") {
    data.in$Rt <- data.in$R_est
  }

  Rt.fit <- glm(Rt ~ -1 + factor(week) + factor(unit) + factor(trt_post),
                family = poisson, data = data.in)
  data.untrt <- data.in %>% filter(trt.unit) %>% mutate(trt_post = FALSE)
  command <- "set matsize 1000
    glm Rt i.unit i.week 1.trt_post, family(poisson) link(log)
    boottest 1.trt_post, cluster(unit) reps(10000) quietly
    gen p = r(p) in 1
    keep p
    keep if _n==1"
  stata.out <- if (is.null(p_value)) {
    my_RStata(src=command, data.in=data.in, data.out=TRUE,
              stata.echo=FALSE, id=parallel.id)
  } else {
    if (length(p_value) != 1L || !is.finite(p_value)) {
      stop("p_value must be NULL or one finite numeric value.")
    }
    data.frame(p = as.numeric(p_value))
  }

  if (is.null(coef)) {
    data.untrt$Rt_fit <- predict(Rt.fit, newdata = data.untrt,
                                 type = "response")
  } else {
    data.untrt$Rt_fit <- data.untrt$Rt / exp(coef)
  }

  data.untrt.deagg <- out.df %>%
    filter(unit %in% trt.IDs, t <= T0 + T1 + burnin) %>%
    mutate(week = ceiling(t/agg), unit = factor(unit)) %>%
    merge(data.untrt %>%
            dplyr::select(unit, week, Rt, Rt_fit) %>%
            mutate(week = week + burnin/agg),
          by = c("unit", "week"), all.x = TRUE) %>%
    mutate(trt.time = t > (T0+burnin),
           Rt_fit = ifelse(trt.time, Rt_fit, Rt)) %>%
    filter(t > burnin)

  if (smearing && smearing_method == "local") {
    data.untrt.deagg <- data.untrt.deagg %>%
      mutate(.trans.trt = Rt / S_frac,
             .trans.untrt = Rt_fit / S_frac)
    Rt.local <- run_local_paired_mc_periods(
      data_deagg = data.untrt.deagg,
      dgp = dgp,
      reps = smearing_reps,
      trt.IDs = trt.IDs,
      treated_rate = ".trans.trt",
      untreated_rate = ".trans.untrt",
      inf_mean = inf_mean,
      delta = delta,
      inf_var = inf_var,
      agg = agg,
      unit_population = unit_population,
      incidence_scale = incidence_scale,
      incidence_aggregation = incidence_aggregation
    )
    Rt.trt <- Rt.local %>% transmute(unit, week, inc = modeled_treated)
    Rt.untrt <- Rt.local %>% transmute(unit, week, inc = modeled_untreated)
  } else if (smearing) {
    Rt.paired <- lapply(trt.IDs, function(ind) {
      pop.ind <- get_unit_population(ind, unit_population)
      I0 <- data.untrt.deagg$I[data.untrt.deagg$unit==ind &
                                 data.untrt.deagg$t==burnin+1]
      recovered <- data.untrt.deagg$R[data.untrt.deagg$unit==ind &
                                        data.untrt.deagg$t==burnin+1]
      S_frac <- data.untrt.deagg$S_frac[data.untrt.deagg$unit==ind]
      trans_prob.trt <- data.untrt.deagg$Rt[data.untrt.deagg$unit==ind] / S_frac
      trans_prob.untrt <- data.untrt.deagg$Rt_fit[data.untrt.deagg$unit==ind] / S_frac
      E0 <- if (dgp == "SEIR") {
        data.untrt.deagg$E[data.untrt.deagg$unit==ind &
                             data.untrt.deagg$t==burnin+1]
      } else NULL

      run_paired_mc_trajectories(
        dgp=dgp, reps=smearing_reps, time_steps=(T0+T1),
        pop.size=pop.ind, I0=I0, E0=E0, recovered=recovered,
        trans_prob_treated=trans_prob.trt,
        trans_prob_untreated=trans_prob.untrt,
        inf_mean=inf_mean, inf_var=inf_var, delta=delta
      )
    })
    Rt.trt <- rbindlist(lapply(Rt.paired, `[[`, "treated"))
    Rt.untrt <- rbindlist(lapply(Rt.paired, `[[`, "untreated"))
  } else {
    Rt.untrt <- rbindlist(lapply(trt.IDs, function(ind) {
      pop.ind <- get_unit_population(ind, unit_population)
      I0 <- data.untrt.deagg$I[data.untrt.deagg$unit==ind &
                                 data.untrt.deagg$t==burnin+1]
      recovered <- data.untrt.deagg$R[data.untrt.deagg$unit==ind &
                                        data.untrt.deagg$t==burnin+1]
      S_frac <- data.untrt.deagg$S_frac[data.untrt.deagg$unit==ind]
      trans_prob <- data.untrt.deagg$Rt_fit[data.untrt.deagg$unit==ind] / S_frac
      if (dgp == "SIR") {
        run_SIR_varying(pop.size=pop.ind, seeds=I0, recovered=recovered,
                        trans_prob=trans_prob, time_steps=(T0+T1),
                        inf_mean=inf_mean, inf_var=inf_var)
      } else if (dgp == "SEIR") {
        E0 <- data.untrt.deagg$E[data.untrt.deagg$unit==ind &
                                   data.untrt.deagg$t==burnin+1]
        run_SEIR_varying(pop.size=pop.ind, I0=I0, E0=E0, recovered=recovered,
                         trans_prob=trans_prob, time_steps=(T0+T1),
                         inf_mean=inf_mean, delta=delta)
      }
    }))

    Rt.trt <- rbindlist(lapply(trt.IDs, function(ind) {
      pop.ind <- get_unit_population(ind, unit_population)
      I0 <- data.untrt.deagg$I[data.untrt.deagg$unit==ind &
                                 data.untrt.deagg$t==burnin+1]
      recovered <- data.untrt.deagg$R[data.untrt.deagg$unit==ind &
                                        data.untrt.deagg$t==burnin+1]
      S_frac <- data.untrt.deagg$S_frac[data.untrt.deagg$unit==ind]
      trans_prob <- data.untrt.deagg$Rt[data.untrt.deagg$unit==ind] / S_frac
      if (dgp == "SIR") {
        run_SIR_varying(pop.size=pop.ind, seeds=I0, recovered=recovered,
                        trans_prob=trans_prob, time_steps=(T0+T1),
                        inf_mean=inf_mean, inf_var=inf_var)
      } else if (dgp == "SEIR") {
        E0 <- data.untrt.deagg$E[data.untrt.deagg$unit==ind &
                                   data.untrt.deagg$t==burnin+1]
        run_SEIR_varying(pop.size=pop.ind, I0=I0, E0=E0, recovered=recovered,
                         trans_prob=trans_prob, time_steps=(T0+T1),
                         inf_mean=inf_mean, delta=delta)
      }
    }))
  }

  if (!(smearing && smearing_method == "local")) {
    Rt.untrt <- aggregate_simulated_incidence(
      Rt.untrt, trt.IDs,
      incidence_aggregation=incidence_aggregation,
      unit_population=unit_population,
      incidence_scale=incidence_scale
    )
    Rt.trt <- aggregate_simulated_incidence(
      Rt.trt, trt.IDs,
      incidence_aggregation=incidence_aggregation,
      unit_population=unit_population,
      incidence_scale=incidence_scale
    )
  }

  # Join by explicit unit/week keys.  Vector assignment could silently misalign
  # trajectories if a unit or week is reordered.
  data.untrt <- data.untrt %>%
    mutate(.unit_key = as.character(unit)) %>%
    left_join(Rt.untrt %>%
                transmute(.unit_key=as.character(unit), week,
                          Y.untrt.Rt=inc),
              by=c(".unit_key", "week")) %>%
    left_join(Rt.trt %>%
                transmute(.unit_key=as.character(unit), week,
                          Y.trt.Rt=inc),
              by=c(".unit_key", "week")) %>%
    dplyr::select(-.unit_key) %>%
    mutate(Y.untrt.Rt.adj1=NA_real_, Y.untrt.Rt.adj2=NA_real_)

  for (i in trt.IDs) {
    for (time in (T0/agg+1):((T0+T1)/agg)) {
      idx <- data.untrt$unit==i & data.untrt$week==time
      Rt.var <- get_var(Rt.fit, time, T0=T0/agg, id=i)
      data.untrt$Y.untrt.Rt.adj1[idx] <- adjust_fitted_untreated(
        data.untrt$Y.untrt.Rt[idx], Rt.var, "approach1")
      data.untrt$Y.untrt.Rt.adj2[idx] <- adjust_fitted_untreated(
        data.untrt$Y.untrt.Rt[idx], Rt.var, "approach2")
    }
  }

  observed <- if (sim) data.untrt$Y.trt.Rt else data.untrt$inc
  smear.out <- NULL
  if (smearing) {
    effect.log <- if (is.null(coef)) as.numeric(tail(coef(Rt.fit), 1)) else as.numeric(coef)
    smear.out <- if (smearing_method == "local") {
      smear_counterfactual_local(
        observed = observed,
        modeled_treated = data.untrt$Y.trt.Rt,
        modeled_untreated = data.untrt$Y.untrt.Rt,
        unit = data.untrt$unit,
        week = data.untrt$week,
        post_index = data.untrt$trt.time,
        log_effect = effect.log
      )
    } else {
      smear_counterfactual_trajectory(
        observed = observed,
        modeled_treated = data.untrt$Y.trt.Rt,
        modeled_untreated = data.untrt$Y.untrt.Rt,
        unit = data.untrt$unit,
        post_index = data.untrt$trt.time,
        log_effect = effect.log
      )
    }
  }

  build_ame_output(
    model = paste0("Rt_", type),
    effect = exp(tail(coef(Rt.fit), 1)),
    p = stata.out$p,
    observed = observed,
    untreated = data.untrt$Y.untrt.Rt,
    untreated_adj1 = data.untrt$Y.untrt.Rt.adj1,
    untreated_adj2 = data.untrt$Y.untrt.Rt.adj2,
    untreated_smear = if (is.null(smear.out)) NULL else smear.out$untreated,
    smear_factor = if (is.null(smear.out)) NA_real_ else mean(smear.out$detail$smear.factor),
    smear_ratio = if (is.null(smear.out)) NA_real_ else smearing_weighted_ratio(smear.out$detail),
    post_index = data.untrt$trt.time
  )
}

run_beta <- function(data.in, out.df, dgp, inf_mean, delta=NULL, inf_var=NULL,
                     trt.IDs=1:N1, coef=NULL, parallel.id=0, sim=FALSE,
                     unit_population=NULL, incidence_scale=NULL,
                     incidence_aggregation=c("sum", "mean"),
                     smearing=FALSE, smearing_reps=500L,
                     smearing_method=c("local", "cumulative"),
                     p_value=NULL) {
  incidence_aggregation <- match.arg(incidence_aggregation)
  smearing_method <- match.arg(smearing_method)
  beta.fit <- glm(beta_est ~ -1 + factor(week) + factor(unit) + factor(trt_post),
                  family = poisson, data = data.in)
  data.untrt <- data.in %>% filter(trt.unit) %>% mutate(trt_post = FALSE)
  command <- "set matsize 600
    glm beta_est i.unit i.week 1.trt_post, family(poisson) link(log)
    boottest 1.trt_post, cluster(unit) reps(10000) quietly
    gen p = r(p) in 1
    keep p
    keep if _n==1"
  stata.out <- if (is.null(p_value)) {
    my_RStata(src=command, data.in=data.in, data.out=TRUE,
              stata.echo=FALSE, id=parallel.id)
  } else {
    if (length(p_value) != 1L || !is.finite(p_value)) {
      stop("p_value must be NULL or one finite numeric value.")
    }
    data.frame(p = as.numeric(p_value))
  }

  if (is.null(coef)) {
    data.untrt$beta_fit <- predict(beta.fit, newdata = data.untrt,
                                   type = "response")
  } else {
    data.untrt$beta_fit <- data.untrt$beta_est / exp(coef)
  }

  data.deagg <- out.df %>%
    filter(unit %in% trt.IDs, t <= T0 + T1 + burnin) %>%
    mutate(week = ceiling(t/agg), unit = factor(unit)) %>%
    merge(data.untrt %>%
            dplyr::select(unit, week, beta_est, beta_fit) %>%
            mutate(week = week + burnin/agg),
          by = c("unit", "week"), all.x = TRUE) %>%
    mutate(trt.time = t > (T0+burnin),
           beta_fit = ifelse(trt.time, beta_fit, beta_est)) %>%
    filter(t > burnin)

  if (smearing && smearing_method == "local") {
    beta.local <- run_local_paired_mc_periods(
      data_deagg = data.deagg,
      dgp = dgp,
      reps = smearing_reps,
      trt.IDs = trt.IDs,
      treated_rate = "beta_est",
      untreated_rate = "beta_fit",
      inf_mean = inf_mean,
      delta = delta,
      inf_var = inf_var,
      agg = agg,
      unit_population = unit_population,
      incidence_scale = incidence_scale,
      incidence_aggregation = incidence_aggregation
    )
    beta.trt <- beta.local %>% transmute(unit, week, inc = modeled_treated)
    beta.untrt <- beta.local %>% transmute(unit, week, inc = modeled_untreated)
  } else if (smearing) {
    beta.paired <- lapply(trt.IDs, function(ind) {
      pop.ind <- get_unit_population(ind, unit_population)
      I0 <- data.deagg$I[data.deagg$unit==ind & data.deagg$t==burnin+1]
      recovered <- data.deagg$R[data.deagg$unit==ind & data.deagg$t==burnin+1]
      trans_prob.trt <- data.deagg$beta_est[data.deagg$unit==ind]
      trans_prob.untrt <- data.deagg$beta_fit[data.deagg$unit==ind]
      E0 <- if (dgp == "SEIR") {
        data.deagg$E[data.deagg$unit==ind & data.deagg$t==burnin+1]
      } else NULL

      run_paired_mc_trajectories(
        dgp=dgp, reps=smearing_reps, time_steps=(T0+T1),
        pop.size=pop.ind, I0=I0, E0=E0, recovered=recovered,
        trans_prob_treated=trans_prob.trt,
        trans_prob_untreated=trans_prob.untrt,
        inf_mean=inf_mean, inf_var=inf_var, delta=delta
      )
    })
    beta.trt <- rbindlist(lapply(beta.paired, `[[`, "treated"))
    beta.untrt <- rbindlist(lapply(beta.paired, `[[`, "untreated"))
  } else {
    beta.untrt <- rbindlist(lapply(trt.IDs, function(ind) {
      pop.ind <- get_unit_population(ind, unit_population)
      I0 <- data.deagg$I[data.deagg$unit==ind & data.deagg$t==burnin+1]
      recovered <- data.deagg$R[data.deagg$unit==ind & data.deagg$t==burnin+1]
      trans_prob <- data.deagg$beta_fit[data.deagg$unit==ind]
      if (dgp == "SIR") {
        run_SIR_varying(pop.size=pop.ind, seeds=I0, recovered=recovered,
                        trans_prob=trans_prob, time_steps=(T0+T1),
                        inf_mean=inf_mean, inf_var=inf_var)
      } else if (dgp == "SEIR") {
        E0 <- data.deagg$E[data.deagg$unit==ind & data.deagg$t==burnin+1]
        run_SEIR_varying(pop.size=pop.ind, I0=I0, E0=E0, recovered=recovered,
                         trans_prob=trans_prob, time_steps=(T0+T1),
                         inf_mean=inf_mean, delta=delta)
      }
    }))

    beta.trt <- rbindlist(lapply(trt.IDs, function(ind) {
      pop.ind <- get_unit_population(ind, unit_population)
      I0 <- data.deagg$I[data.deagg$unit==ind & data.deagg$t==burnin+1]
      recovered <- data.deagg$R[data.deagg$unit==ind & data.deagg$t==burnin+1]
      trans_prob <- data.deagg$beta_est[data.deagg$unit==ind]
      if (dgp == "SIR") {
        run_SIR_varying(pop.size=pop.ind, seeds=I0, recovered=recovered,
                        trans_prob=trans_prob, time_steps=(T0+T1),
                        inf_mean=inf_mean, inf_var=inf_var)
      } else if (dgp == "SEIR") {
        E0 <- data.deagg$E[data.deagg$unit==ind & data.deagg$t==burnin+1]
        run_SEIR_varying(pop.size=pop.ind, I0=I0, E0=E0, recovered=recovered,
                         trans_prob=trans_prob, time_steps=(T0+T1),
                         inf_mean=inf_mean, delta=delta)
      }
    }))
  }

  if (!(smearing && smearing_method == "local")) {
    beta.untrt <- aggregate_simulated_incidence(
      beta.untrt, trt.IDs,
      incidence_aggregation=incidence_aggregation,
      unit_population=unit_population,
      incidence_scale=incidence_scale
    )
    beta.trt <- aggregate_simulated_incidence(
      beta.trt, trt.IDs,
      incidence_aggregation=incidence_aggregation,
      unit_population=unit_population,
      incidence_scale=incidence_scale
    )
  }

  data.untrt <- data.untrt %>%
    mutate(.unit_key = as.character(unit)) %>%
    left_join(beta.untrt %>%
                transmute(.unit_key=as.character(unit), week,
                          Y.untrt.beta=inc),
              by=c(".unit_key", "week")) %>%
    left_join(beta.trt %>%
                transmute(.unit_key=as.character(unit), week,
                          Y.trt.beta=inc),
              by=c(".unit_key", "week")) %>%
    dplyr::select(-.unit_key) %>%
    mutate(Y.untrt.beta.adj1=NA_real_, Y.untrt.beta.adj2=NA_real_)

  for (i in trt.IDs) {
    for (time in (T0/agg+1):((T0+T1)/agg)) {
      idx <- data.untrt$unit==i & data.untrt$week==time
      beta.var <- get_var(beta.fit, time, T0=T0/agg, id=i)
      data.untrt$Y.untrt.beta.adj1[idx] <- adjust_fitted_untreated(
        data.untrt$Y.untrt.beta[idx], beta.var, "approach1")
      data.untrt$Y.untrt.beta.adj2[idx] <- adjust_fitted_untreated(
        data.untrt$Y.untrt.beta[idx], beta.var, "approach2")
    }
  }

  observed <- if (sim) data.untrt$Y.trt.beta else data.untrt$inc
  smear.out <- NULL
  if (smearing) {
    effect.log <- if (is.null(coef)) as.numeric(tail(coef(beta.fit), 1)) else as.numeric(coef)
    smear.out <- if (smearing_method == "local") {
      smear_counterfactual_local(
        observed = observed,
        modeled_treated = data.untrt$Y.trt.beta,
        modeled_untreated = data.untrt$Y.untrt.beta,
        unit = data.untrt$unit,
        week = data.untrt$week,
        post_index = data.untrt$trt.time,
        log_effect = effect.log
      )
    } else {
      smear_counterfactual_trajectory(
        observed = observed,
        modeled_treated = data.untrt$Y.trt.beta,
        modeled_untreated = data.untrt$Y.untrt.beta,
        unit = data.untrt$unit,
        post_index = data.untrt$trt.time,
        log_effect = effect.log
      )
    }
  }

  build_ame_output(
    model = "beta",
    effect = exp(tail(coef(beta.fit), 1)),
    p = stata.out$p,
    observed = observed,
    untreated = data.untrt$Y.untrt.beta,
    untreated_adj1 = data.untrt$Y.untrt.beta.adj1,
    untreated_adj2 = data.untrt$Y.untrt.beta.adj2,
    untreated_smear = if (is.null(smear.out)) NULL else smear.out$untreated,
    smear_factor = if (is.null(smear.out)) NA_real_ else mean(smear.out$detail$smear.factor),
    smear_ratio = if (is.null(smear.out)) NA_real_ else smearing_weighted_ratio(smear.out$detail),
    post_index = data.untrt$trt.time
  )
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
           diff = ame_observed_minus_untreated(y, y.ctl)) %>% # observed minus untreated
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
    mutate(diff = ame_observed_minus_untreated(PosPer1K, Pos.ctl)) %>%
    group_by(ID) %>% summarise(diff = sum(diff)) # difference to get marginal effect for each unit
  mean(df$diff) # average over all units to get the AME
}
