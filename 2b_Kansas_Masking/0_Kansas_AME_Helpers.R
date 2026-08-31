# Reconstruct Kansas SEIR states on the observed incidence scale.

kansas_unit_id <- function(x) {
  if (is.factor(x)) return(as.character(x))
  if (inherits(x, "haven_labelled")) x <- unclass(x)
  as.character(x)
}

reconstruct_kansas_case_states <- function(data, units, start_date, end_date,
                                            inf_mean, delta, incidence_scale=100000) {
  required <- c(
    "ncounty", "date", "coestpop2019", "sus_frac", "stnnewcases7davg"
  )
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("Kansas SEIR state reconstruction requires column(s): ",
         paste(missing, collapse=", "), ".")
  }
  if (!is.numeric(inf_mean) || length(inf_mean)!=1L || !is.finite(inf_mean) ||
      inf_mean<=0) stop("inf_mean must be one finite positive number.")
  if (!is.numeric(delta) || length(delta)!=1L || !is.finite(delta) ||
      delta<=0) stop("delta must be one finite positive number.")

  units <- as.character(units)
  state.data <- data %>%
    mutate(unit=kansas_unit_id(ncounty)) %>%
    filter(unit %in% units) %>%
    group_by(unit) %>%
    arrange(date, .by_group=TRUE) %>%
    mutate(
      .case_count=pmax(as.numeric(stnnewcases7davg), 0) *
        as.numeric(coestpop2019) / incidence_scale,
      # E_t is delta times incidence at t+1.
      E=delta * lead(.case_count),
      I={
        z <- numeric(n())
        for (k in seq_len(n())) {
          z[k] <- (if (k==1L) 0 else (1-1/inf_mean)*z[k-1L]) + .case_count[k]
        }
        z
      },
      S=pmin(pmax(as.numeric(sus_frac) * as.numeric(coestpop2019), 0),
             as.numeric(coestpop2019)),
      R=as.numeric(coestpop2019) - S - E - I,
      S_frac=S/as.numeric(coestpop2019)
    ) %>%
    ungroup() %>%
    filter(date >= as.Date(start_date), date <= as.Date(end_date)) %>%
    group_by(unit) %>%
    arrange(date, .by_group=TRUE) %>%
    mutate(t=row_number()) %>%
    ungroup()

  counts <- table(state.data$unit)
  if (length(counts)!=length(units) || any(names(counts)!=sort(units)) ||
      length(unique(as.integer(counts)))!=1L) {
    stop("Each treated Kansas county must have the same complete daily AME window.")
  }
  state.columns <- c("S", "E", "I", "R")
  if (any(!is.finite(as.matrix(state.data[state.columns]))) ||
      any(as.matrix(state.data[state.columns]) < 0)) {
    stop("Reconstructed Kansas S, E, I, and R states must be finite and nonnegative.")
  }
  state.total <- rowSums(state.data[state.columns])
  tolerance <- pmax(1, state.data$coestpop2019) * 1e-8
  if (any(abs(state.total-state.data$coestpop2019) > tolerance)) {
    stop("Reconstructed Kansas SEIR states do not sum to the county population.")
  }

  state.data %>% dplyr::select(unit, date, t, S, E, I, R, S_frac, coestpop2019)
}

kansas_ame_monte_carlo <- function(data, state.data, fit_column,
                                   transmission_column, coefficients,
                                   trt.IDs, unit_population, inf_mean, delta,
                                   T0_days, agg=7, incidence_scale=100000,
                                   n_rep=50000, seed=2026) {
  required <- c("unit", "week", "inc", "trt.unit", "trt.time", "trt_post",
                fit_column, transmission_column)
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop("Kansas AME Monte Carlo input is missing column(s): ",
         paste(missing, collapse=", "), ".")
  }
  coefficient.names <- names(coefficients)
  coefficients <- as.numeric(coefficients)
  if (is.null(coefficient.names)) coefficient.names <- seq_along(coefficients)
  n_rep <- as.integer(n_rep)
  if (length(n_rep)!=1L || is.na(n_rep) || n_rep<2L) {
    stop("n_rep must be an integer of at least two.")
  }
  T0_days <- as.integer(T0_days)
  agg <- as.integer(agg)
  if (T0_days%%agg!=0L) stop("T0_days must be divisible by agg.")
  T0.week <- T0_days/agg

  fit.formula <- stats::as.formula(paste0(
    fit_column, " ~ -1 + factor(week) + factor(unit) + factor(trt_post)"))
  model.fit <- stats::glm(fit.formula, family=stats::poisson(), data=data)
  post.data <- data %>%
    filter(trt.unit, trt.time, as.character(unit) %in% as.character(trt.IDs)) %>%
    mutate(.unit_key=as.character(unit)) %>%
    arrange(.unit_key, week)
  post.weeks <- sort(unique(post.data$week))
  n.post <- length(post.weeks)
  if (!n.post || nrow(post.data)!=length(trt.IDs)*n.post) {
    stop("Kansas AME requires one balanced post-intervention panel.")
  }

  results <- vector("list", length(coefficients))
  for (coef.index in seq_along(coefficients)) {
    set.seed(as.integer(seed)+coef.index, kind="L'Ecuyer-CMRG")
    raw.sum <- adj1.sum <- adj2.sum <- numeric(n_rep)

    for (unit.id in as.character(trt.IDs)) {
      initial <- state.data[
        as.character(state.data$unit)==unit.id & state.data$t==T0_days,
        , drop=FALSE]
      unit.post <- post.data[post.data$.unit_key==unit.id, , drop=FALSE]
      pop.ind <- unname(unit_population[unit.id])
      if (nrow(initial)!=1L || nrow(unit.post)!=n.post ||
          length(pop.ind)!=1L || !is.finite(pop.ind) || pop.ind<=0) {
        stop("Invalid Kansas AME initial state or population for unit ", unit.id, ".")
      }

      beta.path <- rep(as.numeric(unit.post[[transmission_column]]) /
          exp(coefficients[coef.index]), each=agg)
      S <- rep(as.numeric(initial$S), n_rep)
      E <- rep(as.numeric(initial$E), n_rep)
      I <- rep(as.numeric(initial$I), n_rep)
      weekly.incidence <- matrix(0, nrow=n_rep, ncol=n.post)

      for (tt in seq_along(beta.path)) {
        incidence <- E/delta
        week.index <- ceiling(tt/agg)
        weekly.incidence[, week.index] <-
          weekly.incidence[, week.index] +
          incidence/agg/pop.ind*incidence_scale
        lambda <- beta.path[tt]*I*S/pop.ind
        if (any(!is.finite(lambda)) || any(lambda<0)) {
          stop("Invalid SEIR transmission mean for unit ", unit.id, ".")
        }
        new.exposures <- stats::rpois(n_rep, lambda)
        S <- S-new.exposures
        E <- (1-1/delta)*E+new.exposures
        I <- (1-1/inf_mean)*I+incidence
      }

      variance <- vapply(post.weeks,
        function(week.i) get_var(model.fit, week.i, T0=T0.week, id=unit.id),
        numeric(1))
      divisor1 <- 1+variance/2
      divisor2 <- exp(variance/2)
      observed <- matrix(as.numeric(unit.post$inc), nrow=n_rep, ncol=n.post, byrow=TRUE)
      raw.sum <- raw.sum+rowSums(observed-weekly.incidence)
      adj1.sum <- adj1.sum+rowSums(observed-sweep(weekly.incidence, 2, divisor1, "/"))
      adj2.sum <- adj2.sum+rowSums(observed-sweep(weekly.incidence, 2, divisor2, "/"))
    }

    denominator <- length(trt.IDs)*n.post
    raw.draws <- raw.sum/denominator
    adj1.draws <- adj1.sum/denominator
    adj2.draws <- adj2.sum/denominator
    results[[coef.index]] <- data.frame(
      parameter=coefficient.names[coef.index],
      coef=coefficients[coef.index], effect=exp(coefficients[coef.index]),
      AME=mean(raw.draws), AME_mcse=sd(raw.draws)/sqrt(n_rep),
      AME.adj1=mean(adj1.draws),
      AME.adj1_mcse=sd(adj1.draws)/sqrt(n_rep),
      AME.adj2=mean(adj2.draws),
      AME.adj2_mcse=sd(adj2.draws)/sqrt(n_rep),
      row.names=NULL)
  }
  bind_rows(results)
}
