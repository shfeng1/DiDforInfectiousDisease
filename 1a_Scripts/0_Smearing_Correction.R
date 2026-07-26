# Empirical smearing correction for SIR/SEIR AME conversion.
#
# When smearing = TRUE, treated and untreated incidence trajectories are
# generated in paired Monte Carlo simulations. Each pair uses identical
# initial states and the same uniform random numbers. The default local method
# re-anchors both trajectories to the observed epidemic state at the start of
# each analysis period, averages each trajectory over Monte Carlo replicates,
# and then calculates the untreated-to-treated incidence ratio.
#
# For treated unit d:
#
#   smear.factor = C.obs / mean(C.model.trt)
#   C.untrt.smear = C.obs * mean(C.model.untrt) / mean(C.model.trt)
#
# This adapts Duan's empirical retransformation intuition: estimate the
# multiplicative original-scale discrepancy from observed versus fitted treated
# incidence, then apply the same factor to the fitted untreated trajectory.
# The original Appendix G.1 correction remains available with smearing = FALSE.

run_paired_mc_trajectories <- function(
    dgp,
    reps,
    time_steps,
    pop.size,
    I0,
    recovered,
    trans_prob_treated,
    trans_prob_untreated,
    inf_mean,
    E0 = NULL,
    delta = NULL,
    inf_var = NULL
) {
  reps <- as.integer(reps)
  time_steps <- as.integer(time_steps)
  if (length(reps) != 1L || is.na(reps) || reps < 1L) {
    stop("smearing_reps must be one positive integer.")
  }
  if (length(time_steps) != 1L || is.na(time_steps) || time_steps < 1L) {
    stop("time_steps must be one positive integer.")
  }
  if (length(trans_prob_treated) < time_steps ||
      length(trans_prob_untreated) < time_steps) {
    stop("Treated and untreated transmission paths must cover all simulation time steps.")
  }

  beta.trt <- as.numeric(trans_prob_treated[seq_len(time_steps)])
  beta.untrt <- as.numeric(trans_prob_untreated[seq_len(time_steps)])
  if (anyNA(beta.trt) || anyNA(beta.untrt) ||
      any(!is.finite(beta.trt)) || any(!is.finite(beta.untrt)) ||
      any(beta.trt < 0) || any(beta.untrt < 0)) {
    stop("Transmission paths must be finite and nonnegative.")
  }

  # Vectorize over Monte Carlo replicates. The old implementation called the
  # complete SIR/SEIR simulator twice for every replicate; local smearing then
  # repeated this for every treated unit-period. Here each epidemic day is one
  # vectorized update over all replicates, with the same uniforms used for the
  # treated and untreated paths.
  uniforms <- matrix(
    pmin(pmax(runif(reps * time_steps), .Machine$double.eps),
         1 - .Machine$double.eps),
    nrow = reps,
    ncol = time_steps
  )
  treated.mean <- numeric(time_steps)
  untreated.mean <- numeric(time_steps)

  if (dgp == "SIR") {
    S.trt <- rep(pop.size - I0 - recovered, reps)
    I.trt <- rep(I0, reps)
    R.trt <- rep(recovered, reps)
    S.untrt <- S.trt
    I.untrt <- I.trt
    R.untrt <- R.trt

    if (!is.null(inf_var)) {
      shape <- inf_mean^2 / inf_var
      scale <- inf_var / inf_mean
      F.gamma <- pgamma(seq_len(time_steps + 1L), shape = shape, scale = scale)
      pmf <- F.gamma - c(0, head(F.gamma, -1L))
      inc.trt <- matrix(0, nrow = reps, ncol = time_steps + 1L)
      inc.untrt <- matrix(0, nrow = reps, ncol = time_steps + 1L)
      inc.trt[, 1L] <- I0
      inc.untrt[, 1L] <- I0
    }

    for (k in seq_len(time_steps)) {
      lambda.trt <- pmax(0, beta.trt[k] * I.trt * S.trt / pop.size)
      lambda.untrt <- pmax(0, beta.untrt[k] * I.untrt * S.untrt / pop.size)
      trans.trt <- qpois(uniforms[, k], lambda.trt)
      trans.untrt <- qpois(uniforms[, k], lambda.untrt)

      treated.mean[k] <- mean(trans.trt)
      untreated.mean[k] <- mean(trans.untrt)
      S.trt <- S.trt - trans.trt
      S.untrt <- S.untrt - trans.untrt

      if (is.null(inf_var)) {
        remove.trt <- I.trt / inf_mean
        remove.untrt <- I.untrt / inf_mean
      } else {
        t.index <- k + 1L
        lag.index <- t.index - seq_len(t.index - 1L)
        weights <- pmf[seq_len(t.index - 1L)]
        remove.trt <- as.vector(inc.trt[, lag.index, drop = FALSE] %*% weights)
        remove.untrt <- as.vector(inc.untrt[, lag.index, drop = FALSE] %*% weights)
        inc.trt[, t.index] <- trans.trt
        inc.untrt[, t.index] <- trans.untrt
      }

      I.trt <- I.trt + trans.trt - remove.trt
      R.trt <- R.trt + remove.trt
      I.untrt <- I.untrt + trans.untrt - remove.untrt
      R.untrt <- R.untrt + remove.untrt
    }
  } else if (dgp == "SEIR") {
    if (is.null(E0) || is.null(delta)) {
      stop("E0 and delta are required for paired SEIR simulations.")
    }
    S.trt <- rep(pop.size - E0 - I0 - recovered, reps)
    E.trt <- rep(E0, reps)
    I.trt <- rep(I0, reps)
    R.trt <- rep(recovered, reps)
    S.untrt <- S.trt
    E.untrt <- E.trt
    I.untrt <- I.trt
    R.untrt <- R.trt

    for (k in seq_len(time_steps)) {
      lambda.trt <- pmax(0, beta.trt[k] * I.trt * S.trt / pop.size)
      lambda.untrt <- pmax(0, beta.untrt[k] * I.untrt * S.untrt / pop.size)
      exposed.trt <- qpois(uniforms[, k], lambda.trt)
      exposed.untrt <- qpois(uniforms[, k], lambda.untrt)

      # run_SEIR_varying() defines observed incidence as E[t-1] / delta.
      inc.trt <- E.trt / delta
      inc.untrt <- E.untrt / delta
      treated.mean[k] <- mean(inc.trt)
      untreated.mean[k] <- mean(inc.untrt)

      S.trt <- S.trt - exposed.trt
      S.untrt <- S.untrt - exposed.untrt
      E.trt <- E.trt - inc.trt + exposed.trt
      E.untrt <- E.untrt - inc.untrt + exposed.untrt
      remove.trt <- I.trt / inf_mean
      remove.untrt <- I.untrt / inf_mean
      I.trt <- I.trt - remove.trt + inc.trt
      I.untrt <- I.untrt - remove.untrt + inc.untrt
      R.trt <- R.trt + remove.trt
      R.untrt <- R.untrt + remove.untrt
    }
  } else {
    stop("dgp must be either 'SIR' or 'SEIR'.")
  }

  list(
    treated = data.table(t = seq_len(time_steps), inc = treated.mean),
    untreated = data.table(t = seq_len(time_steps), inc = untreated.mean)
  )
}

smear_counterfactual_trajectory <- function(
    observed,
    modeled_treated,
    modeled_untreated,
    unit,
    post_index,
    log_effect,
    tolerance = 1e-10
) {
  n <- length(observed)
  inputs <- list(
    modeled_treated = modeled_treated,
    modeled_untreated = modeled_untreated,
    unit = unit,
    post_index = post_index
  )
  bad_length <- vapply(inputs, length, integer(1)) != n
  if (any(bad_length)) {
    stop(
      "Smearing inputs must have equal lengths; mismatched input(s): ",
      paste(names(inputs)[bad_length], collapse = ", ")
    )
  }
  if (length(log_effect) != 1L || !is.finite(log_effect)) {
    stop("log_effect must be one finite fitted treatment coefficient.")
  }

  selected <- post_index %in% TRUE
  if (!any(selected)) {
    stop("No post-intervention observations were supplied for smearing.")
  }
  if (anyNA(observed[selected]) || any(!is.finite(observed[selected])) ||
      anyNA(modeled_treated[selected]) || any(!is.finite(modeled_treated[selected])) ||
      anyNA(modeled_untreated[selected]) || any(!is.finite(modeled_untreated[selected]))) {
    stop("Observed and modeled post-intervention incidence must be finite for smearing.")
  }
  if (any(observed[selected] < 0) || any(modeled_treated[selected] < 0) ||
      any(modeled_untreated[selected] < 0)) {
    stop("Incidence cannot be negative in the smearing correction.")
  }

  unit_key <- as.character(unit)
  untreated_smear <- modeled_untreated
  detail <- vector("list", length(unique(unit_key[selected])))
  names(detail) <- unique(unit_key[selected])

  for (id in names(detail)) {
    idx <- selected & unit_key == id
    C.obs <- sum(observed[idx])
    C.model.trt <- sum(modeled_treated[idx])
    C.model.untrt <- sum(modeled_untreated[idx])

    if (!is.finite(C.model.trt) || C.model.trt <= 0) {
      stop("Modeled treated cumulative incidence must be positive for unit ", id, ".")
    }
    if (!is.finite(C.model.untrt) || C.model.untrt <= 0) {
      stop("Modeled untreated cumulative incidence must be positive for unit ", id, ".")
    }

    model_ratio.raw <- C.model.untrt / C.model.trt

    # Monotonicity is assumed for cumulative incidence and imposed directly on
    # the mechanistic ratio; it is not applied to the AME itself.
    model_ratio <- if (log_effect < -tolerance) {
      max(1, model_ratio.raw)
    } else if (log_effect > tolerance) {
      min(1, model_ratio.raw)
    } else {
      1
    }

    target_C.untrt <- C.obs * model_ratio
    untreated_shape <- modeled_untreated[idx] / C.model.untrt
    untreated_smear[idx] <- target_C.untrt * untreated_shape

    detail[[id]] <- data.frame(
      unit = id,
      C.obs = C.obs,
      C.model.trt = C.model.trt,
      C.model.untrt = C.model.untrt,
      smear.factor = C.obs / C.model.trt,
      model.ratio.raw = model_ratio.raw,
      model.ratio = model_ratio,
      C.untrt.smear = target_C.untrt,
      stringsAsFactors = FALSE
    )
  }

  list(
    untreated = untreated_smear,
    detail = do.call(rbind, detail)
  )
}


# Paired Monte Carlo smearing re-anchored at the observed epidemic state at the
# beginning of each unit-period. This prevents small differences in beta from
# accumulating through a single long SIR/SEIR recursion over the entire
# post-intervention horizon.
run_local_paired_mc_periods <- function(
    data_deagg,
    dgp,
    reps,
    trt.IDs,
    treated_rate,
    untreated_rate,
    inf_mean,
    delta = NULL,
    inf_var = NULL,
    agg = 7L,
    unit_population = NULL,
    incidence_scale = NULL,
    incidence_aggregation = c("sum", "mean")
) {
  incidence_aggregation <- match.arg(incidence_aggregation)
  agg <- as.integer(agg)
  if (length(agg) != 1L || is.na(agg) || agg < 1L) {
    stop("agg must be one positive integer.")
  }
  required <- c("unit", "week", "t", "trt.time", "I", "R",
                treated_rate, untreated_rate)
  if (dgp == "SEIR") required <- c(required, "E")
  missing <- setdiff(required, names(data_deagg))
  if (length(missing)) {
    stop("Missing local-smearing input column(s): ", paste(missing, collapse = ", "))
  }

  out <- lapply(trt.IDs, function(ind) {
    unit_key <- as.character(ind)
    unit_data <- data_deagg[as.character(data_deagg$unit) == unit_key, , drop = FALSE]
    unit_data <- unit_data[order(unit_data$t), , drop = FALSE]
    post_weeks <- sort(unique(unit_data$week[unit_data$trt.time %in% TRUE]))
    pop.ind <- get_unit_population(ind, unit_population)

    do.call(rbind, lapply(post_weeks, function(w) {
      period_data <- unit_data[unit_data$week == w & unit_data$trt.time %in% TRUE, , drop = FALSE]
      period_data <- period_data[order(period_data$t), , drop = FALSE]
      if (!nrow(period_data)) stop("No post-intervention days found for unit ", unit_key, ", week ", w, ".")

      start_t <- min(period_data$t)
      state <- unit_data[unit_data$t == start_t - 1L, , drop = FALSE]
      if (nrow(state) != 1L) {
        stop("Exactly one observed state immediately before unit ", unit_key,
             ", week ", w, " is required for local smearing.")
      }

      treated_path <- period_data[[treated_rate]]
      untreated_path <- period_data[[untreated_rate]]
      if (anyNA(treated_path) || any(!is.finite(treated_path)) ||
          anyNA(untreated_path) || any(!is.finite(untreated_path)) ||
          any(treated_path < 0) || any(untreated_path < 0)) {
        stop("Invalid transmission path for unit ", unit_key, ", week ", w, ".")
      }

      paired <- run_paired_mc_trajectories(
        dgp = dgp,
        reps = reps,
        time_steps = nrow(period_data),
        pop.size = pop.ind,
        I0 = state$I,
        E0 = if (dgp == "SEIR") state$E else NULL,
        recovered = state$R,
        trans_prob_treated = treated_path,
        trans_prob_untreated = untreated_path,
        inf_mean = inf_mean,
        inf_var = inf_var,
        delta = delta
      )

      aggregate_one <- function(x) {
        value <- if (incidence_aggregation == "sum") sum(x) else mean(x)
        if (!is.null(incidence_scale)) value <- value / pop.ind * incidence_scale
        value
      }

      data.frame(
        unit = unit_key,
        week = w,
        modeled_treated = aggregate_one(paired$treated$inc),
        modeled_untreated = aggregate_one(paired$untreated$inc),
        stringsAsFactors = FALSE
      )
    }))
  })

  data.table::as.data.table(do.call(rbind, out))
}

# Apply a unit-period-specific modeled untreated-to-treated ratio to the
# observed treated outcome. Monotonicity is maintained as part of the model:
# a negative log effect implies an untreated-to-treated ratio of at least one,
# and a positive log effect implies a ratio of at most one.
smear_counterfactual_local <- function(
    observed,
    modeled_treated,
    modeled_untreated,
    unit,
    week,
    post_index,
    log_effect,
    tolerance = 1e-10
) {
  n <- length(observed)
  inputs <- list(
    modeled_treated = modeled_treated,
    modeled_untreated = modeled_untreated,
    unit = unit,
    week = week,
    post_index = post_index
  )
  bad_length <- vapply(inputs, length, integer(1)) != n
  if (any(bad_length)) {
    stop("Local smearing inputs must have equal lengths; mismatched input(s): ",
         paste(names(inputs)[bad_length], collapse = ", "))
  }
  if (length(log_effect) != 1L || !is.finite(log_effect)) {
    stop("log_effect must be one finite fitted treatment coefficient.")
  }

  selected <- post_index %in% TRUE
  if (!any(selected)) stop("No post-intervention observations were supplied for local smearing.")
  if (anyNA(observed[selected]) || any(!is.finite(observed[selected])) ||
      anyNA(modeled_treated[selected]) || any(!is.finite(modeled_treated[selected])) ||
      anyNA(modeled_untreated[selected]) || any(!is.finite(modeled_untreated[selected]))) {
    stop("Observed and modeled post-intervention outcomes must be finite for local smearing.")
  }
  if (any(observed[selected] < 0) || any(modeled_treated[selected] <= 0) ||
      any(modeled_untreated[selected] <= 0)) {
    stop("Observed incidence must be nonnegative and modeled incidence must be positive for local smearing.")
  }

  key <- paste(as.character(unit[selected]), week[selected], sep = "::")
  if (anyDuplicated(key)) stop("Local smearing requires one observation per unit-period.")

  ratio_raw <- modeled_untreated[selected] / modeled_treated[selected]
  ratio <- if (log_effect < -tolerance) {
    pmax(1, ratio_raw)
  } else if (log_effect > tolerance) {
    pmin(1, ratio_raw)
  } else {
    rep(1, length(ratio_raw))
  }

  untreated_smear <- observed
  untreated_smear[selected] <- observed[selected] * ratio

  detail <- data.frame(
    unit = as.character(unit[selected]),
    week = week[selected],
    C.obs = observed[selected],
    C.model.trt = modeled_treated[selected],
    C.model.untrt = modeled_untreated[selected],
    model.ratio.raw = ratio_raw,
    model.ratio = ratio,
    C.untrt.smear = untreated_smear[selected],
    stringsAsFactors = FALSE
  )
  detail$smear.factor <- detail$C.obs / detail$C.model.trt

  list(untreated = untreated_smear, detail = detail)
}

smearing_weighted_ratio <- function(detail) {
  if (is.null(detail) || !nrow(detail)) return(NA_real_)
  denominator <- sum(detail$C.obs)
  if (!is.finite(denominator) || denominator <= 0) return(mean(detail$model.ratio))
  sum(detail$C.obs * detail$model.ratio) / denominator
}
