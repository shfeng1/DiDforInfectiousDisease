school_loginc_ame <- function(coef, event.times, data, subset=NULL) {
  ind <- if (is.null(subset)) event.times >= 0 else event.times %in% subset
  ATT <- mean(coef[ind])

  out <- data %>%
    mutate(event_time=time-group) %>%
    filter(if (is.null(subset)) event_time >= 0 else event_time %in% subset) %>%
    mutate(diff=y-y/exp(ATT)) %>%
    group_by(event_time) %>%
    summarise(diff=mean(diff), .groups="drop")

  sum(out$diff)
}

school_growth_variance_model <- function(data) {
  variance.data <- data %>%
    transmute(growth=y, week=time, unit=ID,
      trt_post=group < 10000 & time >= group, wt=wt)

  glm(growth ~ -1 + factor(week) + factor(unit) + factor(trt_post), 
      family=poisson, weights=wt, data=variance.data)
}

school_growth_ame <- function(coef, event.times, data, variance.model, subset=NULL) {
  coef <- as.matrix(coef)
  ind <- if (is.null(subset)) event.times >= 0 else event.times %in% subset
  ATT <- colMeans(coef[ind, , drop=FALSE])

  paths <- data %>%
    mutate(event_time=time-group) %>%
    filter(group < 10000, event_time >= -1) %>%
    arrange(ID, time)
  units <- unique(paths$ID)
  totals <- matrix(0, nrow=3, ncol=ncol(coef), dimnames=list(c("AME", "AME.adj1", "AME.adj2"), NULL))

  for (unit in units) {
    unit.data <- paths[paths$ID == unit, ]
    baseline <- unit.data[unit.data$event_time == -1, ]
    post <- unit.data[unit.data$event_time >= 0, ]
    complete.path <- nrow(post) > 0L && all(diff(post$time) == 1)

    untreated <- rep(baseline$PosPer1K, ncol(coef))
    for (i in seq_len(nrow(post))) {
      untreated <- untreated*post$y[i]/exp(ATT)
      if (is.null(subset) || post$event_time[i] %in% subset) {
        variance <- get_var(
          variance.model, post$time[i], T0=post$time[i]-1, id=unit
        )
        untreated.adj1 <- adjust_fitted_untreated(
          untreated, variance, "approach1"
        )
        untreated.adj2 <- adjust_fitted_untreated(
          untreated, variance, "approach2"
        )
        totals["AME", ] <- totals["AME", ]+post$PosPer1K[i]-untreated
        totals["AME.adj1", ] <- totals["AME.adj1", ]+post$PosPer1K[i]-untreated.adj1
        totals["AME.adj2", ] <- totals["AME.adj2", ]+post$PosPer1K[i]-untreated.adj2
      }
      if (!is.null(subset) && post$event_time[i] >= max(subset)) break
    }
  }

  totals/length(units)
}
