# Custom Callaway & Sant'Anna-style group-time GLM estimator used for Table A2.
# Never-treated units must be coded g = 0; treated cohorts must have g > 0.
# The function intentionally returns the weighted calendar-time ATT contributions
# used in the original workflow; sum(attgt.glm(...)) gives the reported Table A2 ATT.

attgt.glm <- function(yname, tname, idname, gname, weightsname,
                      control_group = "notyettreated", data_in,
                      family) {
  if (!identical(control_group, "notyettreated")) {
    stop("attgt.glm() only implements control_group = 'notyettreated'.")
  }

  vars <- c(yname, tname, idname, gname, weightsname)
  vars <- vars[!is.na(vars) & nzchar(vars)]
  missing.vars <- setdiff(vars, names(data_in))
  if (length(missing.vars) > 0) {
    stop("Missing required variable(s): ", paste(missing.vars, collapse = ", "), ".")
  }

  y <- data_in[[yname]]
  time <- data_in[[tname]]
  ID <- data_in[[idname]]
  group <- data_in[[gname]]

  if (!is.numeric(time) || !is.numeric(group)) {
    stop("The time and treatment-cohort variables must be numeric.")
  }
  if (anyNA(time) || anyNA(group) || any(!is.finite(time)) || any(!is.finite(group))) {
    stop("The time and treatment-cohort variables cannot contain NA or non-finite values.")
  }

  glist <- sort(unique(group[group > 0]))
  if (length(glist) == 0) stop("No treated cohorts (g > 0) were found.")
  if (any(!(glist - 1) %in% time)) {
    stop("Each treated cohort must have its reference period g - 1 in the data.")
  }

  observed.time <- sort(unique(time))
  required.time <- seq(min(glist) - 1, max(time))
  if (!all(required.time %in% observed.time)) {
    stop("attgt.glm() requires consecutive integer time periods from min(g) - 1 through max(time).")
  }

  if (is.null(weightsname)) {
    wt <- rep(1, length(y))
  } else {
    wt <- data_in[[weightsname]]
    if (!is.numeric(wt) || anyNA(wt) || any(!is.finite(wt)) || any(wt < 0)) {
      stop("Weights must be finite, nonnegative numeric values.")
    }
  }

  df <- data.frame(ID = ID, group = group, time = time, y = y, wt = wt)

  ATT_gt <- matrix(
    NA_real_,
    nrow = length(glist),
    ncol = length(min(glist):max(time)),
    dimnames = list(as.character(glist), as.character(min(glist):max(time)))
  )

  for (g in glist) {
    for (t in g:max(time)) {
      df.gt <- df[
        ((df$group == g) | (df$group == 0) | (df$group > t)) &
          ((df$time == t) | (df$time == g - 1)),
        , drop = FALSE
      ]
      df.gt$trt_post <- as.integer(df.gt$group == g & df.gt$time == t)

      fit.gt <- stats::glm(
        y ~ -1 + factor(ID) + factor(time) + trt_post,
        data = df.gt,
        weights = df.gt$wt,
        family = family
      )

      beta.gt <- stats::coef(fit.gt)["trt_post"]
      if (length(beta.gt) != 1 || is.na(beta.gt)) {
        stop("trt_post is not estimable for g = ", g, ", t = ", t, ".")
      }
      ATT_gt[as.character(g), as.character(t)] <- unname(beta.gt)
    }
  }

  # Preserve the original Table A2 aggregation: cohort shares are calculated
  # among all ever-treated observations and are not renormalized at each time.
  weight <- vapply(
    glist,
    function(g) mean(df$group[df$group > 0] == g),
    numeric(1)
  )

  out <- vapply(seq_len(ncol(ATT_gt)), function(col) {
    ind <- which(!is.na(ATT_gt[, col]))
    if (length(ind) == 0) return(NA_real_)
    as.numeric(weight[ind] %*% ATT_gt[ind, col])
  }, numeric(1))

  unname(out)
}
