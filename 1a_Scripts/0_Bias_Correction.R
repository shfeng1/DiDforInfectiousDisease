get_var <- function(mod, # fitted Poisson regression model
                    time, # some post-intervention time > T0
                    T0, id) {
  # Appendix G.1 defines the cumulative variance of
  #   (time - T0) * unit_FE + sum(post-period time_FE).
  # Every covariance-matrix extraction must therefore use the exact coefficient
  # names for the requested unit and post-period week effects.
  dispersion <- sum(residuals(mod, type = "pearson")^2) / df.residual(mod)
  vcov.mat <- vcov(mod, dispersion = dispersion)

  unit.name <- paste0("factor(unit)", as.character(id))
  time.names <- paste0("factor(week)", (T0 + 1):time)

  # All requested post-period week effects must be estimable.
  missing.time.names <- setdiff(time.names, rownames(vcov.mat))
  if (length(missing.time.names) > 0L) {
    stop(
      "Cannot compute cumulative AME variance; missing time coefficient(s): ",
      paste(missing.time.names, collapse = ", ")
    )
  }

  # NECESSARY FIX [BC-REFERENCE-LEVEL]:
  # The model formula contains factor(unit), so the model frame normally stores
  # that evaluated variable under the name "factor(unit)", NOT under "unit".
  # Using mod$model$unit therefore returns NULL and falsely reports every unit
  # (including unit 1) as absent.  Read the fitted factor levels from xlevels,
  # with a model-frame fallback for compatibility with alternative model fits.
  unit.levels <- mod$xlevels[["factor(unit)"]]
  if (is.null(unit.levels)) {
    unit.column <- intersect(c("factor(unit)", "unit"), names(mod$model))
    if (length(unit.column) != 1L) {
      stop(
        "Cannot identify the fitted unit factor in the model object. ",
        "Expected xlevels[['factor(unit)']] or one model-frame column named ",
        "'factor(unit)'/'unit'."
      )
    }
    unit.levels <- levels(factor(mod$model[[unit.column]]))
  }

  id.character <- as.character(id)
  if (!id.character %in% as.character(unit.levels)) {
    stop(
      "Cannot compute cumulative AME variance; unit is genuinely absent ",
      "from the fitted model: ", id.character
    )
  }

  # Under -1 + factor(week) + factor(unit) + ..., the design matrix contains all
  # week indicators but omits one unit indicator to resolve the fixed-effect
  # linear dependence.  The omitted unit is the reference unit.  Its unit effect
  # is fixed at zero in this parameterization, so Var(unit_FE) and its covariance
  # with week effects are both zero.  Absence of that ONE coefficient is expected,
  # not a fitting error.  The week coefficients still carry the fitted trajectory
  # for the reference unit.
  if (unit.name %in% rownames(vcov.mat)) {
    var.alpha <- vcov.mat[unit.name, unit.name]
    cov1 <- sum(vcov.mat[unit.name, time.names, drop = TRUE])
  } else {
    reference.unit <- as.character(unit.levels)[1L]
    if (!identical(id.character, reference.unit)) {
      stop(
        "Cannot compute cumulative AME variance; coefficient for non-reference ",
        "unit ", id.character, " is missing from the fitted covariance matrix."
      )
    }
    var.alpha <- 0
    cov1 <- 0
  }

  # Sum_j Var(tau_j).
  var.gamma <- sum(diag(vcov.mat)[time.names])

  if (length(time.names) >= 2L) {
    vcov.time <- vcov.mat[time.names, time.names, drop = FALSE]
    # Sum once over each distinct pair i < j; the formula multiplies by 2 below.
    cov2 <- sum(vcov.time[upper.tri(vcov.time)])
  } else {
    cov2 <- 0
  }

  # Appendix G.1 cumulative variance:
  # (t-T0)^2 Var(eta_d) + sum Var(tau_j)
  # + 2(t-T0) sum Cov(eta_d,tau_j) + 2 sum_{i<j} Cov(tau_i,tau_j).
  variance <- (time - T0)^2 * var.alpha + var.gamma +
    2 * (time - T0) * cov1 + 2 * cov2

  if (!is.finite(variance)) {
    stop(
      "Computed cumulative AME variance is not finite for unit ",
      id.character, " at week ", time, "."
    )
  }

  # A materially negative value indicates an invalid covariance calculation.
  # Tiny negative values can arise from floating-point roundoff and are clamped.
  if (variance < -1e-10) {
    stop(
      "Computed cumulative AME variance is negative for unit ",
      id.character, " at week ", time, ": ", variance
    )
  }

  as.numeric(max(variance, 0))
}
