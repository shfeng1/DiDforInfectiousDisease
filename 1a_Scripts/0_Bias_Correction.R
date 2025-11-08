get_var <- function(mod, # fitted Poisson regression model
                    time, # some post-intervention time > T0
                    T0, id) { 
  # extract covariance matrix
  dispersion <- sum(residuals(mod, type = "pearson")^2) / df.residual(mod)
  vcov <- vcov(mod, dispersion = dispersion)
  
  # unit-level variance
  var.alpha <- diag(vcov)[paste0("factor(unit)", id)]
  
  # time-level variance
  var.gamma <- sum( diag(vcov)[paste0("factor(week)", (T0+1):time)] )
  
  # covariance between unit and each time period
  cov1 <- sum( vcov[paste0("factor(unit)", id), paste0("factor(week)", (T0+1):time)] )
  
  if (time >= (T0+2)) {
    # covariance matrix subsetted to between current time and all previous time periods
    vcov.time <- vcov[(T0+1):time, (T0+1):time]
    
    # add up off-diagonal terms for the covariances between all gammas
    cov2 <- sum( vcov.time[col(vcov.time) != row(vcov.time)] )/2
  } else { # if time == T0+1, then there is no post-intervention time covariance yet
    cov2 <- 0
  }
  
  # total variances
  var <- (time-T0)^2*var.alpha + var.gamma + 2*(time-T0)*cov1 + 2*cov2
  
  return(as.numeric(var))
}
