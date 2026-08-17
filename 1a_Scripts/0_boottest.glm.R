boottest.glm <- function(fit, M=1000, gweight, model, seed=2026) {
  coef <- fit$coefficients
  score <- fit$scores
  H <- fit$hessian
  V <- fit$cov.scaled
  
  set.seed(seed, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
  coef.boot <- foreach(s = 1:M, 
                     .combine = "cbind",
                     .errorhandling = "remove") %dopar% 
    { 
      df.boot <- df_sunab %>%
        group_by(ID) %>%
        mutate(err_sgn = sample(c(-1, 1), 1)) %>%
        ungroup()

      score.boot <- score * df.boot$err_sgn
      
      as.numeric(coef + solve(H) %*% colSums(score.boot))
    }
  
  ATT_gt <- matrix(NA, nrow = length(time_to_trt), ncol = ncol(coef.boot))
  for (time in time_to_trt) {
    R <- matrix(rep(0,  length(fit$coefficients)), nrow = 1)
    
    if (model %in% c("inc", "loginc")) {
      if (time == -27 | time == 14) {
        R[,grepl(paste0("time::", time, ":"), names(coef))] <- 1
      } else if (time == -26) { # only cohort 27 and 28 has -28 time to treatment
        R[,grepl(paste0("time::", time, ":"), names(coef))] <- gweight[2:3]/sum(gweight[2:3])
      } else if (time == 13) { # only cohort 26 and 27 has 13 time to treatment
        R[,grepl(paste0("time::", time, ":"), names(coef))] <- gweight[1:2]/sum(gweight[1:2])
      } else {
        R[,grepl(paste0("time::", time, ":"), names(coef))] <- gweight
      }
    } else if (model == "growth") {
      if (time == -26 | time == 14) {
        R[,grepl(paste0("time::", time, ":"), names(coef))] <- 1
      } else if (time == -25) { # only cohort 27 and 28 has -28 time to treatment
        R[,grepl(paste0("time::", time, ":"), names(coef))] <- gweight[2:3]/sum(gweight[2:3])
      } else if (time == 13) { # only cohort 26 and 27 has 13 time to treatment
        R[,grepl(paste0("time::", time, ":"), names(coef))] <- gweight[1:2]/sum(gweight[1:2])
      } else {
        R[,grepl(paste0("time::", time, ":"), names(coef))] <- gweight
      }
    }
    
    ATT_gt[(time_to_trt==time),] <- as.numeric(R %*% coef.boot)
  }
  
  return(ATT_gt)
}
