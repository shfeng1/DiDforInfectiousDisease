boottest.glm <- function(fit, M=5000, gweight, model, seed=12345) {
  coef <- fit$coefficients
  score <- fit$scores
  H <- fit$hessian
  V <- fit$cov.scaled
  
  set.seed(seed, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
  coef.boot <- foreach(s = 1:M,
                     .combine = "cbind",
                     # NECESSARY CHANGE [BOOT-2]: never silently discard an
                     # unknown bootstrap failure.  Removing failed replicates
                     # changes the requested bootstrap distribution and can
                     # conceal a systematic coding or numerical error.
                     .errorhandling = "stop") %dopar%
    { 
      # NECESSARY CHANGE [BOOT-1]: Appendix K defines the cluster score as
      # the SUM of observation-level score vectors within a cluster, multiplied
      # by one cluster-level mean-zero/variance-one random weight.  Dividing
      # every score contribution by n() changes that sum into a cluster average
      # and therefore implements a different bootstrap.
      cluster.sign <- df_sunab %>%
        dplyr::distinct(ID) %>%
        dplyr::mutate(err_sgn = ifelse(rbinom(dplyr::n(), 1, 0.5) == 1, -1, 1))

      score.boot <- do.call(cbind, lapply(1:ncol(score), function(i) {
        df_sunab %>%
          dplyr::mutate(score = score[, i]) %>%
          dplyr::left_join(cluster.sign, by = "ID") %>%
          dplyr::mutate(score.boot = err_sgn * score) %>%
          dplyr::pull(score.boot)
      }))
      
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
