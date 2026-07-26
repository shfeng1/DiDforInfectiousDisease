rm(list=ls())
here::i_am("2b_Kansas_Masking/4_Kansas_Beta_Exposure.R")
source("./global_options.R")
source("./1a_Scripts/0_SEIR.R")
inf_days <- 5; delta <- 3; burnin <- 0; agg <- 7

# Local paired-Monte-Carlo smearing re-anchors the SIR/SEIR
# trajectories at the observed state at the start of every week.
# FALSE retains the original Appendix G.1 log-normal trajectory correction.
smearing <- TRUE
smearing.method <- "local"
smearing.reps <- 500L
ame.column <- if (smearing) "AME.smear" else "AME.adj2"
n.ame.sim <- if (smearing) 1L else 100L

df.model <- readRDS("./0_Data/Kansas.rds")
df.in <- readRDS("./0_Data/Kansas_Cleaned.rds")
df.first <- df.model %>% filter(dayssincefirstcase == 1)
county.trt <- as.character(sort(unique(df.in$ncounty[df.in$trt_post])))

incidence.scale <- 100000
unit.population.df <- df.in %>%
  filter(as.character(ncounty) %in% county.trt) %>%
  transmute(unit = as.character(ncounty), population = coestpop2019) %>%
  distinct()
if (anyDuplicated(unit.population.df$unit) || nrow(unit.population.df) != length(county.trt)) {
  stop("Each treated county must have exactly one population value.")
}
unit.population <- setNames(unit.population.df$population, unit.population.df$unit)
beta_exposure.fit <- glm(beta_exposure ~ -1 + factor(week) + factor(ncounty) + factor(trt_post), data=df.in, family=poisson())
beta_exposure.coef <- tail(beta_exposure.fit$coefficients, 1) # -0.08564578
beta_exposure.model.command <-
  "glm beta_exposure trt_post i.ncounty i.week, family(poisson) link(log)"
beta_exposure.p.value <- run_stata_boottest_p(
  model_command = beta_exposure.model.command,
  hypothesis = "trt_post",
  cluster = "ncounty",
  data.in = df.in,
  reps = 10000L,
  id = "kansas_beta_exposure_point"
)
####################################################################################################################################
# Get confidence interval
beta_exposure.p <- data.frame(b0 = as.numeric(beta_exposure.coef), p = beta_exposure.p.value)
for (b0 in c(seq(-0.1655, -0.1654, 0.00001), seq(-0.0030, -0.0029, 0.00001))) {
  # print(b0)
  tmp.p <- run_stata_boottest_p(
    model_command = beta_exposure.model.command,
    hypothesis = paste0("trt_post=", format_stata_number(b0)),
    cluster = "ncounty",
    data.in = df.in,
    reps = 10000L,
    id = "kansas_beta_exposure_ci"
  )
  beta_exposure.p <- rbind(beta_exposure.p, data.frame(b0 = b0, p = tmp.p))
}
lower.bound <- min(beta_exposure.p$b0[beta_exposure.p$p >= 0.05 & beta_exposure.p$b0<beta_exposure.coef])
upper.bound <- max(beta_exposure.p$b0[beta_exposure.p$p >= 0.05 & beta_exposure.p$b0>beta_exposure.coef])

print("------------------ LOG BETAt (INSTANTANEOUS ON EXPOSURE) MODEL ------------------")
print(paste0("Treatment effect: ", round(exp(beta_exposure.coef), 2), " with CI: (", 
             format(round(exp(lower.bound), 2), nsmall=2), ", ", 
             format(round(exp(upper.bound), 2), nsmall=2), ")", " and p-value = ",
             format(beta_exposure.p.value, digits = 4)))
# Treatment effect: 0.92 with CI: (0.85, 1.00) and p-value = 0.0429
####################################################################################################################################
## CONVERT TO AME
data.model <- df.in %>%
  mutate(unit = ncounty,
         beta_est = beta_exposure,
         inc = infections / coestpop2019 * incidence.scale,
         S_frac = sus_frac,
         week = week - min(df.in$week) + 1)
T0 <- length(unique(data.model$start_date[! data.model$trt.time]))*agg
T1 <- length(unique(data.model$start_date))*agg - T0
out.df <- df.model %>% filter(date >= "2020-06-05", date < "2020-12-11",
                              (ncounty %in% county.trt), # treated units
                              ! ncounty %in% df.first$ncounty[df.first$date >= "2020-06-24"]) %>%
  group_by(ncounty) %>% arrange(time) %>%
  mutate(unit = ncounty, S_frac = sus_frac, Rt = NULL,
         I = I_est,
         E = delta * lead(infections),
         R = coestpop2019 * (1 - sus_frac) - I - E,
         t = 1:n())

beta_exposure.AME <- data.frame(
  type = c("point estimate", "lower bound", "upper bound"),
  coef = c(beta_exposure.coef, lower.bound, upper.bound),
  Y.obs = NA_real_, Y.untrt = NA_real_,
  Y.untrt.adj1 = NA_real_, Y.untrt.adj2 = NA_real_,
  Y.untrt.smear = NA_real_,
  smear.factor = NA_real_, smear.ratio = NA_real_
)
for (coef in beta_exposure.AME$coef) {
  set.seed(12345, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
  sim.out <- foreach(s = 1:n.ame.sim, 
                     .combine = "rbind",
                     .errorhandling = "stop") %dopar% 
    { 
      tryCatch(run_beta(data.in=data.model, out.df, dgp="SEIR", inf_mean=inf_days, delta=delta,
                        trt.IDs=county.trt, coef=coef, parallel.id=s, sim=FALSE,
                        unit_population=unit.population,
                        incidence_scale=incidence.scale,
                        incidence_aggregation="mean",
                      smearing=smearing, smearing_reps=smearing.reps,
                      smearing_method=smearing.method,
                      p_value=beta_exposure.p.value),
               
               error = function(e) {
                 msg <- conditionMessage(e)
                 
                 stale.nfs <- grepl("unable to open file for writing", msg, fixed = TRUE) &&
                   grepl("Stale NFS file handle", msg, fixed = TRUE)
                 timed.out <- grepl("unable to open file for writing", msg, fixed = TRUE) &&
                   grepl("Operation timed out", msg, fixed = TRUE)
                 if (stale.nfs || timed.out) {
                   return(NULL)
                 }

                 stop(e)
               })
    }
  idx <- beta_exposure.AME$coef == coef
  beta_exposure.AME[idx, c("Y.obs", "Y.untrt", "Y.untrt.adj1", "Y.untrt.adj2")] <- c(
    mean(sim.out$Y.obs),
    mean(sim.out$Y.untrt),
    mean(sim.out$Y.untrt.adj1),
    mean(sim.out$Y.untrt.adj2)
  )
  if (smearing) {
    beta_exposure.AME$Y.untrt.smear[idx] <- mean(sim.out$Y.untrt.smear)
    beta_exposure.AME$smear.factor[idx] <- mean(sim.out$smear.factor)
    beta_exposure.AME$smear.ratio[idx] <- mean(sim.out$smear.ratio)
  }
}
beta_exposure.AME <- beta_exposure.AME %>%
  mutate(AME = Y.obs - Y.untrt,
         AME.adj1 = Y.obs - Y.untrt.adj1,
         AME.adj2 = Y.obs - Y.untrt.adj2,
         AME.smear = Y.obs - Y.untrt.smear) %>%
  mutate(AME.selected = .data[[ame.column]])
if (smearing) {
  print(beta_exposure.AME %>% dplyr::select(type, coef, Y.obs, smear.ratio, Y.untrt.smear, AME.selected))
}
print(paste0("AME: ", format(round(beta_exposure.AME$AME.selected[1], 1), nsmall=1), " with CI: (", 
             format(round(beta_exposure.AME$AME.selected[2], 1), nsmall=1), ", ", 
             format(round(beta_exposure.AME$AME.selected[3], 1), nsmall=1), ")", " and p-value = ",
             format(beta_exposure.p.value, digits = 4)))
# AME: -6.5 with CI: (-13.3, -0.2) and p-value = 0.0429
# Result depends on smearing: TRUE uses AME.smear; FALSE uses AME.adj2.
