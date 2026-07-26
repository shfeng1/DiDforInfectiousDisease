rm(list=ls())
here::i_am("2b_Kansas_Masking/3_Kansas_Rt.R")
source("./global_options.R")
source("./1a_Scripts/0_SEIR.R")
inf_days <- 5; delta <- 3

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

# NECESSARY CHANGE [E13]: the manuscript reports Kansas AMEs as the weekly mean
# of daily cases per 100,000 residents.  Simulate each county at its own Census
# population, then standardize both observed and simulated incidence to 100,000.
incidence.scale <- 100000
unit.population.df <- df.in %>%
  filter(as.character(ncounty) %in% county.trt) %>%
  transmute(unit = as.character(ncounty), population = coestpop2019) %>%
  distinct()
if (anyDuplicated(unit.population.df$unit) || nrow(unit.population.df) != length(county.trt)) {
  stop("Each treated county must have exactly one population value.")
}
unit.population <- setNames(unit.population.df$population, unit.population.df$unit)
Rt_est.fit <- glm(Rt_est ~ -1 + factor(week) + factor(ncounty) + factor(trt_post), data=df.in, family=poisson())
Rt_est.coef <- tail(Rt_est.fit$coefficients, 1) # -0.06769345
# NECESSARY FIX [STATA-BATCH/PVALUE]: return r(p) directly as data.
Rt_est.model.command <-
  "glm Rt_est trt_post i.ncounty i.week, family(poisson) link(log)"
Rt_est.p.value <- run_stata_boottest_p(
  model_command = Rt_est.model.command,
  hypothesis = "trt_post",
  cluster = "ncounty",
  data.in = df.in,
  reps = 10000L,
  id = "kansas_Rt_prevalence_point"
)
####################################################################################################################################
# Get confidence interval
Rt_est.p <- data.frame(b0 = as.numeric(Rt_est.coef), p = Rt_est.p.value)
for (b0 in c(seq(-0.1660, -0.1659, 0.00001), seq(0.0365, 0.0366, 0.00001))) {
  # print(b0)
  tmp.p <- run_stata_boottest_p(
    model_command = Rt_est.model.command,
    hypothesis = paste0("trt_post=", format_stata_number(b0)),
    cluster = "ncounty",
    data.in = df.in,
    reps = 10000L,
    id = "kansas_Rt_prevalence_ci"
  )
  Rt_est.p <- rbind(Rt_est.p, data.frame(b0 = b0, p = tmp.p))
}
lower.bound <- min(Rt_est.p$b0[Rt_est.p$p >= 0.05 & Rt_est.p$b0<Rt_est.coef])
upper.bound <- max(Rt_est.p$b0[Rt_est.p$p >= 0.05 & Rt_est.p$b0>Rt_est.coef])

print("------------------ LOG Rt (PREVALENCE) MODEL ------------------")
print(paste0("Treatment effect: ", round(exp(Rt_est.coef), 2), " with CI: (", 
             round(exp(lower.bound), 2), ", ", round(exp(upper.bound), 2), ")", " and p-value = ",
             format(Rt_est.p.value, digits = 4)))
# Treatment effect: 0.93 with CI: (0.85, 1.04) and p-value = 0.1895
####################################################################################################################################
## CONVERT TO AME
# NECESSARY CHANGE [E13]: align observed incidence with the manuscript's
# per-100,000 daily-mean outcome before comparing it with simulated incidence.
data.model <- df.in %>%
  mutate(unit = ncounty,
         inc = infections / coestpop2019 * incidence.scale,
         S_frac = sus_frac,
         R_est = Rt_est,
         week = week - min(df.in$week) + 1)
burnin <- 0; agg <- 7
T0 <- length(unique(data.model$start_date[! data.model$trt.time]))*agg
T1 <- length(unique(data.model$start_date))*agg - T0
out.df <- df.model %>% filter(date >= "2020-06-05", date < "2020-12-11",
                              (ncounty %in% county.trt), # treated units
                              ! ncounty %in% df.first$ncounty[df.first$date >= "2020-06-24"]) %>%
  group_by(ncounty) %>% arrange(time) %>%
  mutate(unit = ncounty, S_frac = sus_frac, Rt = NULL,
         I = prevalence,
         # NECESSARY CHANGE [STATE-1]: run_SEIR_varying() requires E0 to be
         # the Exposed-compartment STOCK.  Under Eq. 10, incidence[t+1] =
         # E[t]/delta, so E[t] = delta * incidence[t+1].  infected_est is the
         # FLOW of new exposures, not the E-compartment stock.
         E = delta * lead(infections),
         R = coestpop2019 * (1 - sus_frac) - I - E,
         t = 1:n())

Rt_est.AME <- data.frame(
  type = c("point estimate", "lower bound", "upper bound"),
  coef = c(Rt_est.coef, lower.bound, upper.bound),
  Y.obs = NA_real_, Y.untrt = NA_real_,
  Y.untrt.adj1 = NA_real_, Y.untrt.adj2 = NA_real_,
  Y.untrt.smear = NA_real_,
  smear.factor = NA_real_, smear.ratio = NA_real_
)
for (coef in Rt_est.AME$coef) {
  set.seed(12345, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
  sim.out <- foreach(s = 1:n.ame.sim, 
                     .combine = "rbind",
                     .errorhandling = "stop") %dopar% 
    { 
      # NECESSARY CHANGE [E11]: use observed treated incidence minus
      # the simulated untreated trajectory, not a simulated treated contrast.
      tryCatch(run_Rt(data.in=data.model, out.df, dgp="SEIR", inf_mean=inf_days, delta=delta, 
                      trt.IDs=county.trt, coef=coef, parallel.id=s, sim=FALSE,
                      unit_population=unit.population,
                      incidence_scale=incidence.scale,
                      incidence_aggregation="mean",
                      smearing=smearing, smearing_reps=smearing.reps,
                      smearing_method=smearing.method,
                      p_value=Rt_est.p.value),
               
               error = function(e) {
                 msg <- conditionMessage(e)
                 
                 # NECESSARY CHANGE [PAR-1]: ignore only the two
                 # explicitly known transient file-system failures.  Every other
                 # error still stops the parallel batch.
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
  idx <- Rt_est.AME$coef == coef
  Rt_est.AME[idx, c("Y.obs", "Y.untrt", "Y.untrt.adj1", "Y.untrt.adj2")] <- c(
    mean(sim.out$Y.obs),
    mean(sim.out$Y.untrt),
    mean(sim.out$Y.untrt.adj1),
    mean(sim.out$Y.untrt.adj2)
  )
  if (smearing) {
    Rt_est.AME$Y.untrt.smear[idx] <- mean(sim.out$Y.untrt.smear)
    Rt_est.AME$smear.factor[idx] <- mean(sim.out$smear.factor)
    Rt_est.AME$smear.ratio[idx] <- mean(sim.out$smear.ratio)
  }
}
Rt_est.AME <- Rt_est.AME %>%
  mutate(AME = Y.obs - Y.untrt,
         AME.adj1 = Y.obs - Y.untrt.adj1,
         AME.adj2 = Y.obs - Y.untrt.adj2,
         AME.smear = Y.obs - Y.untrt.smear) %>%
  mutate(AME.selected = .data[[ame.column]])
if (smearing) {
  print(Rt_est.AME %>% dplyr::select(type, coef, Y.obs, smear.ratio, Y.untrt.smear, AME.selected))
}
print(paste0("AME: ", format(round(Rt_est.AME$AME.selected[1], 1), nsmall=1), " with CI: (", 
             format(round(Rt_est.AME$AME.selected[2], 1), nsmall=1), ", ", 
             format(round(Rt_est.AME$AME.selected[3], 1), nsmall=1), ")", " and p-value = ",
             format(Rt_est.p.value, digits = 4)))
# Result depends on smearing: TRUE uses AME.smear; FALSE uses AME.adj2.
