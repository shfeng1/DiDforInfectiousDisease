rm(list=ls())
here::i_am("2b_Kansas_Masking/3_Kansas_Rt_Exposure.R")
source("./global_options.R")
source("./1a_Scripts/0_SEIR.R")
source("./2b_Kansas_Masking/0_Kansas_AME_Helpers.R")
inf_days <- 5; delta <- 3

df.model <- readRDS("./0_Data/Kansas.rds")
df.in <- readRDS("./0_Data/Kansas_Cleaned.rds")
df.first <- df.model %>% filter(dayssincefirstcase == 1)
county.trt <- as.character(sort(unique(df.in$ncounty[df.in$trt_post])))
Rt_exposure.fit <- glm(Rt_exposure ~ -1 + factor(week) + factor(ncounty) + factor(trt_post), data=df.in, family=poisson())
Rt_exposure.coef <- tail(Rt_exposure.fit$coefficients, 1) # -0.06194844
Rt_exposure.out <- capture.output(stata("glm Rt_exposure trt_post i.ncounty i.week, family(poisson) link(log)
    boottest trt_post, cluster(ncounty) reps(10000)", stata.echo = T, data.in = df.in))
####################################################################################################################################
# Get confidence interval
Rt_exposure.p <- data.frame(b0=as.numeric(Rt_exposure.coef), p=0.1172)
for (b0 in c(seq(-0.1377, -0.1376, 0.00001), seq(0.0171, 0.0172, 0.00001))) {
  # print(b0)
  tmp.p <- stata(paste0("glm Rt_exposure trt_post i.ncounty i.week, family(poisson) link(log)
    boottest trt_post=", b0, ", cluster(ncounty)  reps(10000)
    gen p=r(p) in 1
    keep p
    keep if _n==1"), stata.echo = F, data.in = df.in, data.out=TRUE)
  Rt_exposure.p <- rbind(Rt_exposure.p, as.numeric(c(b0, tmp.p)))
}
lower.bound <- min(Rt_exposure.p$b0[Rt_exposure.p$p >= 0.05 & Rt_exposure.p$b0<Rt_exposure.coef])
upper.bound <- max(Rt_exposure.p$b0[Rt_exposure.p$p >= 0.05 & Rt_exposure.p$b0>Rt_exposure.coef])

print("------------------ LOG Rt (INSTANTANEOUS ON EXPOSURE) MODEL ------------------")
print(paste0("Treatment effect: ", round(exp(Rt_exposure.coef), 2), " with CI: (", 
             round(exp(lower.bound), 2), ", ", round(exp(upper.bound), 2), ")", " and p-value = ",
             strsplit(trimws(tail(Rt_exposure.out, 1)), "     ")[[1]][2]))
####################################################################################################################################
## CONVERT TO AME
unit.population.df <- df.in %>%
  filter(as.character(ncounty) %in% county.trt) %>%
  transmute(unit=as.character(ncounty), population=coestpop2019) %>%
  distinct()
unit.population <- setNames(unit.population.df$population, unit.population.df$unit)
data.model <- df.in %>% mutate(
  unit=ncounty, R_est=Rt_exposure, inc=stnnewcases7davg,
  S_frac=sus_frac, transmission=Rt_exposure/sus_frac,
  week=week-min(df.in$week)+1)
burnin <- 0; agg <- 7
T0 <- length(unique(data.model$start_date[!data.model$trt.time]))*agg
T1 <- length(unique(data.model$start_date))*agg-T0
out.df <- reconstruct_kansas_case_states(
  df.model, county.trt, start_date="2020-06-05", end_date="2020-12-10",
  inf_mean=inf_days, delta=delta)

Rt_exposure.AME <- kansas_ame_monte_carlo(
  data=data.model, state.data=out.df, fit_column="R_est",
  transmission_column="transmission",
  coefficients=c(`point estimate`=Rt_exposure.coef,
                 `lower bound`=lower.bound, `upper bound`=upper.bound),
  trt.IDs=county.trt, unit_population=unit.population,
  inf_mean=inf_days, delta=delta, T0_days=T0)
print(paste0("AME: ", format(round(Rt_exposure.AME$AME[1],1), nsmall=1), " with CI: (",
             format(round(Rt_exposure.AME$AME[2],1), nsmall=1), ", ",
             format(round(Rt_exposure.AME$AME[3],1), nsmall=1), ")"))
