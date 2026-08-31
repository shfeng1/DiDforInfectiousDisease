rm(list=ls())
here::i_am("2b_Kansas_Masking/4_Kansas_Beta_Exposure.R")
source("./global_options.R")
source("./1a_Scripts/0_SEIR.R")
source("./2b_Kansas_Masking/0_Kansas_AME_Helpers.R", local=TRUE)
inf_days <- 5; delta <- 3; burnin <- 0; agg <- 7

df.model <- readRDS("./0_Data/Kansas.rds")
df.in <- readRDS("./0_Data/Kansas_Cleaned.rds")
df.first <- df.model %>% filter(dayssincefirstcase == 1)
county.trt <- as.character(sort(unique(df.in$ncounty[df.in$trt_post])))
beta_exposure.fit <- glm(beta_exposure ~ -1 + factor(week) + factor(ncounty) + factor(trt_post), data=df.in, family=poisson())
beta_exposure.coef <- tail(beta_exposure.fit$coefficients, 1) # -0.08564578
####################################################################################################################################
# Get confidence interval
beta_exposure.p <- data.frame(b0=as.numeric(beta_exposure.coef), p=0.0429)
for (b0 in c(seq(-0.1655, -0.1654, 0.00001), seq(-0.0030, -0.0029, 0.00001))) {
  tmp.p <- stata(paste0("glm beta_exposure trt_post i.ncounty i.week, family(poisson) link(log)
    boottest trt_post=", b0, ", cluster(ncounty)  reps(10000)
    gen p=r(p) in 1
    keep p
    keep if _n==1"), stata.echo = F, data.in = df.in, data.out=TRUE)
  beta_exposure.p <- rbind(beta_exposure.p, as.numeric(c(b0, tmp.p)))
}
lower.bound <- min(beta_exposure.p$b0[beta_exposure.p$p >= 0.05 & beta_exposure.p$b0<beta_exposure.coef])
upper.bound <- max(beta_exposure.p$b0[beta_exposure.p$p >= 0.05 & beta_exposure.p$b0>beta_exposure.coef])
beta_exposure_effect_20 <- c(
  estimate=as.numeric(exp(beta_exposure.coef)),
  lower=exp(lower.bound),
  upper=exp(upper.bound)
)
beta_exposure_p_20 <- beta_exposure.p$p[1]
####################################################################################################################################
## CONVERT TO AME
unit.population.df <- df.in %>%
  filter(as.character(ncounty) %in% county.trt) %>%
  transmute(unit=as.character(ncounty), population=coestpop2019) %>%
  distinct()
unit.population <- setNames(unit.population.df$population, unit.population.df$unit)
data.model <- df.in %>% mutate(
  unit=ncounty, transmission=beta_exposure,
  inc=stnnewcases7davg, S_frac=sus_frac,
  week=week-min(df.in$week)+1)
T0 <- length(unique(data.model$start_date[!data.model$trt.time]))*agg
T1 <- length(unique(data.model$start_date))*agg-T0
out.df <- reconstruct_kansas_case_states(
  df.model, county.trt, start_date="2020-06-05", end_date="2020-12-10",
  inf_mean=inf_days, delta=delta)

beta_exposure.AME <- kansas_ame_monte_carlo(
  data=data.model, state.data=out.df, fit_column="beta_exposure",
  transmission_column="transmission",
  coefficients=c(`point estimate`=beta_exposure.coef,
                 `lower bound`=lower.bound, `upper bound`=upper.bound),
  trt.IDs=county.trt, unit_population=unit.population,
  inf_mean=inf_days, delta=delta, T0_days=T0)
beta_exposure_AME_20 <- setNames(
  beta_exposure.AME$AME,
  c("estimate", "lower", "upper")
)
