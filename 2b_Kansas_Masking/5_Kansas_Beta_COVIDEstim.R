rm(list=ls())
here::i_am("2b_Kansas_Masking/5_Kansas_Beta_COVIDEstim.R")
source("./global_options.R")
source("./1a_Scripts/0_SEIR.R")
source("./2b_Kansas_Masking/0_Kansas_AME_Helpers.R", local=TRUE)
inf_days <- 5; delta <- 3; burnin <- 0; agg <- 7

df.model <- readRDS("./0_Data/Kansas.rds")
df.in <- readRDS("./0_Data/Kansas_Cleaned.rds")
df.first <- df.model %>% filter(dayssincefirstcase == 1)
county.trt <- as.character(sort(unique(df.in$ncounty[df.in$trt_post])))
beta.fit <- glm(beta ~ -1 + factor(week) + factor(ncounty) + factor(trt_post), data=df.in, family=poisson())
beta.coef <- tail(beta.fit$coefficients, 1) # -0.0627032
####################################################################################################################################
# Get confidence interval
beta.p <- data.frame(b0=as.numeric(beta.coef), p=0.1123)
for (b0 in c(seq(-0.1390, -0.1389, 0.00001), seq(0.0162, 0.0163, 0.00001))) {
  tmp.p <- stata(paste0("glm beta trt_post i.ncounty i.week, family(poisson) link(log)
    boottest trt_post=", b0, ", cluster(ncounty)  reps(10000)
    gen p=r(p) in 1
    keep p
    keep if _n==1"), stata.echo = F, data.in = df.in, data.out=TRUE)
  beta.p <- rbind(beta.p, as.numeric(c(b0, tmp.p)))
}
lower.bound <- min(beta.p$b0[beta.p$p >= 0.05 & beta.p$b0<beta.coef])
upper.bound <- max(beta.p$b0[beta.p$p >= 0.05 & beta.p$b0>beta.coef])
beta_COVIDEstim_effect_20 <- c(
  estimate=as.numeric(exp(beta.coef)),
  lower=exp(lower.bound),
  upper=exp(upper.bound)
)
beta_COVIDEstim_p_20 <- beta.p$p[1]
####################################################################################################################################
## CONVERT TO AME
unit.population.df <- df.in %>%
  filter(as.character(ncounty) %in% county.trt) %>%
  transmute(unit=as.character(ncounty), population=coestpop2019) %>%
  distinct()
unit.population <- setNames(unit.population.df$population, unit.population.df$unit)
data.model <- df.in %>% mutate(unit=ncounty, beta_covidestim=beta/inf_days,
  transmission=beta_covidestim, inc=stnnewcases7davg, S_frac=sus_frac,
  week=week-min(df.in$week)+1)
T0 <- length(unique(data.model$start_date[!data.model$trt.time]))*agg
T1 <- length(unique(data.model$start_date))*agg-T0
out.df <- reconstruct_kansas_case_states(
  df.model, county.trt, start_date="2020-06-05", end_date="2020-12-10",
  inf_mean=inf_days, delta=delta)

beta.AME <- kansas_ame_monte_carlo(
  data=data.model, state.data=out.df, fit_column="beta_covidestim",
  transmission_column="transmission",
  coefficients=c(`point estimate`=beta.coef,
                 `lower bound`=lower.bound, `upper bound`=upper.bound),
  trt.IDs=county.trt, unit_population=unit.population,
  inf_mean=inf_days, delta=delta, T0_days=T0)
beta_COVIDEstim_AME_20 <- setNames(
  beta.AME$AME,
  c("estimate", "lower", "upper")
)
