rm(list=ls())
here::i_am("2b_Kansas_Masking/3_Kansas_Rt.R")
source("./global_options.R")
source("./1a_Scripts/0_SEIR.R")
inf_days <- 5; delta <- 3

df.model <- readRDS("./0_Data/Kansas.rds")
df.in <- readRDS("./0_Data/Kansas_Cleaned.rds")
df.first <- df.model %>% filter(dayssincefirstcase == 1)
county.trt <- as.character(sort(unique(df.in$ncounty[df.in$trt_post])))

Rt_est.fit <- glm(Rt_est ~ -1 + factor(week) + factor(ncounty) + factor(trt_post), data=df.in, family=poisson())
tail(Rt_est.fit$coefficients, 1) # -0.06769345

stata("glm Rt_est trt_post i.ncounty i.week, family(poisson) link(log)
    boottest trt_post, cluster(ncounty) reps(10000)", stata.echo = T, data.in = df.in)
# z = -1.3702; p = 0.1895; b = -0.06769345; RR = 0.9345469
####################################################################################################################################
# Get confidence interval
Rt_est.p <- data.frame(b0=0, p=0.1895)
for (b0 in c(seq(-0.1660, -0.1659, 0.00001), seq(0.0365, 0.0366, 0.00001))) {
  print(b0)
  tmp.p <- stata(paste0("glm Rt_est trt_post i.ncounty i.week, family(poisson) link(log)
    boottest trt_post=", b0, ", cluster(ncounty)  reps(10000)
    gen p=r(p) in 1
    keep p
    keep if _n==1"), stata.echo = F, data.in = df.in, data.out=TRUE)
  Rt_est.p <- rbind(Rt_est.p, as.numeric(c(b0, tmp.p)))
}

lower.bound <- min(Rt_est.p$b0[Rt_est.p$p >= 0.05 & Rt_est.p$b0<tail(Rt_est.fit$coefficients, 1)])
upper.bound <- max(Rt_est.p$b0[Rt_est.p$p >= 0.05 & Rt_est.p$b0>tail(Rt_est.fit$coefficients, 1)])
c(lower.bound, upper.bound) # (-0.16597, 0.03655)
exp(c(lower.bound, upper.bound)) # RR: (0.8470716, 1.0372262)
####################################################################################################################################
## CONVERT TO AME
data.model <- df.in %>% mutate(unit=ncounty, inc=infections, S_frac=sus_frac, R_est=Rt_est, week=week - min(df.in$week)+1)
burnin <- 0; agg <- 7
T0 <- length(unique(data.model$start_date[! data.model$trt.time]))*agg
T1 <- length(unique(data.model$start_date))*agg - T0
out.df <- df.model %>% filter(date >= "2020-06-05", date < "2020-12-11",
                              (ncounty %in% county.trt), # treated units
                              ! ncounty %in% df.first$ncounty[df.first$date >= "2020-06-24"]) %>%
  group_by(ncounty) %>% arrange(time) %>%
  mutate(unit = ncounty, S_frac = sus_frac, Rt = NULL, I = infections, R = (immune+deaths), E = infected_est, t = 1:n())

# Observed Y (treated)
Y.obs <- mean(data.model$inc[data.model$trt_post==1])
Rt_est.out <- data.frame(type = c("point estimate", "lower bound", "upper bound"),
                     coef = c(tail(Rt_est.fit$coefficients, 1), lower.bound, upper.bound),
                     AME.fit = NA, AME.adj1 = NA, AME.adj2 = NA)
for (coef in Rt_est.out$coef) {
  set.seed(12345, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
  sim.out <- foreach(s = 1:100, 
                     .combine = "rbind",
                     .errorhandling = "remove") %dopar% 
    { 
      run_Rt(data.in=data.model, out.df, type="est", dgp="SEIR", inf_mean=inf_days, delta=delta, 
             trt.IDs=county.trt, coef=coef, parallel.id=s)
    }
  Rt_est.out$AME.fit[Rt_est.out$coef==coef] <- Y.obs - mean(sim.out$Y.untrt)
  Rt_est.out$AME.adj1[Rt_est.out$coef==coef] <- Y.obs - mean(sim.out$Y.untrt.adj1)
  Rt_est.out$AME.adj2[Rt_est.out$coef==coef] <- Y.obs - mean(sim.out$Y.untrt.adj2)
}
Rt_est.out
#             type        coef   AME.fit   AME.adj1   AME.adj2
# 1 point estimate -0.06769345 -16.02926  -4.801676  -3.804949
# 2    lower bound -0.16597000 -65.87541 -52.136033 -51.083223
# 3    upper bound  0.03655000  30.84589  37.912504  38.605293
