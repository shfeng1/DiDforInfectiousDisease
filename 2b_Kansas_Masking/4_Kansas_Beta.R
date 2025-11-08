rm(list=ls())
here::i_am("2b_Kansas_Masking/4_Kansas_Beta.R")
source("./global_options.R")
source("./1a_Scripts/0_SEIR.R")
inf_days <- 5; delta <- 3

df.model <- readRDS("./0_Data/Kansas.rds")
df.in <- readRDS("./0_Data/Kansas_Cleaned.rds")
df.first <- df.model %>% filter(dayssincefirstcase == 1)
county.trt <- as.character(sort(unique(df.in$ncounty[df.in$trt_post])))

beta_est.fit <- glm(beta_est ~ -1 + factor(week) + factor(ncounty) + factor(trt_post), data=df.in, family=poisson())
tail(beta_est.fit$coefficients, 1) # -0.08883859

stata("glm beta_est trt_post i.ncounty i.week, family(poisson) link(log)
    boottest trt_post, cluster(ncounty) reps(10000)", stata.echo = T, data.in = df.in)
# z = -1.7209; p = 0.0964; b = -0.08883859; RR = 0.9149933
####################################################################################################################################
# Get confidence interval
beta_est.p <- data.frame(b0=0, p=0.0996)
for (b0 in c(seq(-0.1904, -0.1903, 0.00001), seq(0.01780, 0.01790, 0.00001))) {
  print(b0)
  tmp.p <- stata(paste0("glm beta_est trt_post i.ncounty i.week, family(poisson) link(log)
    boottest trt_post=", b0, ", cluster(ncounty)  reps(10000)
    gen p=r(p) in 1
    keep p
    keep if _n==1"), stata.echo = F, data.in = df.in, data.out=TRUE)
  beta_est.p <- rbind(beta_est.p, as.numeric(c(b0, tmp.p)))
}

lower.bound <- min(beta_est.p$b0[beta_est.p$p >= 0.05 & beta_est.p$b0<tail(beta_est.fit$coefficients, 1)])
upper.bound <- max(beta_est.p$b0[beta_est.p$p >= 0.05 & beta_est.p$b0>tail(beta_est.fit$coefficients, 1)])
c(lower.bound, upper.bound) # (-0.19037, 0.01789)
exp(c(lower.bound, upper.bound)) # (0.8266532, 1.0180510)
####################################################################################################################################
## CONVERT TO AME
data.model <- df.in %>% mutate(unit=ncounty, inc=infections, S_frac=sus_frac, week=week - min(df.in$week)+1)
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
beta.out <- data.frame(type = c("point estimate", "lower bound", "upper bound"),
                         coef = c(tail(beta_est.fit$coefficients, 1), lower.bound, upper.bound),
                         AME.fit = NA, AME.adj1 = NA, AME.adj2 = NA)
for (coef in beta.out$coef) {
  set.seed(12345, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
  sim.out <- foreach(s = 1:100, 
                     .combine = "rbind",
                     .errorhandling = "remove") %dopar% 
    { 
      run_beta(data.in=data.model, out.df, dgp="SEIR", inf_mean=inf_days, delta=delta,
             trt.IDs=county.trt, coef=coef, parallel.id=s)
    }
  beta.out$AME.fit[beta.out$coef==coef] <- Y.obs - mean(sim.out$Y.untrt)
  beta.out$AME.adj1[beta.out$coef==coef] <- Y.obs - mean(sim.out$Y.untrt.adj1)
  beta.out$AME.adj2[beta.out$coef==coef] <- Y.obs - mean(sim.out$Y.untrt.adj2)
}

beta.out
#             type        coef   AME.fit  AME.adj1  AME.adj2
# 1 point estimate -0.08883859 -28.80563 -17.19021 -16.26130
# 2    lower bound -0.19037000 -80.33327 -66.78184 -65.86729
# 3    upper bound  0.01789000  21.82733  29.45702  30.13702
