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
Rt_est.coef <- tail(Rt_est.fit$coefficients, 1) # -0.06769345
Rt_est.out <- capture.output(stata("glm Rt_est trt_post i.ncounty i.week, family(poisson) link(log)
    boottest trt_post, cluster(ncounty) reps(10000)", stata.echo = T, data.in = df.in))
# z = -1.3702; p = 0.1895; b = -0.06769345; RR = 0.9345469
####################################################################################################################################
# Get confidence interval
Rt_est.p <- data.frame(b0=as.numeric(Rt_est.coef), p=0.1895)
for (b0 in c(seq(-0.1660, -0.1659, 0.00001), seq(0.0365, 0.0366, 0.00001))) {
  # print(b0)
  tmp.p <- stata(paste0("glm Rt_est trt_post i.ncounty i.week, family(poisson) link(log)
    boottest trt_post=", b0, ", cluster(ncounty)  reps(10000)
    gen p=r(p) in 1
    keep p
    keep if _n==1"), stata.echo = F, data.in = df.in, data.out=TRUE)
  Rt_est.p <- rbind(Rt_est.p, as.numeric(c(b0, tmp.p)))
}
lower.bound <- min(Rt_est.p$b0[Rt_est.p$p >= 0.05 & Rt_est.p$b0<Rt_est.coef])
upper.bound <- max(Rt_est.p$b0[Rt_est.p$p >= 0.05 & Rt_est.p$b0>Rt_est.coef])

print("------------------ LOG Rt (PREVALENCE) MODEL ------------------")
print(paste0("Treatment effect: ", round(exp(Rt_est.coef), 2), " with CI: (", 
             round(exp(lower.bound), 2), ", ", round(exp(upper.bound), 2), ")", " and p-value = ",
             strsplit(trimws(tail(Rt_est.out, 1)), "     ")[[1]][2]))
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
  mutate(unit = ncounty, S_frac = sus_frac, Rt = NULL, I = prevalence, R = (immune+deaths), E = infected_est, t = 1:n())

Y.obs <- mean(data.model$inc[data.model$trt_post==1])
Rt_est.AME <- data.frame(type = c("point estimate", "lower bound", "upper bound"),
                         coef = c(Rt_est.coef, lower.bound, upper.bound), AME = NA)
for (coef in Rt_est.AME$coef) {
  set.seed(2025, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
  sim.out <- foreach(s = 1:100, 
                     .combine = "rbind",
                     .errorhandling = "remove") %dopar% 
    { 
      run_Rt(data.in=data.model, out.df, type="est", dgp="SEIR", inf_mean=inf_days, delta=delta, 
             trt.IDs=county.trt, coef=coef, parallel.id=s)
    }
  Rt_est.AME$AME[Rt_est.AME$coef==coef] <- Y.obs - mean(sim.out$Y.untrt)
}
print(paste0("AME: ", format(round(Rt_est.AME$AME[1], 1), nsmall=1), " with CI: (", 
             format(round(Rt_est.AME$AME[2], 1), nsmall=1), ", ", 
             format(round(Rt_est.AME$AME[3], 1), nsmall=1), ")", " and p-value = ",
             strsplit(trimws(tail(Rt_est.out, 1)), "     ")[[1]][2]))
# AME: -31.3 with CI: (-81.1, 19.1) and p-value = 0.1895
