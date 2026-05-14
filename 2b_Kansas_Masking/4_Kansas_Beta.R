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
beta_est.coef <- tail(beta_est.fit$coefficients, 1) # -0.08883859
beta_est.out <- capture.output(stata("glm beta_est trt_post i.ncounty i.week, family(poisson) link(log)
    boottest trt_post, cluster(ncounty) reps(10000)", stata.echo = T, data.in = df.in))
# z = -1.7209; p = 0.0964; b = -0.08883859; RR = 0.9149933
####################################################################################################################################
# Get confidence interval
beta_est.p <- data.frame(b0=as.numeric(beta_est.coef), p=0.0964)
for (b0 in c(seq(-0.1904, -0.1903, 0.00001), seq(0.0178, 0.0179, 0.00001))) {
  # print(b0)
  tmp.p <- stata(paste0("glm beta_est trt_post i.ncounty i.week, family(poisson) link(log)
    boottest trt_post=", b0, ", cluster(ncounty)  reps(10000)
    gen p=r(p) in 1
    keep p
    keep if _n==1"), stata.echo = F, data.in = df.in, data.out=TRUE)
  beta_est.p <- rbind(beta_est.p, as.numeric(c(b0, tmp.p)))
}
lower.bound <- min(beta_est.p$b0[beta_est.p$p >= 0.05 & beta_est.p$b0<beta_est.coef])
upper.bound <- max(beta_est.p$b0[beta_est.p$p >= 0.05 & beta_est.p$b0>beta_est.coef])

print("------------------ LOG BETAt (PREVALENCE) MODEL ------------------")
print(paste0("Treatment effect: ", round(exp(beta_est.coef), 2), " with CI: (", 
             round(exp(lower.bound), 2), ", ", round(exp(upper.bound), 2), ")", " and p-value = ",
             strsplit(trimws(tail(beta_est.out, 1)), "     ")[[1]][2]))
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
  mutate(unit = ncounty, S_frac = sus_frac, Rt = NULL, I = prevalence, R = (immune+deaths), E = infected_est, t = 1:n())

Y.obs <- mean(data.model$inc[data.model$trt_post==1])
beta_est.AME <- data.frame(type = c("point estimate", "lower bound", "upper bound"),
                           coef = c(beta_est.coef, lower.bound, upper.bound), AME = NA)
for (coef in beta_est.AME$coef) {
  set.seed(2025, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
  sim.out <- foreach(s = 1:100, 
                     .combine = "rbind",
                     .errorhandling = "remove") %dopar% 
    { 
      run_beta(data.in=data.model, out.df, dgp="SEIR", inf_mean=inf_days, delta=delta,
             trt.IDs=county.trt, coef=coef, parallel.id=s)
    }
  beta_est.AME$AME[beta_est.AME$coef==coef] <- Y.obs - mean(sim.out$Y.untrt)
}
print(paste0("AME: ", format(round(beta_est.AME$AME[1], 1), nsmall=1), " with CI: (", 
             format(round(beta_est.AME$AME[2], 1), nsmall=1), ", ", 
             format(round(beta_est.AME$AME[3], 1), nsmall=1), ")", " and p-value = ",
             strsplit(trimws(tail(beta_est.out, 1)), "     ")[[1]][2]))
# -45.8 with CI: (-95.5, 10.8) and p-value = 0.0964