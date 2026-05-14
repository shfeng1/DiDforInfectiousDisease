rm(list=ls())
here::i_am("2b_Kansas_Masking/4_Kansas_Beta_Exposure.R")
source("./global_options.R")
source("./1a_Scripts/0_SEIR.R")
inf_days <- 5; delta <- 3

df.model <- readRDS("./0_Data/Kansas.rds")
df.in <- readRDS("./0_Data/Kansas_Cleaned.rds")
df.first <- df.model %>% filter(dayssincefirstcase == 1)
county.trt <- as.character(sort(unique(df.in$ncounty[df.in$trt_post])))
beta_exposure.fit <- glm(beta_exposure ~ -1 + factor(week) + factor(ncounty) + factor(trt_post), data=df.in, family=poisson())
beta_exposure.coef <- tail(beta_exposure.fit$coefficients, 1) # -0.07247426
beta_exposure.out <- capture.output(stata("glm beta_exposure trt_post i.ncounty i.week, family(poisson) link(log)
    boottest trt_post, cluster(ncounty) reps(10000)", stata.echo = T, data.in = df.in))
# z = -2.1214; p = 0.0396; b = -0.0724743; RR = 0.9300897
####################################################################################################################################
# Get confidence interval
beta_exposure.p <- data.frame(b0=as.numeric(beta_exposure.coef), p=0.0396)
for (b0 in c(seq(-0.1390, -0.1389, 0.00001), seq(-0.0038, -0.0037, 0.00001))) {
  # print(b0)
  tmp.p <- stata(paste0("glm beta_exposure trt_post i.ncounty i.week, family(poisson) link(log)
    boottest trt_post=", b0, ", cluster(ncounty)  reps(10000)
    gen p=r(p) in 1
    keep p
    keep if _n==1"), stata.echo = F, data.in = df.in, data.out=TRUE)
  beta_exposure.p <- rbind(beta_exposure.p, as.numeric(c(b0, tmp.p)))
}
lower.bound <- min(beta_exposure.p$b0[beta_exposure.p$p >= 0.05 & beta_exposure.p$b0<beta_exposure.coef])
upper.bound <- max(beta_exposure.p$b0[beta_exposure.p$p >= 0.05 & beta_exposure.p$b0>beta_exposure.coef])

print("------------------ LOG BETAt (INSTANTANEOUS ON EXPOSURE) MODEL ------------------")
print(paste0("Treatment effect: ", round(exp(beta_exposure.coef), 2), " with CI: (", 
             format(round(exp(lower.bound), 2), nsmall=2), ", ", 
             format(round(exp(upper.bound), 2), nsmall=2), ")", " and p-value = ",
             strsplit(trimws(tail(beta_exposure.out, 1)), "     ")[[1]][2]))
# Treatment effect: 0.93 with CI: (0.87, 1.00) and p-value = 0.0396
####################################################################################################################################
## CONVERT TO AME
data.model <- df.in %>% mutate(unit=ncounty, inc=infections, S_frac=sus_frac, week=week - min(df.in$week)+1)
burnin <- 0; agg <- 7
T0 <- length(unique(data.model$start_date[! data.model$trt.time]))*agg
T1 <- length(unique(data.model$start_date))*agg - T0
out.df <- df.model %>% filter(date >= "2020-06-05", date <= "2020-12-11",
                              (ncounty %in% county.trt), # treated units
                              ! ncounty %in% df.first$ncounty[df.first$date >= "2020-06-24"]) %>%
  group_by(ncounty) %>% arrange(time) %>%
  mutate(unit = ncounty, S_frac = sus_frac, Rt = NULL, I = prevalence, R = (immune+deaths), E = infected_est, t = 1:n())

Y.obs <- mean(data.model$inc[data.model$trt_post==1])
beta_exposure.AME <- data.frame(type = c("point estimate", "lower bound", "upper bound"),
                                coef = c(beta_exposure.coef, lower.bound, upper.bound), AME = NA)
for (coef in beta_exposure.AME$coef) {
  set.seed(2025, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
  sim.out <- foreach(s = 1:100, 
                     .combine = "rbind",
                     .errorhandling = "remove") %dopar% 
    { 
      run_beta(data.in=data.model, out.df, dgp="SEIR", inf_mean=inf_days, delta=delta,
             trt.IDs=county.trt, coef=coef, parallel.id=s)
    }
  beta_exposure.AME$AME[beta_exposure.AME$coef==coef] <- Y.obs - mean(sim.out$Y.untrt)
}
print(paste0("AME: ", format(round(beta_exposure.AME$AME[1], 1), nsmall=1), " with CI: (", 
             format(round(beta_exposure.AME$AME[2], 1), nsmall=1), ", ", 
             format(round(beta_exposure.AME$AME[3], 1), nsmall=1), ")", " and p-value = ",
             strsplit(trimws(tail(beta_exposure.out, 1)), "     ")[[1]][2]))
# AME: -38.9 with CI: (-70.5, -0.5) and p-value = 0.0396