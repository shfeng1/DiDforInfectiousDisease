rm(list=ls())
here::i_am("2b_Kansas_Masking/5_Kansas_Beta_COVIDEstim.R")
source("./global_options.R")
source("./1a_Scripts/0_SEIR.R")
inf_days <- 5; delta <- 3

df.model <- readRDS("./0_Data/Kansas.rds")
df.in <- readRDS("./0_Data/Kansas_Cleaned.rds")
df.first <- df.model %>% filter(dayssincefirstcase == 1)
county.trt <- as.character(sort(unique(df.in$ncounty[df.in$trt_post])))
beta.fit <- glm(beta ~ -1 + factor(week) + factor(ncounty) + factor(trt_post), data=df.in, family=poisson())
beta.coef <- tail(beta.fit$coefficients, 1) # -0.0627032
beta.out <- capture.output(stata("glm beta trt_post i.ncounty i.week, family(poisson) link(log)
    boottest trt_post, cluster(ncounty) reps(10000)", stata.echo = T, data.in = df.in))
# z = -1.6445; p = 0.1123; b = -0.0627032; RR = 0.9392222
####################################################################################################################################
# Get confidence interval
beta.p <- data.frame(b0=as.numeric(beta.coef), p=0.1123)
for (b0 in c(seq(-0.1390, -0.1389, 0.00001), seq(0.0162, 0.0163, 0.00001))) {
  # print(b0)
  tmp.p <- stata(paste0("glm beta trt_post i.ncounty i.week, family(poisson) link(log)
    boottest trt_post=", b0, ", cluster(ncounty)  reps(10000)
    gen p=r(p) in 1
    keep p
    keep if _n==1"), stata.echo = F, data.in = df.in, data.out=TRUE)
  beta.p <- rbind(beta.p, as.numeric(c(b0, tmp.p)))
}
lower.bound <- min(beta.p$b0[beta.p$p >= 0.05 & beta.p$b0<beta.coef])
upper.bound <- max(beta.p$b0[beta.p$p >= 0.05 & beta.p$b0>beta.coef])

print("------------------ LOG BETAt (COVIDestim) MODEL ------------------")
print(paste0("Treatment effect: ", round(exp(beta.coef), 2), " with CI: (", 
             round(exp(lower.bound), 2), ", ", round(exp(upper.bound), 2), ")", " and p-value = ",
             strsplit(trimws(tail(beta.out, 1)), "     ")[[1]][2]))
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
beta.AME <- data.frame(type = c("point estimate", "lower bound", "upper bound"),
                       coef = c(beta.coef, lower.bound, upper.bound), AME = NA)
for (coef in beta.AME$coef) {
  set.seed(2025, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
  sim.out <- foreach(s = 1:100,
                     .combine = "rbind",
                     .errorhandling = "remove") %dopar% 
    { 
      run_beta(data.in=data.model, out.df, dgp="SEIR", inf_mean=inf_days, delta=delta,
             trt.IDs=county.trt, coef=coef, parallel.id=s)
    }
  beta.AME$AME[beta.AME$coef==coef] <- Y.obs - mean(sim.out$Y.untrt)
}
print(paste0("AME: ", format(round(beta.AME$AME[1], 1), nsmall=1), " with CI: (", 
             format(round(beta.AME$AME[2], 1), nsmall=1), ", ", 
             format(round(beta.AME$AME[3], 1), nsmall=1), ")", " and p-value = ",
             strsplit(trimws(tail(beta.out, 1)), "     ")[[1]][2]))
# AME: -30.7 with CI: (-73.0, 8.2) and p-value = 0.1123
