rm(list=ls())
here::i_am("2b_Kansas_Masking/4_Kansas_Beta_Exposure.R")
source("./global_options.R")
source("./1a_Scripts/0_SEIR.R")
inf_days <- 5; delta <- 3; burnin <- 0; agg <- 7

df.model <- readRDS("./0_Data/Kansas.rds")
df.in <- readRDS("./0_Data/Kansas_Cleaned.rds")
df.first <- df.model %>% filter(dayssincefirstcase == 1)
county.trt <- as.character(sort(unique(df.in$ncounty[df.in$trt_post])))
beta_exposure.fit <- glm(beta_exposure ~ -1 + factor(week) + factor(ncounty) + factor(trt_post), data=df.in, family=poisson())
beta_exposure.coef <- tail(beta_exposure.fit$coefficients, 1) # -0.08564578
beta_exposure.out <- capture.output(stata("glm beta_exposure trt_post i.ncounty i.week, family(poisson) link(log)
    boottest trt_post, cluster(ncounty) reps(10000)", stata.echo = T, data.in = df.in))
####################################################################################################################################
# Get confidence interval
beta_exposure.p <- data.frame(b0=as.numeric(beta_exposure.coef), p=0.0429)
for (b0 in c(seq(-0.1655, -0.1654, 0.00001), seq(-0.0030, -0.0029, 0.00001))) {
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
# Treatment effect: 0.92 with CI: (0.85, 1.00) and p-value = 0.0429
####################################################################################################################################
## CONVERT TO AME
data.model <- df.in %>% mutate(unit=ncounty, beta_est=beta_exposure*(agg/delta), inc=infections, S_frac=sus_frac, week=week-min(df.in$week)+1)
T0 <- length(unique(data.model$start_date[! data.model$trt.time]))*agg
T1 <- length(unique(data.model$start_date))*agg - T0
out.df <- df.model %>% filter(date >= "2020-06-05", date < "2020-12-11",
                              (ncounty %in% county.trt), # treated units
                              ! ncounty %in% df.first$ncounty[df.first$date >= "2020-06-24"]) %>%
  group_by(ncounty) %>% arrange(time) %>%
  mutate(unit = ncounty, S_frac = sus_frac, Rt = NULL, I = I_est, R = 0, E = infections, t = 1:n())

beta_exposure.AME <- data.frame(type = c("point estimate", "lower bound", "upper bound"),
                         coef = c(beta_exposure.coef, lower.bound, upper.bound), AME = NA,
                         AME.adj1 = NA, AME.adj2 = NA)
for (coef in beta_exposure.AME$coef) {
  set.seed(12345, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
  sim.out <- foreach(s = 1:100, 
                     .combine = "rbind",
                     .errorhandling = "stop") %dopar% 
    { 
      tryCatch(run_beta(data.in=data.model, out.df, dgp="SEIR", inf_mean=inf_days, delta=delta,
                        trt.IDs=county.trt, coef=coef, parallel.id=s, sim=TRUE),
               
               error = function(e) {
                 msg <- conditionMessage(e)
                 
                 # Allow the known transient NFS write failure
                 if (grepl("unable to open file for writing", msg, fixed = TRUE) &&
                     grepl("Stale NFS file handle", msg, fixed = TRUE)) {
                   return(NULL)
                 }
                 
                 # Any other error stops the parallel batch.
                 stop(e)
               })
    }
  beta_exposure.AME$AME[beta_exposure.AME$coef==coef] <- mean(sim.out$AME)
  beta_exposure.AME$AME.adj1[beta_exposure.AME$coef==coef] <- mean(sim.out$AME.adj1)
  beta_exposure.AME$AME.adj2[beta_exposure.AME$coef==coef] <- mean(sim.out$AME.adj2)
}
print(paste0("AME: ", format(round(beta_exposure.AME$AME.adj2[1], 1), nsmall=1), " with CI: (", 
             format(round(beta_exposure.AME$AME.adj2[2], 1), nsmall=1), ", ", 
             format(round(beta_exposure.AME$AME.adj2[3], 1), nsmall=1), ")", " and p-value = ",
             strsplit(trimws(tail(beta_exposure.out, 1)), "     ")[[1]][2]))
# AME: -33.8 with CI: (-77.0, -0.3) and p-value = 0.0429