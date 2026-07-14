rm(list=ls())
here::i_am("2b_Kansas_Masking/3_Kansas_Rt.R")
source("./global_options.R")
source("./1a_Scripts/0_SEIR.R")
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
# Treatment effect: 0.94 with CI: (0.87, 1.02) and p-value = 0.1172
####################################################################################################################################
## CONVERT TO AME
data.model <- df.in %>% mutate(unit=ncounty, Rt_est=Rt_exposure*(agg/delta), inc=infections, S_frac=sus_frac, R_est=Rt_est, week=week - min(df.in$week)+1)
burnin <- 0; agg <- 7
T0 <- length(unique(data.model$start_date[! data.model$trt.time]))*agg
T1 <- length(unique(data.model$start_date))*agg - T0
out.df <- df.model %>% filter(date >= "2020-06-05", date < "2020-12-11",
                              (ncounty %in% county.trt), # treated units
                              ! ncounty %in% df.first$ncounty[df.first$date >= "2020-06-24"]) %>%
  group_by(ncounty) %>% arrange(time) %>%
  mutate(unit = ncounty, S_frac = sus_frac, Rt = NULL, I = prevalence, R = (immune+deaths), E = infected_est, t = 1:n())

Rt_exposure.AME <- data.frame(type = c("point estimate", "lower bound", "upper bound"),
                         coef = c(Rt_exposure.coef, lower.bound, upper.bound), AME = NA,
                         AME.adj1 = NA, AME.adj2 = NA)
for (coef in Rt_exposure.AME$coef) {
  set.seed(12345, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
  sim.out <- foreach(s = 1:100, 
                     .combine = "rbind",
                     .errorhandling = "stop") %dopar% 
    { 
      tryCatch(run_Rt(data.in=data.model, out.df, dgp="SEIR", inf_mean=inf_days, delta=delta, 
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
  Rt_exposure.AME$AME[Rt_exposure.AME$coef==coef] <- mean(sim.out$AME)
  Rt_exposure.AME$AME.adj1[Rt_exposure.AME$coef==coef] <- mean(sim.out$AME.adj1)
  Rt_exposure.AME$AME.adj2[Rt_exposure.AME$coef==coef] <- mean(sim.out$AME.adj2)
}
print(paste0("AME: ", format(round(Rt_exposure.AME$AME[1], 1), nsmall=1), " with CI: (", 
             format(round(Rt_exposure.AME$AME[2], 1), nsmall=1), ", ", 
             format(round(Rt_exposure.AME$AME[3], 1), nsmall=1), ")", " and p-value = ",
             strsplit(trimws(tail(Rt_exposure.out, 1)), "     ")[[1]][2]))
# AME: -33.3 with CI: (-73.1, 2.5) and p-value = 0.1172