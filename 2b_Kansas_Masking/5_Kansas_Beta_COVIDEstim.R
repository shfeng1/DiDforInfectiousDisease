rm(list=ls())
here::i_am("2b_Kansas_Masking/5_Kansas_Beta_COVIDEstim.R")
source("./global_options.R")
source("./1a_Scripts/0_SEIR.R")
inf_days <- 5; delta <- 3; burnin <- 0; agg <- 7

df.model <- readRDS("./0_Data/Kansas.rds")
df.in <- readRDS("./0_Data/Kansas_Cleaned.rds")
df.first <- df.model %>% filter(dayssincefirstcase == 1)
county.trt <- as.character(sort(unique(df.in$ncounty[df.in$trt_post])))
beta.fit <- glm(beta ~ -1 + factor(week) + factor(ncounty) + factor(trt_post), data=df.in, family=poisson())
beta.coef <- tail(beta.fit$coefficients, 1) # -0.0627032
beta.out <- capture.output(stata("glm beta trt_post i.ncounty i.week, family(poisson) link(log)
    boottest trt_post, cluster(ncounty) reps(10000)", stata.echo = T, data.in = df.in))
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
incidence.scale <- 100000
unit.population.df <- df.in %>%
  filter(as.character(ncounty) %in% county.trt) %>%
  transmute(unit=as.character(ncounty), population=coestpop2019) %>%
  distinct()
unit.population <- setNames(unit.population.df$population, unit.population.df$unit)
data.model <- df.in %>% mutate(unit=ncounty, beta_est=beta/(agg/inf_days),
  inc=stnnewcases7davg, S_frac=sus_frac,
  week=week-min(df.in$week)+1)
T0 <- length(unique(data.model$start_date[!data.model$trt.time]))*agg
T1 <- length(unique(data.model$start_date))*agg-T0
out.df <- df.model %>% filter(date >= "2020-06-05", date < "2020-12-11",
  ncounty %in% county.trt, !ncounty %in% df.first$ncounty[df.first$date >= "2020-06-24"]) %>%
  group_by(ncounty) %>% arrange(time) %>%
  mutate(unit=ncounty, S=sus_frac*coestpop2019, S_frac=sus_frac, Rt=NULL,
         I=I_est, E=infections, R=0, t=1:n())

beta.AME <- data.frame(type=c("point estimate","lower bound","upper bound"),
  coef=c(beta.coef, lower.bound, upper.bound), AME=NA, AME.adj1=NA, AME.adj2=NA)

for (k in seq_len(nrow(beta.AME))) {
  coef_i <- as.numeric(beta.AME$coef[k])
  sim.out <- foreach(s=1:50, .combine="rbind", .errorhandling="stop", .export="coef_i") %dopar% {
    set.seed(12345, kind="L'Ecuyer-CMRG")
    tryCatch(
      run_beta(data.in=data.model, out.df=out.df, dgp="SEIR", inf_mean=inf_days, delta=delta,
               trt.IDs=county.trt, coef=coef_i, parallel.id=s,
               unit_population=unit.population, incidence_scale=incidence.scale,
               simulate_from_trt=FALSE, difference=TRUE),
      error=function(e) {
        msg <- conditionMessage(e)
        allowed_write_error <- grepl("unable to open file for writing", msg, fixed=TRUE) &&
          (grepl("Stale NFS file handle", msg, fixed=TRUE) || grepl("Operation timed out", msg, fixed=TRUE))
        if (allowed_write_error) return(NULL)
        stop(e)
      }
    )
  }
  beta.AME$AME[k] <- mean(sim.out$AME)
  beta.AME$AME.adj1[k] <- mean(sim.out$AME.adj1)
  beta.AME$AME.adj2[k] <- mean(sim.out$AME.adj2)
}
print(paste0("AME: ", format(round(beta.AME$AME[1],1), nsmall=1), " with CI: (",
             format(round(beta.AME$AME[2],1), nsmall=1), ", ",
             format(round(beta.AME$AME[3],1), nsmall=1), ")"))
