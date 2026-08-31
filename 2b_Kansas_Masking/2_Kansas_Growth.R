rm(list=ls())
here::i_am("2b_Kansas_Masking/2_Kansas_Growth.R")
source("./global_options.R")
source("./1a_Scripts/0_Run_Estimators.R", local=TRUE)

df.in <- readRDS("./0_Data/Kansas_Cleaned.rds")
county.trt <- sort(unique(df.in$ncounty[df.in$trt_post]))
growth.fit <- glm(growth ~ -1 + factor(week) + factor(ncounty) + factor(trt_post), data=df.in, family=poisson())
growth.coef <- tail(growth.fit$coefficients, 1) # -0.08332984
####################################################################################################################################
# Get confidence interval
growth.p <- data.frame(b0=as.numeric(growth.coef), p=0.1073)
for (b0 in c(seq(-0.1815, -0.1814, 0.00001), seq(0.0199, 0.0200, 0.00001))) {
  tmp.p <- stata(paste0("glm growth trt_post i.ncounty i.week, family(poisson) link(log)
    boottest trt_post=", b0, ", cluster(ncounty)  reps(10000)
    gen p=r(p) in 1
    keep p
    keep if _n==1"), stata.echo = F, data.in = df.in, data.out=TRUE)
  growth.p <- rbind(growth.p, as.numeric(c(b0, tmp.p)))
}
lower.bound <- min(growth.p$b0[growth.p$p >= 0.05 & growth.p$b0<growth.coef])
upper.bound <- max(growth.p$b0[growth.p$p >= 0.05 & growth.p$b0>growth.coef])
growth_effect_20 <- c(
  estimate=as.numeric(exp(growth.coef)),
  lower=exp(lower.bound),
  upper=exp(upper.bound)
)
growth_p_20 <- growth.p$p[1]
####################################################################################################################################
## CONVERT TO AME
data.model <- df.in %>% mutate(unit = ncounty, inc = infections, S_frac = sus_frac, week = week - min(df.in$week)+1)
T0 <- length(unique(data.model$start_date[! data.model$trt.time]))
T1 <- length(unique(data.model$start_date[data.model$trt.time])); burnin <- 0; agg <- 1

Y.obs <- mean(data.model$inc[data.model$trt_post==1])
growth.AME <- data.frame(type = c("point estimate", "lower bound", "upper bound"),
                         coef = c(tail(growth.fit$coefficients, 1), lower.bound, upper.bound),
                         AME = NA, AME.adj1 = NA, AME.adj2 = NA)
for (coef in growth.AME$coef) {
  tmp <- run_growth(data.model, trt.IDs=county.trt, coef=coef)
  growth.AME$AME[growth.AME$coef==coef] <- mean(tmp$AME)
  growth.AME$AME.adj1[growth.AME$coef==coef] <- mean(tmp$AME.adj1)
  growth.AME$AME.adj2[growth.AME$coef==coef] <- mean(tmp$AME.adj2)
}
growth_AME_20 <- setNames(growth.AME$AME, c("estimate", "lower", "upper"))
