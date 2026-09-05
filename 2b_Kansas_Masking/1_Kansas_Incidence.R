rm(list=ls())
here::i_am("2b_Kansas_Masking/1_Kansas_Incidence.R")
source("./global_options.R")

df.in <- readRDS("./0_Data/Kansas_Cleaned.rds")
inc.fit <- lm(stnnewcases7davg ~ factor(ncounty) + factor(week) + factor(trt_post), data=df.in)
inc.coef <- tail(inc.fit$coefficients, 1) # -20.39995
####################################################################################################################################
# Get confidence interval from inverting null hypothesis
inc.p <- data.frame(b0 = as.numeric(inc.coef), p = 0.0000)
for (b0 in c(seq(-27.57, -27.56, 0.001), seq(-13.20, -13.19, 0.001))) {
  tmp.p <- stata(paste0("glm stnnewcases7davg trt_post i.ncounty i.week, family(gaussian) link(identity)
    boottest trt_post=", b0, ", cluster(ncounty)  reps(10000)
    gen p=r(p) in 1
    keep p
    keep if _n==1"), stata.echo = F, data.in = df.in, data.out=TRUE)
  inc.p <- rbind(inc.p, as.numeric(c(b0, tmp.p)))
}

lower.bound <- min(inc.p$b0[inc.p$p >= 0.05 & inc.p$b0<inc.coef])
upper.bound <- max(inc.p$b0[inc.p$p >= 0.05 & inc.p$b0>inc.coef])
inc_effect_20 <- c(estimate=as.numeric(inc.coef),
                   lower=lower.bound,
                   upper=upper.bound)
inc_AME_20 <- inc_effect_20
inc_p_20 <- inc.p$p[1]
####################################################################################################################################
loginc.fit <- glm(stnnewcases7davg ~ factor(ncounty) + factor(week) + factor(trt_post), data=df.in, family="poisson")
loginc.coef <- tail(loginc.fit$coef, 1) # -0.9643996
####################################################################################################################################
# Get confidence interval
loginc.p <- data.frame(b0 = as.numeric(loginc.coef), p = 0.0028)
for (b0 in c(seq(-1.4396, -1.4395, 0.00001), seq(-0.3557, -0.3556, 0.00001))) {
  tmp.p <- stata(paste0("glm stnnewcases7davg trt_post i.ncounty i.week, family(poisson) link(log)
    boottest trt_post=", b0, ", cluster(ncounty)  reps(10000)
    gen p=r(p) in 1
    keep p
    keep if _n==1"), stata.echo = F, data.in = df.in, data.out=TRUE)
  loginc.p <- rbind(loginc.p, as.numeric(c(b0, tmp.p)))
}
lower.bound <- min(loginc.p$b0[loginc.p$p >= 0.05 & loginc.p$b0<loginc.coef])
upper.bound <- max(loginc.p$b0[loginc.p$p >= 0.05 & loginc.p$b0>loginc.coef])
loginc_effect_20 <- c(estimate=as.numeric(exp(loginc.coef)),
                      lower=exp(lower.bound),
                      upper=exp(upper.bound))
loginc_p_20 <- loginc.p$p[1]
####################################################################################################################################
# calculate AMEs
loginc.AME <- mean(df.in$stnnewcases7davg[df.in$trt_post] - df.in$stnnewcases7davg[df.in$trt_post]/exp(tail(loginc.fit$coef, 1)))
loginc.AME.lo <- mean(df.in$stnnewcases7davg[df.in$trt_post] - df.in$stnnewcases7davg[df.in$trt_post]/exp(lower.bound))
loginc.AME.hi <- mean(df.in$stnnewcases7davg[df.in$trt_post] - df.in$stnnewcases7davg[df.in$trt_post]/exp(upper.bound))
loginc_AME_20 <- c(estimate=loginc.AME,
                   lower=loginc.AME.lo,
                   upper=loginc.AME.hi)
