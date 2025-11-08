rm(list=ls())
here::i_am("2b_Kansas_Masking/1_Kansas_Incidence.R")
source("./global_options.R")

df.in <- readRDS("./0_Data/Kansas_Cleaned.rds")
inc.fit <- lm(stnnewcases7davg ~ factor(ncounty) + factor(week) + factor(trt_post), data=df.in)
tail(inc.fit$coefficients, 1) # -20.39995
stata("glm stnnewcases7davg trt_post i.ncounty i.week, family(gaussian) link(identity)
    boottest trt_post, cluster(ncounty) reps(10000)", stata.echo = T, data.in = df.in)
# z = -4.2559; p = 0.0000; b = -20.39995
####################################################################################################################################
# Get confidence interval from inverting null hypothesis
inc.p <- data.frame(b0 = 0, p = 0.0000)
for (b0 in c(seq(-27.59, -27.55, 0.01), seq(-13.25, -13.15, 0.01))) {
  print(b0)
  tmp.p <- stata(paste0("glm stnnewcases7davg trt_post i.ncounty i.week, family(gaussian) link(identity)
    boottest trt_post=", b0, ", cluster(ncounty)  reps(10000)
    gen p=r(p) in 1
    keep p
    keep if _n==1"), stata.echo = F, data.in = df.in, data.out=TRUE)
  inc.p <- rbind(inc.p, as.numeric(c(b0, tmp.p)))
}

lower.bound <- min(inc.p$b0[inc.p$p >= 0.05 & inc.p$b0<tail(inc.fit$coefficients, 1)])
upper.bound <- max(inc.p$b0[inc.p$p >= 0.05 & inc.p$b0>tail(inc.fit$coefficients, 1)])
c(lower.bound, upper.bound) # (-27.56, -13.20)
####################################################################################################################################
loginc.fit <- glm(stnnewcases7davg ~ factor(ncounty) + factor(week) + factor(trt_post), data=df.in, family="poisson")
tail(loginc.fit$coef, 1) # -0.9643996 
stata("glm stnnewcases7davg trt_post i.ncounty i.week, family(poisson) link(log)
    boottest trt_post, cluster(ncounty) reps(10000)", stata.echo = T, data.in = df.in)
# z = -2.4242; p = 0.0028; b = -0.9643996; RR = 0.381212
####################################################################################################################################
# Get confidence interval
loginc.p <- data.frame(b0 = 0, p = 0.0028)
for (b0 in c(seq(-1.4396, -1.4395, 0.00001), seq(-0.3557, -0.3556, 0.00001))) {
  print(b0)
  tmp.p <- stata(paste0("glm stnnewcases7davg trt_post i.ncounty i.week, family(poisson) link(log)
    boottest trt_post=", b0, ", cluster(ncounty)  reps(10000)
    gen p=r(p) in 1
    keep p
    keep if _n==1"), stata.echo = F, data.in = df.in, data.out=TRUE)
  loginc.p <- rbind(loginc.p, as.numeric(c(b0, tmp.p)))
}
lower.bound <- min(loginc.p$b0[loginc.p$p >= 0.05 & loginc.p$b0<tail(loginc.fit$coefficients, 1)])
upper.bound <- max(loginc.p$b0[loginc.p$p >= 0.05 & loginc.p$b0>tail(loginc.fit$coefficients, 1)])
c(lower.bound, upper.bound) # (-1.43955, -0.35562)
exp(c(lower.bound, upper.bound)) # RR: (0.2370344, 0.7007389)
####################################################################################################################################
# calculate AMEs
mean(loginc.fit$fitted.values[df.in$trt_post]*exp(tail(loginc.fit$coef, 1))-loginc.fit$fitted.values[df.in$trt_post]) # -19.93338
mean(loginc.fit$fitted.values[df.in$trt_post]*exp(lower.bound)-loginc.fit$fitted.values[df.in$trt_post]) # -24.57786
mean(loginc.fit$fitted.values[df.in$trt_post]*exp(upper.bound)-loginc.fit$fitted.values[df.in$trt_post]) # -9.640277
