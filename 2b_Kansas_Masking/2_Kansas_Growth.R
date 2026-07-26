rm(list=ls())
here::i_am("2b_Kansas_Masking/2_Kansas_Growth.R")
source("./global_options.R")
source("./1a_Scripts/0_Run_Estimators.R")

df.in <- readRDS("./0_Data/Kansas_Cleaned.rds")
county.trt <- sort(unique(df.in$ncounty[df.in$trt_post]))
growth.fit <- glm(growth ~ -1 + factor(week) + factor(ncounty) + factor(trt_post), data=df.in, family=poisson())
growth.coef <- tail(growth.fit$coefficients, 1) # -0.005761539
growth.model.command <- "glm growth trt_post i.ncounty i.week, family(poisson) link(log)"
growth.p.value <- run_stata_boottest_p(
  model_command = growth.model.command,
  hypothesis = "trt_post",
  cluster = "ncounty",
  data.in = df.in,
  reps = 10000L,
  id = "kansas_growth_point"
)
####################################################################################################################################
# Get confidence interval
growth.p <- data.frame(b0 = as.numeric(growth.coef), p = growth.p.value)
for (b0 in c(seq(-0.1815, -0.1814, 0.00001), seq(0.0199, 0.0200, 0.00001))) {
  tmp.p <- run_stata_boottest_p(
    model_command = growth.model.command,
    hypothesis = paste0("trt_post=", format_stata_number(b0)),
    cluster = "ncounty",
    data.in = df.in,
    reps = 10000L,
    id = "kansas_growth_ci"
  )
  growth.p <- rbind(growth.p, data.frame(b0 = b0, p = tmp.p))
}
lower.bound <- min(growth.p$b0[growth.p$p >= 0.05 & growth.p$b0<growth.coef])
upper.bound <- max(growth.p$b0[growth.p$p >= 0.05 & growth.p$b0>growth.coef])

print("------------------ LOG GROWTH MODEL ------------------")
print(paste0("Treatment effect: ", round(exp(growth.coef), 2), " with CI: (", 
             round(exp(lower.bound), 2), ", ", round(exp(upper.bound), 2), ")", " and p-value = ",
             format(growth.p.value, digits = 4)))
####################################################################################################################################
## CONVERT TO AME
data.model <- df.in %>%
  mutate(unit = ncounty,
         inc = infections / coestpop2019 * 100000,
         S_frac = sus_frac,
         week = week - min(df.in$week) + 1)
T0 <- length(unique(data.model$start_date[! data.model$trt.time]))
T1 <- length(unique(data.model$start_date[data.model$trt.time])); burnin <- 0; agg <- 1

Y.obs <- mean(data.model$inc[data.model$trt_post==1])
growth.AME <- data.frame(
  type = c("point estimate", "lower bound", "upper bound"),
  coef = c(tail(growth.fit$coefficients, 1), lower.bound, upper.bound),
  Y.obs = NA_real_, Y.untrt = NA_real_,
  Y.untrt.adj1 = NA_real_, Y.untrt.adj2 = NA_real_
)
for (coef in growth.AME$coef) {
  tmp <- run_growth(data.model, trt.IDs=county.trt, coef=coef, sim=FALSE)
  idx <- growth.AME$coef == coef
  growth.AME[idx, c("Y.obs", "Y.untrt", "Y.untrt.adj1", "Y.untrt.adj2")] <-
    tmp[1, c("Y.obs", "Y.untrt", "Y.untrt.adj1", "Y.untrt.adj2")]
}
growth.AME <- growth.AME %>%
  mutate(AME = Y.obs - Y.untrt,
         AME.adj1 = Y.obs - Y.untrt.adj1,
         AME.adj2 = Y.obs - Y.untrt.adj2)
print(paste0("AME: ", format(round(growth.AME$AME.adj2[1], 1), nsmall=1), " with CI: (", 
             format(round(growth.AME$AME.adj2[2], 1), nsmall=1), ", ", 
             format(round(growth.AME$AME.adj2[3], 1), nsmall=1), ")", " and p-value = ",
             format(growth.p.value, digits = 4)))
# AME: -183.3 with CI: (-1348.1, 42.1) and p-value = 0.1073