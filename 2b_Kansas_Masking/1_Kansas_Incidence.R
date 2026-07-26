rm(list=ls())
here::i_am("2b_Kansas_Masking/1_Kansas_Incidence.R")
source("./global_options.R")
source("./1a_Scripts/0_RStata.R")

df.in <- readRDS("./0_Data/Kansas_Cleaned.rds")
inc.fit <- lm(stnnewcases7davg ~ factor(ncounty) + factor(week) + factor(trt_post), data=df.in)
inc.coef <- tail(inc.fit$coefficients, 1) # -20.39995
# Return the wild-bootstrap p-value directly instead of scraping printed output.
inc.model.command <-
  "glm stnnewcases7davg trt_post i.ncounty i.week, family(gaussian) link(identity)"
inc.p.value <- run_stata_boottest_p(
  model_command = inc.model.command,
  hypothesis = "trt_post",
  cluster = "ncounty",
  data.in = df.in,
  reps = 10000L,
  id = "kansas_incidence_point"
)
####################################################################################################################################
# Get confidence interval from inverting null hypothesis
inc.p <- data.frame(b0 = as.numeric(inc.coef), p = inc.p.value)
for (b0 in c(seq(-27.57, -27.56, 0.001), seq(-13.21, -13.20, 0.01))) {
  # print(b0)
  tmp.p <- run_stata_boottest_p(
    model_command = inc.model.command,
    hypothesis = paste0("trt_post=", format_stata_number(b0)),
    cluster = "ncounty",
    data.in = df.in,
    reps = 10000L,
    id = "kansas_incidence_ci"
  )
  inc.p <- rbind(inc.p, data.frame(b0 = b0, p = tmp.p))
}

lower.bound <- min(inc.p$b0[inc.p$p >= 0.05 & inc.p$b0<inc.coef])
upper.bound <- max(inc.p$b0[inc.p$p >= 0.05 & inc.p$b0>inc.coef])

# Calculate the incidence AME from the fitted untreated counterfactual rather
# than assigning the treatment coefficient to the AME by definition.
inc.post <- df.in$trt_post
inc.untrt.data <- df.in
inc.untrt.data$trt_post <- FALSE
inc.Y.obs <- mean(df.in$stnnewcases7davg[inc.post])
inc.Y.untrt <- mean(predict(inc.fit, newdata=inc.untrt.data)[inc.post])
inc.AME <- inc.Y.obs - inc.Y.untrt

print("------------------ INCIDENCE MODEL ------------------")
print(paste0("Treatment effect: ", round(inc.coef, 1), " with CI: (",
             round(lower.bound, 1), ", ", round(upper.bound, 1), ")",
             "; AME (Y.obs - Y.untrt): ", round(inc.AME, 1),
             " and p-value = ",
             format(inc.p.value, digits = 4)))
# Treatment effect: -20.4 with CI: (-27.6, -13.2); AME (Y.obs - Y.untrt): -20.4 and p-value = 0
####################################################################################################################################
loginc.fit <- glm(stnnewcases7davg ~ factor(ncounty) + factor(week) + factor(trt_post), data=df.in, family="poisson")
loginc.coef <- tail(loginc.fit$coef, 1) # -0.9643996 
loginc.model.command <- "glm stnnewcases7davg trt_post i.ncounty i.week, family(poisson) link(log)"
loginc.p.value <- run_stata_boottest_p(
  model_command = loginc.model.command,
  hypothesis = "trt_post",
  cluster = "ncounty",
  data.in = df.in,
  reps = 10000L,
  id = "kansas_logincidence_point"
)
####################################################################################################################################
# Get confidence interval
loginc.p <- data.frame(b0 = as.numeric(loginc.coef), p = loginc.p.value)
for (b0 in c(seq(-1.4396, -1.4395, 0.00001), seq(-0.3557, -0.3556, 0.00001))) {
  # print(b0)
  tmp.p <- run_stata_boottest_p(
    model_command = loginc.model.command,
    hypothesis = paste0("trt_post=", format_stata_number(b0)),
    cluster = "ncounty",
    data.in = df.in,
    reps = 10000L,
    id = "kansas_logincidence_ci"
  )
  loginc.p <- rbind(loginc.p, data.frame(b0 = b0, p = tmp.p))
}
lower.bound <- min(loginc.p$b0[loginc.p$p >= 0.05 & loginc.p$b0<loginc.coef])
upper.bound <- max(loginc.p$b0[loginc.p$p >= 0.05 & loginc.p$b0>loginc.coef])

print("------------------ LOG INCIDENCE MODEL ------------------")
print(paste0("Treatment effect: ", round(exp(loginc.coef), 2), " with CI: (", 
             round(exp(lower.bound), 2), ", ", round(exp(upper.bound), 2), ")", " and p-value = ",
             format(loginc.p.value, digits = 4)))
# Treatment effect: 0.38 with CI: (0.24, 0.7) and p-value = 0.0028
####################################################################################################################################
# Calculate each AME only after constructing the corresponding untreated counterfactual. 
# No recursive log-normal correction applies to log incidence, so Y.untrt.adj1 and Y.untrt.adj2 equal Y.untrt.
loginc.AME <- data.frame(
  type = c("point estimate", "lower bound", "upper bound"),
  coef = c(tail(loginc.fit$coef, 1), lower.bound, upper.bound),
  Y.obs = mean(df.in$stnnewcases7davg[df.in$trt_post])
) %>%
  mutate(Y.untrt = Y.obs / exp(coef),
         Y.untrt.adj1 = Y.untrt,
         Y.untrt.adj2 = Y.untrt,
         AME = Y.obs - Y.untrt,
         AME.adj1 = Y.obs - Y.untrt.adj1,
         AME.adj2 = Y.obs - Y.untrt.adj2)
print(paste0("AME: ", round(loginc.AME$AME[1], 1), " with CI: (",
             round(loginc.AME$AME[2], 1), ", ", round(loginc.AME$AME[3], 1), ")",
             " and p-value = ", format(loginc.p.value, digits = 4)))
# AME: -52.3 with CI: (-103.7, -13.8) and p-value = 0.0028