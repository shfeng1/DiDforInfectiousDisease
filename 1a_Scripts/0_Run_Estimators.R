# The helper functions below run each of the estimators in the following steps
# 1. Fit lm or glm model
# 2. get p-value from calling Stata using wild score cluster bootstrap
# 3. Reover AMEs
# 4. Bias correction if needed
# 5. Summarize and return results

source("./1a_Scripts/0_RStata.R")
source("./1a_Scripts/0_Bias_Correction.R")

run_inc <- function(data.in, parallel.id) {
  inc.fit <- lm(inc ~ -1 + factor(week) + factor(unit) + factor(trt_post), data = data.in)
  
  command <- "glm inc i.unit i.week i.trt_post, family(gaussian) link(identity)
    boottest 1.trt_post, cluster(unit) reps(10000) quietly
    gen p = r(p) in 1
    keep p
    keep if _n==1"
  stata.out <- my_RStata(src=command, data.in=data.in, data.out=TRUE, stata.echo=FALSE, id=parallel.id)
  
  out <- data.frame(model="inc", effect=tail(coef(inc.fit), 1), p=stata.out$p, 
                    Y.untrt=(mean(data.in$inc[data.in$trt_post])-tail(coef(inc.fit), 1)),
                    Y.untrt.adj1=NA, Y.untrt.adj2=NA)
  rownames(out) <- NULL
  return(out)
}

run_loginc <- function(data.in, parallel.id) {
  loginc.fit <- glm(inc ~ -1 + factor(week) + factor(unit) + factor(trt_post), family = poisson, data = data.in)
  
  command <- "glm inc i.unit i.week 1.trt_post, family(poisson) link(log)
    boottest 1.trt_post, cluster(unit) reps(10000) quietly
    gen p = r(p) in 1
    keep p
    keep if _n==1"
  stata.out <- my_RStata(src=command, data.in=data.in, data.out=TRUE, stata.echo=FALSE, id=parallel.id)
  
  # For log incidence models, predict based on the fitted model to directly recover the untreated Y
  data.untrt <- data.in %>% filter(trt.unit) %>% mutate(trt_post = FALSE)
  data.untrt$loginc_fit <- predict(loginc.fit, newdata = data.untrt, type = "response")
  
  out <- data.frame(model="loginc", effect=exp(tail(coef(loginc.fit), 1)), p=stata.out$p, 
                    Y.untrt=mean(data.untrt$loginc_fit[data.untrt$trt.time]),
                    Y.untrt.adj1=NA, Y.untrt.adj2=NA)
  rownames(out) <- NULL
  return(out)
}

run_growth <- function(data.in, parallel.id=0, trt.IDs=1:N1, coef=NULL) {
  growth.fit <- glm(growth ~ -1 + factor(week) + factor(unit) + factor(trt_post), family = poisson, data = data.in)
  command <- "glm growth i.unit i.week 1.trt_post, family(poisson) link(log)
    boottest 1.trt_post, cluster(unit) reps(10000) quietly
    gen p = r(p) in 1
    keep p
    keep if _n==1"
  stata.out <- my_RStata(src=command, data.in=data.in, data.out=TRUE, stata.echo=FALSE, id=parallel.id)
  
  data.untrt <- data.in %>% filter(trt.unit) %>% mutate(trt_post = FALSE)
  if (is.null(coef)) {
    # For log growth models:
    # first calculate what the untreated potential outcome would be in the scale of growth rate for each treated unit
    # then use the untreated growth rates to recover the trajectory of untreated potential outcome on the case scale
    data.untrt$growth_fit <- predict(growth.fit, newdata = data.untrt, type = "response") 
  } else {
    data.untrt$growth_fit <- growth.fit$fitted.values[df.in$trt.unit] / exp(coef)
  }
  
  growth.df <- data.untrt %>% # at individual unit level for each t
    filter(trt.unit, # filter to treated unit
           week >= T0/agg) %>% # post-intervention + the last period before intervention (which will be used as baseline)
    mutate(trt_post = (trt.unit & trt.time),
           # initialize the last period before intervention to baseline value and build on that, i.e., the observed incidence
           Y.untrt.growth = ifelse(trt_post, NA, inc))  %>%
    dplyr::select(unit, week, S_frac, trt.unit, trt.time, trt_post, inc, Y.untrt.growth, growth, growth_fit)
  
  # Translate to untreated potential outcomes on the case scale
  for (unit in unique(growth.df$unit)) {
    for (time in (T0/agg+1):((T0+T1)/agg)) { # for each of the post-intervention period
      growth.df$Y.untrt.growth[growth.df$unit==unit & growth.df$week==time] <- growth.df$Y.untrt.growth[growth.df$unit==unit & growth.df$week==time-1] *
        growth.df$growth_fit[growth.df$unit==unit & growth.df$week==time]
    }
  }
  data.untrt$Y.untrt.growth <- NA
  data.untrt$Y.untrt.growth[data.untrt$week >= T0/agg] <- growth.df$Y.untrt.growth
  
  # Bias correction
  for (i in trt.IDs) {
    for (time in (T0/agg+1):((T0+T1)/agg)) { # correct for each post intervention time period
      growth.var <- get_var(growth.fit, time, T0=T0/agg, id=i)
      data.untrt$Y.growth.adj1[data.untrt$unit==i & data.untrt$week==time] <- data.untrt$Y.untrt.growth[data.untrt$unit==i & data.untrt$week==time] / (1 + growth.var/2)
      data.untrt$Y.growth.adj2[data.untrt$unit==i & data.untrt$week==time] <- data.untrt$Y.untrt.growth[data.untrt$unit==i & data.untrt$week==time] / exp(growth.var/2)
    }
  }
  
  out <- data.frame(model="growth", effect=exp(tail(coef(growth.fit), 1)), p=stata.out$p, 
                    Y.untrt=mean(data.untrt$Y.untrt.growth[data.untrt$trt.time]),
                    Y.untrt.adj1=mean(data.untrt$Y.growth.adj1[data.untrt$trt.time]), 
                    Y.untrt.adj2=mean(data.untrt$Y.growth.adj2[data.untrt$trt.time]))
  rownames(out) <- NULL
  return(out)
}

# type is either "wt" for cohort-based or "est" for prevalence-based
# dgp is either "SIR or "SEIR"
run_Rt <- function(data.in, out.df, type, dgp, inf_mean, delta=NULL, inf_var=NULL, trt.IDs=1:N1, coef=NULL, parallel.id=0) {
  if (type=="wt") {
    data.in$Rt <- data.in$R_wt
  } else if (type=="est") {
    data.in$Rt <- data.in$R_est
  } # otherwise do not change Rt
  
  Rt.fit <- glm(Rt ~ -1 + factor(week) + factor(unit) + factor(trt_post), family = poisson, data = data.in)
  data.untrt <- data.in %>% filter(trt.unit) %>% mutate(trt_post = FALSE)
  command <- "glm Rt i.unit i.week 1.trt_post, family(poisson) link(log)
    boottest 1.trt_post, cluster(unit) reps(10000) quietly
    gen p = r(p) in 1
    keep p
    keep if _n==1"
  stata.out <- my_RStata(src=command, data.in=data.in, data.out=TRUE, stata.echo=FALSE, id=parallel.id)
  
  if (is.null(coef)) {
    # For log Rt models, first predict untreated Rt's using the fitted model
    data.untrt$Rt_fit <- predict(Rt.fit, newdata = data.untrt, type = "response") 
  } else {
    data.untrt$Rt_fit <- Rt.fit$fitted.values[df.in$trt.unit] / exp(coef)
  }
  
  # then de-aggregate the weekly data to daily, assuming constant transmission rate across the weekly window
  data.untrt.deagg <- out.df %>% filter(unit %in% trt.IDs, t <= T0 + T1 + burnin) %>%
    mutate(week = ceiling(t/agg), unit = factor(unit)) %>% # start with de-aggregated daily data
    merge(data.untrt %>% dplyr::select(unit, week, Rt, Rt_fit) %>% mutate(week = week + burnin/agg), 
          by = c("unit", "week"), all.x = T) %>%
    mutate(trt.time = (t > (T0+burnin)),
           # use the observed rate for before intervention
           Rt_fit = ifelse(trt.time, Rt_fit, Rt)) %>%
    filter(t > burnin)
  
  # Simulate case trajectories using untreated Rt
  Rt.untrt <- rbindlist(lapply(trt.IDs, function(ind) {
    I0 <- data.untrt.deagg$I[data.untrt.deagg$unit==ind & data.untrt.deagg$t==burnin+1]
    recovered <- data.untrt.deagg$R[data.untrt.deagg$unit==ind & data.untrt.deagg$t==burnin+1]
    S_frac <- data.untrt.deagg$S_frac[data.untrt.deagg$unit==ind]
    trans_prob <- data.untrt.deagg$Rt_fit[data.untrt.deagg$unit==ind] / S_frac # divide out S_frac to get back to beta scale
    if (dgp=="SIR") {
      run_SIR_varying(pop.size=pop.size, seeds=I0, recovered=recovered, trans_prob = trans_prob,
                      time_steps=(T0+T1), inf_mean=inf_mean, inf_var=inf_var)
    } else if (dgp=="SEIR") {
      E0 <- data.untrt.deagg$E[data.untrt.deagg$unit==ind & data.untrt.deagg$t==burnin+1]
      run_SEIR_varying(pop.size=pop.size, I0=I0, E0=E0, recovered=recovered, trans_prob = trans_prob,
                       time_steps=(T0+T1), inf_mean=inf_mean, delta=delta)
    }
  })) %>%
    mutate(unit=rep(trt.IDs, each=(T0+T1)), week = ceiling(t/agg)) %>%
    group_by(unit, week) %>%
    summarise(inc = sum(inc))
  data.untrt$Y.untrt.Rt <- Rt.untrt$inc
  
  # Bias correction
  for (i in trt.IDs) {
    for (time in (T0/agg+1):((T0+T1)/agg)) { # correct for each post intervention time period
      Rt.var <- get_var(Rt.fit, time, T0=T0/agg, id=i)
      data.untrt$Y.Rt.adj1[data.untrt$unit==i & data.untrt$week==time] <- data.untrt$Y.untrt.Rt[data.untrt$unit==i & data.untrt$week==time] / (1 + Rt.var/2)
      data.untrt$Y.Rt.adj2[data.untrt$unit==i & data.untrt$week==time] <- data.untrt$Y.untrt.Rt[data.untrt$unit==i & data.untrt$week==time] / exp(Rt.var/2)
    }
  }
  
  out <- data.frame(model=paste0("Rt_", type), effect=exp(tail(coef(Rt.fit), 1)), p=stata.out$p, 
                    Y.untrt=mean(data.untrt$Y.untrt.Rt[data.untrt$trt.time]),
                    Y.untrt.adj1=mean(data.untrt$Y.Rt.adj1[data.untrt$trt.time]), 
                    Y.untrt.adj2=mean(data.untrt$Y.Rt.adj2[data.untrt$trt.time]))
  rownames(out) <- NULL
  return(out)
}

run_beta <- function(data.in, out.df, dgp, inf_mean, delta=NULL, inf_var=NULL, trt.IDs=1:N1, coef=NULL, parallel.id=0) {
  beta.fit <- glm(beta_est ~ -1 + factor(week) + factor(unit) + factor(trt_post), family = poisson, data = data.in)
  data.untrt <- data.in %>% filter(trt.unit) %>% mutate(trt_post = FALSE)
  command <- "glm beta_est i.unit i.week 1.trt_post, family(poisson) link(log)
    boottest 1.trt_post, cluster(unit) reps(10000) quietly
    gen p = r(p) in 1
    keep p
    keep if _n==1"
  stata.out <- my_RStata(src=command, data.in=data.in, data.out=TRUE, stata.echo=FALSE, id=parallel.id)
  
  if (is.null(coef)) {
    # For log beta models, first de-aggregate the weekly data to daily, assuming constant transmission rate across the weekly window
    data.untrt$beta_fit <- predict(beta.fit, newdata = data.untrt, type = "response")
  } else {
    data.untrt$beta_fit <- beta.fit$fitted.values[df.in$trt.unit] / exp(coef)
  }
  
  data.untrt.deagg <- out.df %>% filter(unit %in% trt.IDs, t <= T0 + T1 + burnin) %>%
    mutate(week = ceiling(t/agg), unit = factor(unit)) %>% # start with de-aggregated daily data
    merge(data.untrt %>% dplyr::select(unit, week, beta_est, beta_fit) %>% mutate(week = week + burnin/agg), 
          by = c("unit", "week"), all.x = T) %>%
    mutate(trt.time = (t > (T0+burnin)),
           # use the observed rate for before intervention
           beta_fit = ifelse(trt.time, beta_fit, beta_est)) %>%
    filter(t > burnin)
  
  # Simulate case trajectories using untreated beta
  beta.untrt <- rbindlist(lapply(trt.IDs, function(ind) {
    I0 <- data.untrt.deagg$I[data.untrt.deagg$unit==ind & data.untrt.deagg$t==burnin+1]
    recovered <- data.untrt.deagg$R[data.untrt.deagg$unit==ind & data.untrt.deagg$t==burnin+1]
    trans_prob <- data.untrt.deagg$beta_fit[data.untrt.deagg$unit==ind]
    if (dgp=="SIR") {
      run_SIR_varying(pop.size=pop.size, seeds=I0, recovered=recovered, trans_prob = trans_prob,
                      time_steps=(T0+T1), inf_mean=inf_mean, inf_var=inf_var)
    } else if (dgp=="SEIR") {
      E0 <- data.untrt.deagg$E[data.untrt.deagg$unit==ind & data.untrt.deagg$t==burnin+1]
      run_SEIR_varying(pop.size=pop.size, I0=I0, E0=E0, recovered=recovered, trans_prob = trans_prob,
                       time_steps=(T0+T1), inf_mean=inf_mean, delta=delta)
    }
  })) %>%
    mutate(unit=rep(trt.IDs, each=(T0+T1)), week = ceiling(t/agg)) %>%
    group_by(unit, week) %>%
    summarise(inc = sum(inc))
  data.untrt$Y.untrt.beta <- beta.untrt$inc
  
  # Bias correction
  for (i in trt.IDs) {
    for (time in (T0/agg+1):((T0+T1)/agg)) { # correct for each post intervention time period
      beta.var <- get_var(beta.fit, time, T0=T0/agg, id=i)
      data.untrt$Y.beta.adj1[data.untrt$unit==i & data.untrt$week==time] <- data.untrt$Y.untrt.beta[data.untrt$unit==i & data.untrt$week==time] / (1 + beta.var/2)
      data.untrt$Y.beta.adj2[data.untrt$unit==i & data.untrt$week==time] <- data.untrt$Y.untrt.beta[data.untrt$unit==i & data.untrt$week==time] / exp(beta.var/2)
    }
  }
  
  out <- data.frame(model="beta", effect=exp(tail(coef(beta.fit), 1)), p=stata.out$p, 
                    Y.untrt=mean(data.untrt$Y.untrt.beta[data.untrt$trt.time]),
                    Y.untrt.adj1=mean(data.untrt$Y.beta.adj1[data.untrt$trt.time]), 
                    Y.untrt.adj2=mean(data.untrt$Y.beta.adj2[data.untrt$trt.time]))
  rownames(out) <- NULL
  return(out)
}

# Calculate the true untreated incidence counterfactual
run_true <- function(out.df, trans_prob.base1, dgp, trt.IDs=1:N1) {
  true.untrt <- rbindlist(lapply(trt.IDs, function(ind) {
    I0 <- out.df$I[out.df$unit==ind & out.df$t==burnin+1]
    recovered <- out.df$R[out.df$unit==ind & out.df$t==burnin+1]
    if (dgp=="SIR") {
      run_SIR_varying(pop.size=pop.size, seeds=I0, recovered=recovered, time_steps=(T0+T1), 
                      inf_mean=inf_mean, trans_prob=rep(trans_prob.base1, (T0+T1)))
    } else if (dgp=="SEIR") {
      E0 <- out.df$E[out.df$unit==ind & out.df$t==burnin+1]
      run_SEIR_varying(pop.size=pop.size, I0=I0, E0=E0, recovered=recovered, time_steps=(T0+T1), 
                       inf_mean=inf_mean, delta=delta, trans_prob=rep(trans_prob.base1, (T0+T1)))
    }
  })) %>%
    mutate(unit=rep(trt.IDs, each=(T0+T1)), week = ceiling(t/agg)) %>%
    group_by(unit, week) %>%
    summarise(inc = sum(inc))
  Y.untrt.true <- mean(true.untrt$inc[true.untrt$week > (T0/agg)])
  return(Y.untrt.true)
}

# Helper functions for School Masking example
loginc_AME <- function(coef, subset=NULL) {
  ATT_gt.tmp <- data.frame(time_to_trt = time_to_trt, coef = coef)
  loginc.obs <- df_sunab %>% mutate(time_to_trt = time - group)
  
  if (!is.null(subset)) {
    loginc.obs <- loginc.obs %>% filter(time_to_trt %in% subset)
  } else {
    loginc.obs <- loginc.obs %>% filter(time_to_trt >= 0)
  }
    
  loginc.obs <- loginc.obs %>%
    merge(ATT_gt.tmp, by = "time_to_trt") %>%
    mutate(y.ctl = y/exp(coef), # fitted control potential outcome on case scale, cannot convert on the original log scale b/c of zeros
           diff = y - y.ctl) %>% # difference to get marginal effect for each unit i to ME
    group_by(time_to_trt) %>%
    summarise(diff = mean(diff)) # average over all units to get the average marginal effect
  sum(loginc.obs$diff)
}

growth_AME <- function(coef, subset=NULL) {
  ATT_gt.tmp <- data.frame(time_to_trt = time_to_trt, coef = coef)
  
  df.all <- df_sunab %>% mutate(time_to_trt = time - group)
  df <- df.all %>%
    filter(time_to_trt >= -1) %>%
    merge(ATT_gt.tmp, by = "time_to_trt", all.x = T) %>%
    mutate(coef = ifelse(time_to_trt==-1, 0, coef),
           growth.ctl = log(y) - coef,
           Pos.ctl = NA)
  
  for (unit in unique(df$ID)) {
    for (time in sort(unique(df$time[df$ID==unit]))) {
      if ((time+1) == unique(df$group[df$ID==unit])) { # last period before intervention
        inc.last <- df$PosPer1K[df$ID==unit & df$time==time]
        if (inc.last==0) {
          df$Pos.ctl[df$ID==unit & df$time==time] <- inc.last
        } else {
          df$Pos.ctl[df$ID==unit & df$time==time] <- mean(df.all$PosPer1K[df.all$ID==unit & df.all$time %in% (time-4):time])
        }
      } else { # recover the untreated trajectory for the treated group
        time.ind <- ifelse(nrow(df[df$ID==unit & df$time==time-1,])==0, 2, 1)
        df$Pos.ctl[df$ID==unit & df$time==time] <- df$Pos.ctl[df$ID==unit & df$time==time-time.ind] * 
          exp(df$growth.ctl[df$ID==unit & df$time==time])
      }
    }
  }
  
  if (!is.null(subset)) df <- df %>% filter(time_to_trt %in% subset)
  
  df <- df %>% filter(time >= group) %>%
    mutate(diff = PosPer1K - Pos.ctl) %>%
    group_by(ID) %>% summarise(diff = sum(diff)) # difference to get marginal effect for each unit
  mean(df$diff) # average over all units to get the AME
}
