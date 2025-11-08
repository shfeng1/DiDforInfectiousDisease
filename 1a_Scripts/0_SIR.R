#### PURPOSE:
#### Define SIR data-generating process
#### and simulation wrapper functions

source("./1a_Scripts/0_RStata.R")
source("./1a_Scripts/0_Bias_Correction.R")
source("./1a_Scripts/0_Estimate_Rt.R")
source("./1a_Scripts/0_Run_Estimators.R")

#### FUNCTION #1 ####
#### Define SIR data-generating process
run_SIR_varying <- function(
    time_steps,      # number of time‐steps you want to return
    pop.size,        # total population
    seeds,           # initial infection
    recovered = 0,   # initial recovered
    trans_prob,      # vector of beta_t, length >= time_steps + 1

    inf_mean,        # mean infectious period (required)
    inf_var  = NULL  # variance (optional; if NULL ⇒ geometric)
) {
  # fix simulation by one extra time step
  Ttot <- time_steps + 1
  trans_prob <- c(NA, trans_prob)
  
  # 1) get pmf for discrete gamma
  if (!is.null(inf_var)) { # if inf_var is specified => use discrete-Gamma
    shape   <- inf_mean^2 / inf_var
    scale   <- inf_var / inf_mean
    F_gamma <- pgamma(1:Ttot, shape = shape, scale = scale)
    F_lag   <- c(0, head(F_gamma, -1))
    pmf     <- F_gamma - F_lag
  }
  
  # 2) initialize states and variables
  S <- I <- R <- inc <- Rt_true <- beta_true <- numeric(Ttot)
  
  # 3) initial conditions at t = 1
  S[1]         <- pop.size - seeds - recovered
  I[1]         <- seeds
  R[1]         <- recovered
  inc[1]       <- seeds
  Rt_true[1]   <- NA
  beta_true[1] <- NA
  
  # 4) simulate t = 2 .. Ttot
  for (t in 2:Ttot) {
    beta_t  <- trans_prob[t]
    
    # draw new infections
    lambda_t <- beta_t * I[t-1] * (S[t-1] / pop.size)
    trans_t  <- rpois(1, lambda_t)
    
    # update susceptibles and incidence
    S[t]   <- S[t-1] - trans_t
    inc[t] <- trans_t
    
    if (is.null(inf_var)) {
      # geometric removal dynamics
      I[t] <- (1 - 1/inf_mean) * I[t-1] + trans_t
      R[t] <- R[t-1] + 1/inf_mean * I[t-1]
    } else {
      # discrete‐Gamma removal dynamics
      prev_incs  <- inc[t - (1:(t-1))]
      new_removals <- sum(pmf[1:(t-1)] * prev_incs)
      
      I[t] <- I[t-1] + trans_t - new_removals
      R[t] <- R[t-1] + new_removals
    }
    
    # compute true R_t and beta_t
    Rt_true[t]   <- if (I[t-1] > 0) inc[t] / I[t-1] else NA
    beta_true[t] <- Rt_true[t] / (S[t-1] / pop.size)
  }
  
  # 5) drop the first (burn‐in) step, return exactly time_steps rows
  data.table(t          = 1:time_steps,
             S          = S[-1],
             I          = I[-1],
             R          = R[-1],
             inc        = inc[-1],
             Rt_true    = Rt_true[-1],
             beta_true  = beta_true[-1]) %>%
    mutate(S_frac = S/pop.size)
}

#### FUNCTION #2 ####
#### 1) Simulate data according to SIR
#### 2) Calculate estimators Rt: W-T + prev;  beta_t: (prev Rt) / (St/N)
#### 3) Fit DiD models 
#### 4) Pull point estimates + calculate AME
#### 5) Bias correction
SIR_sim <- function(pop.size, N, N1, T0, T1, burnin, seed1, seed2, inf_mean,
                    trans_prob.base1, trans_prob.base2, eff.multi1, parallel.id=0) {
  parallel.id <- paste0("SIR", parallel.id)
  out.df <- gen_SIR(trans_prob.base1, trans_prob.base2, eff.multi1, inf_mean) # simulate data according SIR
  data.in <- process_data(out.df, inf_mean, agg, dgp="SIR") # estimate R_t, beta_t, aggregate to weekly level
  ################################################################################################################################
  inc.out <- run_inc(data.in, parallel.id)
  loginc.out <- run_loginc(data.in, parallel.id)
  growth.out <- run_growth(data.in, parallel.id)
  Rt_wt.out <- run_Rt(data.in, out.df, type="wt", dgp="SIR", inf_mean=inf_mean, parallel.id=parallel.id)
  Rt_est.out <- run_Rt(data.in, out.df, type="est", dgp="SIR", inf_mean=inf_mean, parallel.id=parallel.id)
  beta.out <- run_beta(data.in, out.df, dgp="SIR", inf_mean=inf_mean, parallel.id=parallel.id)
  Y.untrt.true <- run_true(out.df, trans_prob.base1, dgp="SIR")
  ############################################################################################################################
  # summarize outputs
  out <- rbind(inc.out, loginc.out, growth.out, Rt_wt.out, Rt_est.out, beta.out) %>% 
    mutate(N=N, N1=N1, trans_prob.base1=trans_prob.base1, trans_prob.base2=trans_prob.base2, pop.size=pop.size, seed=seed1,
           eff.multi=eff.multi1, burnin=burnin, T0=T0, T1=T1, 
           S_frac.mean=mean(data.in$S_frac[data.in$trt.time]), S_frac.min=min(data.in$S_frac),
           Y.trt=mean(data.in$inc[data.in$trt_post]), Y.untrt.true=Y.untrt.true)
  # true AME * T1 / pop.size is the effect size as % population
  
  return(out)
}

#### FUNCTION #3 ####
#### 1) Simulate data according to SIR
#### 2) Calculate estimators Rt: W-T + prev;  beta_t: prev Rt/St/N
#### 3) Get p-values from wild score bootstrap
#### 4) Get p-values from normal-based standard error approach
inference_sim <- function(pop.size, N, N1, T0, T1, burnin, seed1, seed2,
                          trans_prob.base1, trans_prob.base2, inf_mean,
                          eff.multi1=1, eff.multi2=1, parallel.id=0) {
  parallel.id <- paste0("Inference", parallel.id)
  out.df <- gen_SIR(trans_prob.base1, trans_prob.base2, eff.multi1, inf_mean)
  data.in <- process_data(out.df, inf_mean, agg, dgp="SIR")
  ############################################################################################################################
  # to get p-value from wild score bootstrap
  stata.out <- data.frame()
  for (var in c("inc", "loginc", "growth", "R_wt", "R_est", "beta_est")) {
    if (var == "inc") {
      command <- paste0("glm inc i.unit i.week i.trt_post, family(gaussian) link(identity)
    boottest 1.trt_post, cluster(unit) reps(10000) quietly
    gen p = r(p) in 1
    keep p
    keep if _n==1")
    } else if (var == "loginc") {
      command <- paste0("glm inc i.unit i.week 1.trt_post, family(poisson) link(log)
    boottest 1.trt_post, cluster(unit) reps(10000) quietly
    gen p = r(p) in 1
    keep p
    keep if _n==1")
    } else {
      command <- paste0("glm ", var, " i.unit i.week 1.trt_post, family(poisson) link(log)
    boottest 1.trt_post, cluster(unit) reps(10000) quietly
    gen p = r(p) in 1
    keep p
    keep if _n==1")
    }
    out.tmp <- my_RStata(src=command, data.in=data.in, data.out=TRUE, stata.echo=FALSE, id=parallel.id)
    out.tmp$var <- var
    stata.out <- rbind(stata.out, out.tmp)
  }
  stata.out <- rename(stata.out, wild = p)
  ############################################################################################################################
  # to get p-value from normal-based clustered standard error
  normal.out <- data.frame()
  for (var in c("inc", "loginc", "growth", "R_wt", "R_est", "beta_est")) {
    if (var == "inc") {
      fit <- lm(inc ~ factor(unit) + factor(week) + trt_post, data = data.in)
    } else if (var == "loginc") {
      fit <- glm(inc ~ factor(unit) + factor(week) + trt_post, family = poisson(), data = data.in)
    } else {
      fit <- glm(as.formula(paste0(var, " ~ factor(unit) + factor(week) + trt_post")), family = poisson(), data = data.in)
    }
    out.tmp <- data.frame(var = var, coef = as.numeric(tail(coef(fit), 1)))
    out.tmp$normal <- tail(coeftest(fit, vcov = vcovCL(fit, cluster = ~unit)), 1)[,4] # p-value
    normal.out <- rbind(normal.out, out.tmp)
  }
  out <- merge(stata.out, normal.out) %>%
    mutate(N=N, N1=N1, trans_prob.base1=trans_prob.base1, trans_prob.base2=trans_prob.base2, pop.size=pop.size, seed=seed1, 
           eff.multi=eff.multi1, burnin=burnin, T0=T0, T1=T1)

  out
}

#### FUNCTION #4 ####
#### 1) Simulate data according to SIR
#### 2) Simulate untreated potential outcomes according to SIR
#### 3) Calculate what the true RR would be
SIR_true_eff <- function(trans_prob.base1, trans_prob.base2, eff.multi1) {
  out.df <- gen_SIR(trans_prob.base1, trans_prob.base2, eff.multi1, inf_mean)
  data.in <- process_data(out.df, inf_mean, agg, dgp="SIR")
  
  # Construct the untreated counterfactual
  untrt.df <- lapply(1:N1, function(ind) { 
    run_SIR_varying(pop.size=pop.size, time_steps=(T0+T1+burnin*3), seeds=seed1, inf_mean=inf_mean,
                    trans_prob = rep(trans_prob.base1, (T0+T1+burnin*3)))}) %>% rbindlist() %>% 
    mutate(unit=rep((1:N1), each=(T0+T1+burnin*3)),
           prevalence = compute_prevalence(inf_mean=inf_mean, ID=unit, inc=inc, time=t, Ttot=T0+T1+burnin),
           R_cohort = compute_Rt_cohort(inc, unit, inf_mean, (T0+T1+burnin*3)))
  untrt.data.in <- process_data(untrt.df, inf_mean, agg, dgp="SIR")
  
  out <- data.frame(trans_prob.base1=trans_prob.base1, trans_prob.base2=trans_prob.base2,
                    eff.multi = eff.multi1, model = c("inc", "loginc", "growth", "Rt_wt", "Rt_cohort", "Rt_est", "Rt_true", "beta"),
                    eff.true = c(mean(data.in$inc[data.in$trt_post])-mean(untrt.data.in$inc[untrt.data.in$trt_post]),
                                 mean(data.in$inc[data.in$trt_post])/mean(untrt.data.in$inc[untrt.data.in$trt_post]),
                                 mean(data.in$growth[data.in$trt_post])/mean(untrt.data.in$growth[untrt.data.in$trt_post]),
                                 mean(data.in$R_wt[data.in$trt_post]) / mean(untrt.data.in$R_wt[untrt.data.in$trt_post]),
                                 mean(data.in$R_cohort[data.in$trt_post]) / mean(untrt.data.in$R_cohort[untrt.data.in$trt_post]),
                                 mean(data.in$R_est[data.in$trt_post]) / mean(untrt.data.in$R_est[untrt.data.in$trt_post]),
                                 mean(data.in$R_true[data.in$trt_post]) / mean(untrt.data.in$R_true[untrt.data.in$trt_post]),
                                 eff.multi1))
  out
}
