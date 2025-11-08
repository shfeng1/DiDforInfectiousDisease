#### PURPOSE:
#### Define SEIR data-generating process
#### and simulation wrapper functions

source("./1a_Scripts/0_Estimate_Rt.R")
source("./1a_Scripts/0_Run_Estimators.R")

#### FUNCTION #1 ####
#### Define SEIR data-generating process
run_SEIR_varying = function(
    time_steps, # number of total timesteps
    pop.size, # population size
    I0, E0=0, # number of initial infections
    recovered = 0, # number of initial recovered
    death = 0, # number of initial death
    death_prob = 0.01, # probability of deaths
    trans_prob, # a vector probability of transmission given contact
    inf_mean, # average days of infectiousness
    delta # average days of incubation period
){
  time_steps <- time_steps + 1
  trans_prob <- c(NA, trans_prob)
  
  # track states
  S = rep(0, time_steps)
  I = rep(0, time_steps)
  E = rep(0, time_steps)
  R = rep(0, time_steps)
  Rt = rep(0, time_steps)
  inc = rep(0, time_steps)
  infected = rep(0, time_steps)
  mean = rep(0, time_steps)
  Rt_true = rep(NA, time_steps) # true Rt
  beta_true = rep(0, time_steps) # true beta t
  
  # initial conditions
  S[1] = pop.size - E0 - I0 - recovered
  E[1] = E0
  I[1] = I0
  R[1] = recovered
  
  for(i in 2:time_steps){
    beta <- trans_prob[i]
    
    # set up random draw
    mean[i] = beta*I[i-1]*S[i-1]/pop.size
    trans_t = rpois(1, lambda = mean[i])
    
    # susceptible
    S[i] = S[i-1] - trans_t
    
    # infected & infectious
    infected[i] = trans_t # individuals become infected on day i; unobservable
    inc[i] = 1/delta*E[i-1] # incidence; observed
    E[i] =  (1-1/delta)*E[i-1] + trans_t
    I[i] = (1-1/inf_mean)*I[i-1] + 1/delta*E[i-1]
    
    # recovered
    R[i] = R[i-1] + 1/inf_mean*I[i-1]
    
    # true Rt
    Rt_true[i] = infected[i] / I[i-1]
    
    # true beta t
    beta_true[i] = Rt_true[i]  / (S[i-1]/pop.size)
  }
  
  d = data.table(trans_prob, S, E, I, R, Rt_true, beta_true, inc, infected, mean)[2:time_steps,] %>%
    mutate(S_frac = S/pop.size, t = row_number())
  return(d)
}

#### FUNCTION #2 ####
#### 1) Simulate data according to SEIR
#### 2) Calculate estimators Rt: W-T + prev;  beta_t: (prev Rt) / (St/N)
#### 3) Fit DiD models 
#### 4) Pull point estimates + calculate AME
#### 5) Bias correction
SEIR_sim <- function(pop.size, N, N1, T0, T1, burnin, seed1, seed2, inf_mean, delta,
                     trans_prob.base1, trans_prob.base2, eff.multi1,  parallel.id=0) {
  parallel.id <- paste0("SEIR", parallel.id)
  out.df <- gen_SEIR(trans_prob.base1, trans_prob.base2, eff.multi1, inf_mean)
  data.in <- process_data(out.df, inf_mean, agg, dgp="SEIR")
  ################################################################################################################################
  inc.out <- run_inc(data.in, parallel.id)
  loginc.out <- run_loginc(data.in, parallel.id)
  growth.out <- run_growth(data.in, parallel.id)
  Rt_wt.out <- run_Rt(data.in, out.df, type="wt", dgp="SEIR", inf_mean=inf_mean, delta=delta, parallel.id=parallel.id)
  Rt_est.out <- run_Rt(data.in, out.df, type="est", dgp="SEIR", inf_mean=inf_mean, delta=delta, parallel.id=parallel.id)
  beta.out <- run_beta(data.in, out.df, dgp="SEIR", inf_mean=inf_mean, delta=delta, parallel.id=parallel.id)
  Y.untrt.true <- run_true(out.df, trans_prob.base1, dgp="SEIR")
  ############################################################################################################################
  # summarize outputs
  out <- rbind(inc.out, loginc.out, growth.out, Rt_wt.out, Rt_est.out, beta.out) %>% 
    mutate(N=N, N1=N1, trans_prob.base1=trans_prob.base1, trans_prob.base2=trans_prob.base2, pop.size=pop.size, seed=seed1,
           eff.multi=eff.multi1, burnin=burnin, T0=T0, T1=T1, 
           S_frac.mean=mean(data.in$S_frac[data.in$trt.time]), S_frac.min=min(data.in$S_frac),
           Y.trt=mean(data.in$inc[data.in$trt_post]), Y.untrt.true=Y.untrt.true)
  return(out)
}

#### FUNCTION #3 ####
#### 1) Simulate data according to SEIR
#### 2) Simulate untreated potential outcomes according to SEIR
#### 3) Calculate what the true RR would be
SEIR_true_eff <- function(trans_prob.base1, trans_prob.base2, eff.multi1) {
  out.df <- gen_SEIR(trans_prob.base1, trans_prob.base2, eff.multi1, inf_mean)
  data.in <- process_data(out.df, inf_mean, agg, dgp="SEIR")
  
  # Construct the untreated counterfactual
  untrt.df <- lapply(1:N1, function(ind) { 
    run_SEIR_varying(pop.size=pop.size, time_steps=(T0+T1+burnin*3), I0=seed1, inf_mean=inf_mean, delta=delta,
                     trans_prob = rep(trans_prob.base1, (T0+T1+burnin*3)))}) %>% rbindlist() %>% 
    mutate(unit=rep((1:N1), each=(T0+T1+burnin*3)),
           prevalence = compute_prevalence(inf_mean=inf_mean, ID=unit, inc=inc, time=t, Ttot=T0+T1+burnin),
           infected_est = compute_infected(delta=delta, ID=unit, inc=inc, time=t, Ttot=T0+T1+burnin),
           R_cohort = compute_Rt_cohort(inc, unit, inf_mean, (T0+T1+burnin*3)))
  untrt.data.in <- process_data(untrt.df, inf_mean, agg, dgp="SEIR")
  ################################################################################################################################
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
