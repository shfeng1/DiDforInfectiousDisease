#### PURPOSE:
#### Define SEIR data-generating process
#### and simulation wrapper functions

source("./1a_Scripts/0_Estimate_Rt.R")
source("./1a_Scripts/0_Run_Estimators.R")

if (!("incubation_mean" %in% names(formals(process_data)))) {
  stop("Incompatible 0_Estimate_Rt.R: process_data() must include incubation_mean.")
}

#### FUNCTION #1 ####
#### Define SEIR data-generating process
run_SEIR_varying = function(
    time_steps, # number of total timesteps
    pop.size, # population size
    I0, E0=0, # number of initial infections
    recovered = 0, # number of initial recovered
    S0 = NULL, # optional susceptible state at simulation start
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
  S[1] = if (is.null(S0)) pop.size - E0 - I0 - recovered else as.numeric(S0)
  if (length(S[1]) != 1L || !is.finite(S[1]) || S[1] < 0 || S[1] > pop.size) {
    stop("S0 must be NULL or one finite value between 0 and pop.size.")
  }
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
#### 2) Calculate estimators Rt_exposure and beta_exposure
#### 3) Fit DiD models 
#### 4) Pull point estimates + calculate AME
#### 5) Bias correction
SEIR_sim <- function(pop.size, N, N1, T0, T1, burnin, seed1, seed2, inf_mean, delta,
                     trans_prob.base1, trans_prob.base2, eff.multi1, parallel.id=0,
                     simulate_from_trt=FALSE, end_buffer=NULL,
                     calculate_p=TRUE, return_data=FALSE) {
  if (is.null(end_buffer)) end_buffer <- get("end_buffer", envir=.GlobalEnv)
  parallel.id <- paste0("SEIR", parallel.id)
  out.df <- gen_SEIR(trans_prob.base1, trans_prob.base2, eff.multi1, inf_mean,
                     delta=delta, N=N, N1=N1, pop.size=pop.size,
                     seed1=seed1, seed2=seed2, T0=T0, T1=T1,
                     burnin=burnin, end_buffer=end_buffer)
  data.in <- process_data(out.df, inf_mean, agg, dgp="SEIR", incubation_mean=delta,
                          discard_start=burnin, discard_end=end_buffer)
  ################################################################################################################################
  inc.out <- run_inc(data.in, parallel.id, calculate_p=calculate_p)
  loginc.out <- run_loginc(data.in, parallel.id, calculate_p=calculate_p)
  growth.out <- run_growth(data.in, parallel.id, calculate_p=calculate_p)
  Rt_exposure.out <- run_Rt(data.in, out.df, dgp="SEIR", inf_mean=inf_mean, delta=delta, parallel.id=parallel.id,
                            calculate_p=calculate_p,
                            simulate_from_trt=simulate_from_trt)
  beta.out <- run_beta(data.in, out.df, dgp="SEIR", inf_mean=inf_mean, delta=delta, parallel.id=parallel.id,
                       calculate_p=calculate_p,
                       simulate_from_trt=simulate_from_trt)
  Y.untrt.true <- run_true(out.df, trans_prob.base1, dgp="SEIR")
  ############################################################################################################################
  # summarize outputs
  out <- rbind(inc.out, loginc.out, growth.out, Rt_exposure.out, beta.out) %>%
    mutate(N=N, N1=N1, trans_prob.base1=trans_prob.base1, trans_prob.base2=trans_prob.base2, pop.size=pop.size, seed=seed1,
           eff.multi=eff.multi1, burnin=burnin, T0=T0, T1=T1, 
           S_frac.mean=mean(data.in$S_frac[data.in$trt.time]), S_frac.min=min(data.in$S_frac),
           Y.trt=mean(data.in$inc[data.in$trt_post]), Y.untrt.true=Y.untrt.true)
  if (return_data) return(list(result=out, bootstrap_data=data.in))
  out
}

#### FUNCTION #3 ####
#### 1) Simulate data according to SEIR
#### 2) Simulate untreated potential outcomes according to SEIR
#### 3) Calculate what the true RR would be
SEIR_true_eff <- function(trans_prob.base1, trans_prob.base2, eff.multi1) {
  assert_scalar_parameter <- function(x, name) {
    if (!is.numeric(x) || length(x)!=1L || is.na(x) || !is.finite(x) || x<=0) {
      stop(name, " must be one finite numeric scalar greater than zero.")
    }
  }
  assert_scalar_parameter(trans_prob.base1, "trans_prob.base1")
  assert_scalar_parameter(trans_prob.base2, "trans_prob.base2")
  assert_scalar_parameter(eff.multi1, "eff.multi1")

  out.df <- gen_SEIR(trans_prob.base1, trans_prob.base2, eff.multi1, inf_mean)
  data.in <- process_data(out.df, inf_mean, agg, dgp="SEIR", incubation_mean=delta)

  total.days <- T0+T1+burnin+end_buffer
  untrt.df <- lapply(1:N1, function(ind) {
    run_SEIR_varying(pop.size=pop.size, time_steps=total.days,
                     I0=seed1, inf_mean=inf_mean, delta=delta,
                     trans_prob=rep(trans_prob.base1, total.days))
  }) %>% rbindlist() %>%
    mutate(unit=rep(1:N1, each=total.days), initial_seed=seed1)
  untrt.df$prevalence <- compute_prevalence_seeded(untrt.df, inf_mean, total.days)
  untrt.df$infected_est <- compute_infected_seeded(untrt.df, delta, total.days)
  untrt.data.in <- process_data(untrt.df, inf_mean, agg, dgp="SEIR", incubation_mean=delta)

  finite_mean <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) NA else mean(x)
  }
  finite_ratio <- function(x, y) {
    x <- finite_mean(x); y <- finite_mean(y)
    if (!is.finite(x) || !is.finite(y) || y==0) NA else x/y
  }
  trt.idx <- data.in$trt_post
  untrt.idx <- untrt.data.in$trt_post
  inc.trt <- finite_mean(data.in$inc[trt.idx])
  inc.untrt <- finite_mean(untrt.data.in$inc[untrt.idx])
  model.names <- c("inc", "loginc", "growth", "Rt_exposure",
                   "Rt_true", "beta_exposure")
  true.effects <- c(
    if (is.finite(inc.trt) && is.finite(inc.untrt)) inc.trt-inc.untrt else NA,
    finite_ratio(data.in$inc[trt.idx], untrt.data.in$inc[untrt.idx]),
    finite_ratio(data.in$growth[trt.idx], untrt.data.in$growth[untrt.idx]),
    finite_ratio(data.in$Rt_exposure[trt.idx], untrt.data.in$Rt_exposure[untrt.idx]),
    finite_ratio(data.in$R_true[trt.idx], untrt.data.in$R_true[untrt.idx]),
    eff.multi1)
  if (length(true.effects)!=length(model.names)) stop("Internal SEIR true-effect output length error.")
  data.frame(trans_prob.base1=rep(trans_prob.base1, length(model.names)),
             trans_prob.base2=rep(trans_prob.base2, length(model.names)),
             eff.multi=rep(eff.multi1, length(model.names)),
             model=model.names, eff.true=true.effects, stringsAsFactors=FALSE)
}
