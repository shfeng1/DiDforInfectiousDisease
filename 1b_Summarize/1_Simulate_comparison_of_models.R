rm(list=ls())
here::i_am("1b_Summarize/1_Simulate_comparison_of_models.R")
source("./global_options.R")
source("./1a_Scripts/0_SIR.R")
##############################################################################################################################
# set up variable inputs
T.total <- 80
burnin <- 25
beta <- 0.105 # baseline beta
beta2 <- 0.20 # to change for different transmission
beta3 <- 0.16 # to change for non-trivial susceptible depletion
inf_days <- 10
pop <- 1e+08 # make it large enough so there is NO susceptible depletion in other cases
pop2 <- 1e+04 # to force a non-trivial susceptible depletion where needed
agg <- 2 # days of aggregation / smoothing
##############################################################################################################################
foreach::getDoParWorkers()

# (a) All parameters match
set.seed(2026, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
test1.tmp <- foreach(k=1:1000, .errorhandling = "stop", .combine=function(x,y) rbindlist(list(x,y))) %dopar% {
  sim_model_compare(pop.size1=pop, pop.size2=pop, T.total=T.total, burnin=burnin, T1=0, inf_days=inf_days,
            trans_prob.base1=beta, trans_prob.base2=beta) %>% mutate(sim=k)
}

# (b) Different numbers of initial infections
set.seed(2026, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
test2.tmp <- foreach(k=1:1000, .errorhandling = "stop", .combine=function(x,y) rbindlist(list(x,y))) %dopar% {
  sim_model_compare(seed1=200, pop.size1=pop, pop.size2=pop, T.total=T.total, burnin=burnin, T1=0, inf_days=inf_days,
            trans_prob.base1=beta, trans_prob.base2=beta) %>% mutate(sim=k)
}

# (c) Different baseline transmission, but ratio was constant over time
set.seed(2026, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
test3.tmp <- foreach(k=1:1000, .errorhandling = "stop", .combine=function(x,y) rbindlist(list(x,y))) %dopar% {
  sim_model_compare(trans_prob.base1=beta2, pop.size1=pop, pop.size2=pop, T.total=T.total, burnin=burnin,
            T1=30,
            inf_days=inf_days, trans_prob.base2=beta) %>% mutate(sim=k)
}

# (d) Different time-varying transmission, doubled at time 45
set.seed(2026, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
test4.tmp <- foreach(k=1:1000, .errorhandling = "stop", .combine=function(x,y) rbindlist(list(x,y))) %dopar% {
  sim_model_compare(trans_prob.base1=beta2, pop.size1=pop, pop.size2=pop, T.total=T.total, burnin=burnin, 
            T1=30,
            inf_days=inf_days, trans_prob.base2=beta,
            eff.multi1 = 2, eff.multi2 = 2) %>% mutate(sim=k)
}

# (e) Different transmission & non-trivial susceptible depletion
set.seed(2026, kind = "L'Ecuyer-CMRG") # set seed properly for %dopar%
test5.tmp <- foreach(k=1:1000, .errorhandling = "stop", .combine=function(x,y) rbindlist(list(x,y))) %dopar% {
  sim_model_compare(pop.size1=pop2, pop.size2=pop, T.total=T.total, burnin=burnin, T1=0, inf_days=inf_days,
            trans_prob.base1=beta3, trans_prob.base2=beta) %>% mutate(sim=k)
}
##############################################################################################################################
# Summarize / Average over simulations
test1 <- test1.tmp %>%
  group_by(t, unit, pop.size) %>%
  reframe(S = mean(S), I = mean(I), R = mean(R), inc = mean(inc), C = mean(C), mean = mean(mean),
          id = "'(a) All parameters match'", lab = "")
test2 <- test2.tmp %>%
  group_by(t, unit, pop.size) %>%
  reframe(S = mean(S), I = mean(I), R = mean(R), inc = mean(inc), C = mean(C), mean = mean(mean),
          id = "'(b) Different initial\ninfections'", lab = "I[1*`,`*0]~`=`~200")
test3 <- test3.tmp %>%
  group_by(t, unit, pop.size) %>%
  reframe(S = mean(S), I = mean(I), R = mean(R), inc = mean(inc), C = mean(C), mean = mean(mean),
          id = "'(c) Different transmission'", lab = paste0("beta[1*`,`*t]~`=`~",beta2))
test4 <- test4.tmp %>%
  group_by(t, unit, pop.size) %>%
  reframe(S = mean(S), I = mean(I), R = mean(R), inc = mean(inc), C = mean(C), mean = mean(mean),
          id = "'(d) Different time-varying\ntransmission'", lab = 'beta[0*`,`*t]~`&`~beta[1*`,`*t]~"doubled at time 45"')
test5 <- test5.tmp %>%
  group_by(t, unit, pop.size) %>%
  reframe(S = mean(S), I = mean(I), R = mean(R), inc = mean(inc),  C = mean(C), mean = mean(mean),
          id = "'(e) Different transmission and\nnon-trivial susceptible depletion'", 
          lab = paste0("N[1]*`=`*10*`,`*0*0*0~`&`~beta[1*`,`*t]~`=`~", beta3))
##############################################################################################################################
# calculate the outcomes in each model
# some numbers are scaled down to make figure look better for illustration purposes; otherwise all other lines are shown as a flat line at 0
inc_df <- rbind(test1, test2, test3, test4, test5) %>%  
  mutate(y = inc, 
         y = ifelse(id=="'(c) Different transmission'", y/500, y),
         y = ifelse(id=="'(d) Different time-varying\ntransmission'", y/1500, y),
         y = ifelse(id=="'(e) Different transmission and\nnon-trivial susceptible depletion'", y/3, y),
         model = "incidence") %>%
  dplyr::select(model, unit, id, lab, y, t)

loginc_df <- rbind(test1, test2, test3, test4, test5) %>% 
  mutate(y = ifelse(inc > 0, log(inc), 0),
         y = ifelse(id=="'(c) Different transmission'", y/2, y),
         y = ifelse(id=="'(d) Different time-varying\ntransmission'", y/2, y),
         model = "'log incidence'") %>%
  dplyr::select(model, unit, id, lab, y, t)

growth_df <- rbind(test1, test2, test3, test4, test5) %>%
  mutate(week = ceiling(t/agg)) %>%
  group_by(unit, id, lab, week) %>%
  summarize(inc = sum(inc)) %>% 
  mutate(y = log(inc) - log(lag(inc, 1)),
         y = ifelse(id=="'(d) Different time-varying\ntransmission'", y/2, y),
         t = week*agg,
         model = "'log growth'") %>%
  dplyr::select(model, unit, id, lab, y, t)

calculate_exposure_outcomes <- function(data, id, lab) {
  data %>%
    group_by(sim, unit) %>% arrange(t, .by_group=TRUE) %>%
    mutate(prevalence=Reduce(function(previous, current) current+(1-1/inf_days)*previous,
                             inc[-1], init=first(I), accumulate=TRUE),
           prevalence_lag=lag(prevalence), S_lag=lag(S), week=ceiling(t/agg)) %>%
    filter(t>=2) %>% group_by(sim, unit, week) %>%
    summarise(Rt_exposure=sum(inc)/sum(prevalence_lag),
              S_frac=sum(S_lag)/(first(pop.size)*agg), .groups="drop") %>%
    mutate(beta_exposure=Rt_exposure/S_frac, t=week*agg, id=id, lab=lab)
}

exposure_df <- bind_rows(
  calculate_exposure_outcomes(test1.tmp, test1$id[[1]], test1$lab[[1]]),
  calculate_exposure_outcomes(test2.tmp, test2$id[[1]], test2$lab[[1]]),
  calculate_exposure_outcomes(test3.tmp, test3$id[[1]], test3$lab[[1]]),
  calculate_exposure_outcomes(test4.tmp, test4$id[[1]], test4$lab[[1]]),
  calculate_exposure_outcomes(test5.tmp, test5$id[[1]], test5$lab[[1]])
) %>%
  group_by(id, lab, unit, t) %>%
  summarise(Rt_exposure=mean(Rt_exposure), beta_exposure=mean(beta_exposure), .groups="drop") %>%
  mutate(Rt_exposure=ifelse(id=="'(d) Different time-varying\ntransmission'", Rt_exposure/1.5, Rt_exposure),
         beta_exposure=ifelse(id=="'(d) Different time-varying\ntransmission'", beta_exposure/1.5, beta_exposure))

Rt_df <- exposure_df %>%
  mutate(y=log(Rt_exposure), model="log~R[t]") %>%
  dplyr::select(model, unit, id, lab, y, t)

beta_df <- exposure_df %>%
  mutate(y=log(beta_exposure), model="log~beta[t]") %>%
  dplyr::select(model, unit, id, lab, y, t)
##############################################################################################################################
# merge all model data together for graphing
p.df <- inc_df %>% rbind(loginc_df) %>% rbind(beta_df) %>% rbind(growth_df) %>% rbind(Rt_df) %>% 
  mutate(model = factor(model, levels = c("incidence", "'log incidence'", "'log growth'",
                                          "log~R[t]", "log~beta[t]"))) %>%
  filter(t <= T.total,
         t > (burnin + inf_days*0.5)) %>%
  mutate(t = t - (burnin + inf_days*0.5)) %>%
  group_by(id, model, unit) %>%
  mutate(first.time = ifelse(t==min(t), 1, 0))

p.df.lab <- p.df %>% filter(unit==1 & first.time==1) %>% mutate(checkmark = NA)
p.df.lab$checkmark[p.df.lab$model=="incidence"] <- c("\u2713", rep("X", 4))
p.df.lab$checkmark[p.df.lab$model=="'log incidence'"] <- c(rep("\u2713", 2), rep("X", 3))
p.df.lab$checkmark[p.df.lab$model=="'log growth'"] <- c(rep("\u2713", 3), rep("X", 2))
p.df.lab$checkmark[p.df.lab$model=="log~beta[t]"] <- c(rep("\u2713", 5))
p.df.lab$checkmark[p.df.lab$model=="log~R[t]"] <- c(rep("\u2713", 4), "X")
##############################################################################################################################
Figure1 <- ggplot(p.df, aes(x=t, y=y, linetype=factor(unit), col=factor(unit))) + 
  geom_line() +
  facet_grid(model~id, scales = "free_y", labeller = label_parsed) +
  scale_y_continuous(expand = c(0.2, 0.2), n.breaks = 4) +
  scale_linetype_manual(values = c(1, 5)) +
  geom_text(data = p.df.lab, size = 5,
            aes(label = lab,  x = 0, y = Inf), parse = T, hjust = 0, vjust = 1.3) +
  geom_text(data = p.df.lab, size = 8, aes(label = checkmark,  x = Inf, y = Inf),
            parse = F, hjust = 1, vjust = 1) +
  labs(x = "", y = "") +
  scale_color_brewer(guide = "none", name = "", palette = "Set1") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), legend.position = "none",
        strip.text = element_text(size = 14))
Figure1
