rm(list=ls())
here::i_am("1b_Summarize/4_Small_N1_summ.R")
source("./global_options.R")
source("./1a_Scripts/0_Format_Table.R")

p_out <- readRDS("4_Output/SIR_N1=5.rds") %>%
  mutate(model.lab = case_when(var=="inc" ~ "incidence", 
                               var=="loginc" ~ "log incidence",
                               var=="growth" ~ "log growth",
                               var=="R_wt" ~ "Rt (Wallinga Teunis)",
                               var=="R_est" ~ "Rt (Prevalence Estimation)",
                               var=="beta_est" ~ "\u03B2t (Prevalence Estimation)"),
         model.lab = factor(model.lab, levels = c("incidence", "log incidence", "log growth",
                                                  "Rt (Wallinga Teunis)", "Rt (Prevalence Estimation)",
                                                  "\u03B2t (Prevalence Estimation)")))
##############################################################################################################################
# Power / type I error rate
power.df2 <- p_out %>%
  group_by(var, model.lab, eff.multi) %>%
  summarise(wild = mean(wild < 0.05)*100, 
            normal = mean(normal < 0.05)*100, 
            nsim = n())
##############################################################################################################################
# Make kable
# Type I error rate
typeIerr <- power.df2 %>% filter(eff.multi==1) %>% group_by(var) %>%
  summarise(error_wild = round(mean(wild)), error_norm = round(mean(normal))) %>% 
  dplyr::select(var, error_wild, error_norm)

# Powers at selected effect size
power90 <- power.df2 %>% filter(eff.multi==0.90) %>% group_by(var) %>%
  summarise(wild90 = round(mean(wild)), norm90 = round(mean(normal))) %>% 
  dplyr::select(var, wild90, norm90)

power95 <- power.df2 %>% filter(eff.multi==0.95) %>% group_by(var) %>%
  summarise(wild95 = round(mean(wild)), norm95 = round(mean(normal))) %>% 
  dplyr::select(var, wild95, norm95)

power105 <- power.df2 %>% filter(eff.multi==1.05) %>% group_by(var) %>%
  summarise(wild105 = round(mean(wild)), norm105 = round(mean(normal))) %>% 
  dplyr::select(var, wild105, norm105)

power110 <- power.df2 %>% filter(eff.multi==1.1) %>% group_by(var) %>%
  summarise(wild110 = round(mean(wild)), norm110 = round(mean(normal))) %>% 
  dplyr::select(var, wild110, norm110)

power.out <- typeIerr %>% merge(power90) %>% merge(power95) %>% merge(power105) %>% merge(power110)

# Put in Latex format for wild score bootstrap
power.out %>% dplyr::select(var, error_wild, wild90, wild95, wild105, wild110) %>%
  arrange(match(var, c("inc", "loginc", "growth", "R_wt", "R_est", "beta_est"))) %>% 
  kable(format = "latex")

# Put in Latex format for normal-based SE
power.out %>% dplyr::select(var, error_norm, norm90, norm95, norm105, norm110) %>%
  arrange(match(var, c("inc", "loginc", "growth", "R_wt", "R_est", "beta_est"))) %>% 
  kable(format = "latex")
