rm(list=ls())
here::i_am("1b_Summarize/3_SEIR_long_time_summ.R")
source("./global_options.R")
source("./1a_Scripts/0_Format_Table.R")

p_out <- readRDS("4_Output/SEIR_long_time.rds") %>%
  mutate(model.lab = case_when(model=="inc" ~ "incidence", 
                               model=="loginc" ~ "log incidence",
                               model=="growth" ~ "log growth",
                               model=="Rt_wt" ~ "Rt (Wallinga Teunis)",
                               model=="Rt_est" ~ "Rt (Prevalence Estimation)",
                               model=="beta" ~ "\u03B2t (Prevalence Estimation)"),
         model.lab = factor(model.lab, levels = c("incidence", "log incidence", "log growth",
                                                  "Rt (Wallinga Teunis)", "Rt (Prevalence Estimation)",
                                                  "\u03B2t (Prevalence Estimation)")))

eff.truth <- readRDS("./4_Output/SEIR_RR.rds") %>%
  filter(model != "Rt_wt", !is.na(eff.true)) %>%
  group_by(trans_prob.base1, trans_prob.base2, eff.multi, model) %>%
  summarise(eff.true = mean(eff.true)) %>%
  mutate(model = ifelse(model=="Rt_cohort", "Rt_wt", model))
##############################################################################################################################
# Power / type I error rate
power.df2 <- p_out %>%
  group_by(pop.size, trans_prob.base1, trans_prob.base2, model, model.lab, eff.multi) %>%
  summarise(p = mean(p < 0.05)*100, nsim = n())
##############################################################################################################################
bias.df <- p_out %>%
  filter(eff.multi %in% c(0.9, 0.95)) %>%  data.frame() %>%
  mutate(trans_prob.base1=as.character(trans_prob.base1), trans_prob.base2=as.character(trans_prob.base2)) %>%
  group_by(pop.size, trans_prob.base1, trans_prob.base2, eff.multi) %>%
  mutate(AME.true = Y.trt - Y.untrt.true,
         AME.fit = Y.trt - Y.untrt,
         AME.adj = Y.trt - Y.untrt.adj2,
         eff.true = ifelse(model=="inc", AME.true, eff.multi))
##############################################################################################################################
bias.AME2 <- bias.df %>%
  group_by(pop.size, trans_prob.base1, trans_prob.base2, model, model.lab, eff.multi) %>%
  summarise(nsim = n(), Y.trt = mean(Y.trt), Y.untrt = mean(Y.untrt), Y.untrt.true = mean(Y.untrt.true), 
            AME.true = mean(AME.true), AME.fit = mean(AME.fit), AME.adj = mean(AME.adj, na.rm = T)) %>%
  mutate(bias.fit = abs(AME.fit - AME.true),
         bias.adj = ifelse(model %in% c("inc", "loginc"), bias.fit, AME.adj - AME.true),
         bias.adj = abs(bias.adj),
         bias.fit.pct = bias.fit / Y.untrt.true,
         bias.adj.pct = bias.adj / Y.untrt.true)
##############################################################################################################################
bias.original2 <- bias.df %>% 
  group_by(trans_prob.base1, trans_prob.base2, eff.multi, model, model.lab) %>%
  summarise(nsim = n(), eff = mean(effect)) %>%
  merge(eff.truth, by = c("trans_prob.base1", "trans_prob.base2", "eff.multi", "model"))
bias.original2$eff.true[bias.original2$model=="inc"] <- bias.AME2$AME.true[bias.AME2$model=="inc"]
bias.original2$eff.bias <- bias.original2$eff - bias.original2$eff.true
bias.original2$eff.bias.pct <- bias.original2$eff.bias / bias.original2$eff.true
##############################################################################################################################
# Make kable
power <- power.df2 %>% filter(eff.multi==1) %>% group_by(model) %>%
  summarise(p = round(mean(p), 1)) %>% dplyr::select(model, p) %>%
  arrange(match(model, c("inc", "loginc", "growth", "Rt_wt", "Rt_est", "beta")))
bias.original <- bias.original2 %>% group_by(model) %>% 
  summarise(mean = format(round(mean(abs(eff.bias)), 2), nsmall=2), 
            min = format(round(min(eff.bias), 2), nsmall=2), 
            max = format(round(max(eff.bias), 2), nsmall=2)) %>%
  mutate(bias_original = paste0(mean, " (", min, ", ", max, ")")) %>%
  arrange(match(model, c("inc", "loginc", "growth", "Rt_wt", "Rt_est", "beta")))
bias.AME <- bias.AME2 %>% group_by(model) %>% 
  filter(model != "true") %>%
  summarise(mean = format(round(mean(abs(bias.fit)), 1), nsmall=1), 
            min = format(round(min(bias.fit), 1), nsmall=1), 
            max = format(round(max(bias.fit), 1), nsmall=1)) %>%
  mutate(bias_AME = paste0(mean, " (", min, ", ", max, ")")) %>%
  arrange(match(model, c("inc", "loginc", "growth", "Rt_wt", "Rt_est", "beta")))
bias.AME.adj <- bias.AME2 %>% group_by(model) %>% 
  filter(model != "true") %>%
  summarise(mean = format(round(mean(abs(bias.adj)), 1), nsmall=1), 
            min = format(round(min(bias.adj), 1), nsmall=1), 
            max = format(round(max(bias.adj), 1), nsmall=1)) %>%
  mutate(bias_AME_correct = paste0(mean, " (", min, ", ", max, ")")) %>%
  arrange(match(model, c("inc", "loginc", "growth", "Rt_wt", "Rt_est", "beta")))

power %>% cbind(bias.original$bias_original) %>% cbind(bias.AME$bias_AME) %>% cbind(bias.AME.adj$bias_AME_correct) %>% kable(format = "latex")

