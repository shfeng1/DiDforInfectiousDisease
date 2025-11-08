rm(list=ls())
here::i_am("1b_Summarize/2b_Misspecify_SEIR_to_SIR_summ.R")
source("./global_options.R")
source("./1a_Scripts/0_Format_Table.R")

p_out <- readRDS("./4_Output/misspecify_SEIR_to_SIR_inf_days.rds") %>%
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
  filter(model != "Rt_wt", !is.na(eff.true), trans_prob.base1=="0.1265", trans_prob.base2=="0.115") %>%
  group_by(trans_prob.base1, trans_prob.base2, eff.multi, model) %>%
  summarise(eff.true = mean(eff.true)) %>%
  mutate(model = ifelse(model=="Rt_cohort", "Rt_wt", model))
##############################################################################################################################
# Power / type I error rate
power.df2 <- p_out %>%
  group_by(model, model.lab, eff.multi) %>%
  summarise(p = mean(p < 0.05)*100, nsim = n())
##############################################################################################################################
bias.df <- p_out %>%
  data.frame() %>%
  group_by(eff.multi) %>%
  mutate(AME.true = Y.trt - Y.untrt.true,
         AME.fit = Y.trt - Y.untrt,
         AME.adj = Y.trt - Y.untrt.adj2,
         eff.true = ifelse(model=="inc", AME.true, eff.multi))
##############################################################################################################################
bias.AME2 <- bias.df %>%
  group_by(model, model.lab, eff.multi) %>%
  summarise(nsim = n(), Y.trt = mean(Y.trt), Y.untrt.fit = mean(Y.untrt), Y.untrt.true = mean(Y.untrt.true), 
            AME.true = mean(AME.true), AME.fit = mean(AME.fit), AME.adj = mean(AME.adj, na.rm = T)) %>%
  mutate(bias.fit = AME.fit - AME.true,
         bias.adj = ifelse(model %in% c("inc", "loginc"), bias.fit, AME.adj - AME.true),
         bias.fit.pct = bias.fit / Y.untrt.true,
         bias.adj.pct = bias.adj / Y.untrt.true)
##############################################################################################################################
bias.original2 <- bias.df %>% 
  group_by(eff.multi, model, model.lab) %>%
  summarise(nsim = n(), eff = mean(effect)) %>%
  merge(eff.truth, by = c("eff.multi", "model"))
bias.original2$eff.true[bias.original2$model=="inc"] <- bias.AME2$AME.true[bias.AME2$model=="inc"]
bias.original2$eff.bias <- bias.original2$eff - bias.original2$eff.true
bias.original2$eff.bias.pct <- bias.original2$eff.bias / bias.original2$eff.true
##############################################################################################################################
# Make kable
format.tbl(power.df2, bias.original2, bias.AME2) %>% kable(format = "latex")
