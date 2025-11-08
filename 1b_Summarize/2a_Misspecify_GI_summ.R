rm(list=ls())
here::i_am("1b_Summarize/2a_Misspecify_GI_summ.R")
source("./global_options.R")
source("./1a_Scripts/0_Format_Table.R")

p_out <- readRDS("./4_Output/misspecify_GI.rds") %>%
  mutate(model = ifelse(model=="beta_est", "beta", model))

eff.truth <- readRDS("./4_Output/SIR_RR.rds") %>%
  filter(model != "Rt_wt", !is.na(eff.true)) %>%
  group_by(trans_prob.base1, trans_prob.base2, eff.multi, model) %>%
  summarise(eff.true = mean(eff.true)) %>%
  mutate(model = ifelse(model=="Rt_cohort", "Rt_wt", model))
##############################################################################################################################
# Power / type I error rate
power.df.mean <- p_out %>%
  filter(var_spe==1) %>%
  group_by(model, eff.multi, mean_spe) %>%
  summarise(p = mean(p < 0.05)*100, nsim = n())

power.df.var <- p_out %>%
  filter(mean_spe==1) %>%
  group_by(model, eff.multi, var_spe) %>%
  summarise(p = mean(p < 0.05)*100, nsim = n())
##############################################################################################################################
bias.AME.mean <- p_out %>% # keep variance at true value, vary mean
  filter(var_spe==1) %>%
  group_by(model, eff.multi, mean_spe) %>%
  summarise(nsim = n(), Y.trt = mean(Y.trt), Y.untrt.true = mean(Y.untrt.true),
            Y.untrt.fit = mean(Y.untrt), Y.untrt.adj = mean(Y.untrt.adj2)) %>%
  mutate(bias.fit = Y.untrt.fit - Y.untrt.true,
         bias.adj = Y.untrt.adj - Y.untrt.true,
         bias.adj.pct = bias.adj / Y.untrt.true)

bias.AME.var <- p_out %>% # keep mean at true value, vary variance
  filter(mean_spe==1) %>%
  group_by(model, eff.multi, var_spe) %>%
  summarise(nsim = n(), Y.trt = mean(Y.trt), Y.untrt.true = mean(Y.untrt.true),
            Y.untrt.fit = mean(Y.untrt), Y.untrt.adj = mean(Y.untrt.adj2)) %>%
  mutate(bias.fit = Y.untrt.fit - Y.untrt.true,
         bias.adj = Y.untrt.adj - Y.untrt.true,
         bias.adj.pct = bias.adj / Y.untrt.true)
##############################################################################################################################
bias.original.mean <- p_out %>% 
  filter(var_spe==1) %>%
  group_by(eff.multi, model) %>%
  summarise(nsim = n(), eff = mean(effect)) %>%
  mutate(trans_prob.base1="0.1265", trans_prob.base2="0.115") %>%
  merge(eff.truth, by = c("trans_prob.base1", "trans_prob.base2", "eff.multi", "model")) %>%
  mutate(eff.bias = eff - eff.true, eff.bias.pct = eff.bias / eff.true)

bias.original.var <- p_out %>% 
  filter(mean_spe==1) %>%
  group_by(eff.multi, model) %>%
  summarise(nsim = n(), eff = mean(effect)) %>%
  mutate(trans_prob.base1="0.1265", trans_prob.base2="0.115") %>%
  merge(eff.truth, by = c("trans_prob.base1", "trans_prob.base2", "eff.multi", "model")) %>%
  mutate(eff.bias = eff - eff.true, eff.bias.pct = eff.bias / eff.true)
##############################################################################################################################
# Make kable
format.tbl(power.df.mean, bias.original.mean, bias.AME.mean) %>% kable(format = "latex")
format.tbl(power.df.var, bias.original.var, bias.AME.var) %>% kable(format = "latex")
