# rm(list=ls())
here::i_am("1b_Summarize/3a_Misspecify_GI_summ.R")
source("./global_options.R")
source("./1a_Scripts/0_Format_Table.R")
source("./1b_Summarize/0_AME_Summary_Helpers.R")

p_out <- read_required_rds("./4_Output/misspecify_GI.rds") %>%
  mutate(model = ifelse(model=="beta_est", "beta", model)) %>%
  add_direct_ame_estimands()

eff.truth <- read_required_rds("./4_Output/SIR_RR.rds") %>%
  filter(model != "Rt_wt", !is.na(eff.true)) %>%
  group_by(trans_prob.base1, trans_prob.base2, eff.multi, model) %>%
  summarise(eff.true = mean(eff.true), .groups = "drop") %>%
  mutate(model = ifelse(model=="Rt_cohort", "Rt_wt", model))
##############################################################################################################################
# Power / type I error rate
power.df.mean <- p_out %>%
  filter(var_spe==1, !is.na(p)) %>%
  group_by(model, eff.multi, mean_spe) %>%
  summarise(p = mean(p < 0.05)*100, nsim = n(), .groups = "drop")

power.df.var <- p_out %>%
  filter(mean_spe==1, !is.na(p)) %>%
  group_by(model, eff.multi, var_spe) %>%
  summarise(p = mean(p < 0.05)*100, nsim = n(), .groups = "drop")
##############################################################################################################################
bias.AME.mean <- p_out %>% # keep variance at true value, vary mean
  filter(var_spe==1) %>%
  group_by(model, eff.multi, mean_spe) %>%
  summarise(nsim = n(),
            Y.obs = mean(Y.obs),
            Y.untrt.true = mean(Y.untrt.true),
            Y.untrt = mean(Y.untrt),
            Y.untrt.adj1 = mean(Y.untrt.adj1),
            Y.untrt.adj2 = mean(Y.untrt.adj2),
            AME.true = mean(AME.true),
            AME = mean(AME),
            AME.adj1 = mean(AME.adj1),
            AME.adj2 = mean(AME.adj2),
            .groups = "drop") %>%
  mutate(bias.fit = AME - AME.true,
         bias.adj1 = AME.adj1 - AME.true,
         bias.adj2 = AME.adj2 - AME.true,
         bias.fit.pct = bias.fit / Y.untrt.true,
         bias.adj1.pct = bias.adj1 / Y.untrt.true,
         bias.adj2.pct = bias.adj2 / Y.untrt.true)

bias.AME.var <- p_out %>% # keep mean at true value, vary variance
  filter(mean_spe==1) %>%
  group_by(model, eff.multi, var_spe) %>%
  summarise(nsim = n(),
            Y.obs = mean(Y.obs),
            Y.untrt.true = mean(Y.untrt.true),
            Y.untrt = mean(Y.untrt),
            Y.untrt.adj1 = mean(Y.untrt.adj1),
            Y.untrt.adj2 = mean(Y.untrt.adj2),
            AME.true = mean(AME.true),
            AME = mean(AME),
            AME.adj1 = mean(AME.adj1),
            AME.adj2 = mean(AME.adj2),
            .groups = "drop") %>%
  mutate(bias.fit = AME - AME.true,
         bias.adj1 = AME.adj1 - AME.true,
         bias.adj2 = AME.adj2 - AME.true,
         bias.fit.pct = bias.fit / Y.untrt.true,
         bias.adj1.pct = bias.adj1 / Y.untrt.true,
         bias.adj2.pct = bias.adj2 / Y.untrt.true)
##############################################################################################################################
# Keep the misspecification level in the grouping.  Averaging over mean_spe or
# var_spe before calculating the reported range would erase the scenario that
# this table is intended to summarize.
bias.original.mean <- p_out %>%
  filter(var_spe==1) %>%
  group_by(eff.multi, mean_spe, model) %>%
  summarise(nsim = n(), eff = mean(effect), .groups = "drop") %>%
  mutate(trans_prob.base1=0.1265, trans_prob.base2=0.115) %>%
  merge(eff.truth, by = c("trans_prob.base1", "trans_prob.base2", "eff.multi", "model")) %>%
  mutate(eff.bias = eff - eff.true,
         eff.bias.pct = eff.bias / eff.true)

bias.original.var <- p_out %>%
  filter(mean_spe==1) %>%
  group_by(eff.multi, var_spe, model) %>%
  summarise(nsim = n(), eff = mean(effect), .groups = "drop") %>%
  mutate(trans_prob.base1=0.1265, trans_prob.base2=0.115) %>%
  merge(eff.truth, by = c("trans_prob.base1", "trans_prob.base2", "eff.multi", "model")) %>%
  mutate(eff.bias = eff - eff.true,
         eff.bias.pct = eff.bias / eff.true)
##############################################################################################################################
# Make kable
# Table 2: Misspecify the mean of SIR generation interval
tbl2.misspecify.mean <- format.tbl(power.df.mean, bias.original.mean, bias.AME.mean) %>% kable(format = "latex")

# Table 2: Misspecify the variance of SIR generation interval
tbl2.misspecify.var <- format.tbl(power.df.var, bias.original.var, bias.AME.var) %>% kable(format = "latex")
