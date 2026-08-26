# rm(list=ls())
here::i_am("1b_Summarize/2b_SEIR_summ.R")
source("./global_options.R")
source("./1a_Scripts/0_Format_Table.R")

p_out <- readRDS("4_Output/SEIR_base_case.rds")
eff.truth <- readRDS("./4_Output/SEIR_RR.rds") %>%
  filter(model != "Rt_wt", !is.na(eff.true)) %>%
  group_by(trans_prob.base1, trans_prob.base2, eff.multi, model) %>%
  summarise(eff.true = mean(eff.true)) %>%
  mutate(model = ifelse(model=="Rt_cohort", "Rt_wt", model))
##############################################################################################################################
# Power / type I error rate
power.df2 <- p_out %>%
  filter(!is.na(p), trans_prob.base1=="0.1265", trans_prob.base2=="0.115") %>%
  group_by(pop.size, trans_prob.base1, trans_prob.base2, model, eff.multi) %>%
  summarise(p = mean(p < 0.05)*100, nsim = n())
##############################################################################################################################
bias.df <- p_out %>%
  data.frame() %>%
  mutate(trans_prob.base1=as.character(trans_prob.base1), trans_prob.base2=as.character(trans_prob.base2)) %>%
  group_by(pop.size, trans_prob.base1, trans_prob.base2, eff.multi) %>%
  mutate(AME.true = Y.trt - Y.untrt.true, AME.fit = AME, AME.adj = AME.adj2,
         eff.true = ifelse(model=="inc", AME.true, eff.multi))
##############################################################################################################################
bias.AME2 <- bias.df %>%
  filter(trans_prob.base1=="0.1265", trans_prob.base2=="0.115") %>%
  # filter(eff.multi != 1) %>%
  group_by(pop.size, trans_prob.base1, trans_prob.base2, model, eff.multi) %>%
  summarise(nsim = n(), Y.untrt.true = mean(Y.untrt.true), AME.true = mean(AME.true), 
            AME.fit = mean(AME.fit), AME.adj = mean(AME.adj, na.rm = T)) %>%
  mutate(bias.fit = AME.fit - AME.true,
         bias.adj = ifelse(model %in% c("inc", "loginc"), bias.fit, AME.adj - AME.true),
         bias.fit.pct = bias.fit / Y.untrt.true,
         bias.adj.pct = bias.adj / Y.untrt.true)
##############################################################################################################################
bias.original2 <- bias.df %>% 
  filter(trans_prob.base1=="0.1265", trans_prob.base2=="0.115") %>%
  group_by(trans_prob.base1, trans_prob.base2, eff.multi, model) %>%
  summarise(nsim = n(), eff = mean(effect)) %>%
  merge(eff.truth, by = c("trans_prob.base1", "trans_prob.base2", "eff.multi", "model"))
bias.original2$eff.true[bias.original2$model=="inc"] <- bias.AME2$AME.true[bias.AME2$model=="inc"]
bias.original2$eff.bias <- bias.original2$eff - bias.original2$eff.true
bias.original2$eff.bias.pct <- bias.original2$eff.bias / bias.original2$eff.true
##############################################################################################################################
# Make kable
# Table 2: SEIR with unequal transmission between groups
tbl2.SEIR <- format.tbl(power.df2, bias.original2, bias.AME2) %>% kable(format = "latex")
