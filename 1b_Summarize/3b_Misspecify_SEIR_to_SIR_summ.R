here::i_am("1b_Summarize/3b_Misspecify_SEIR_to_SIR_summ.R")
source("./global_options.R")
source("./1a_Scripts/0_Format_Table.R")

p_out <- readRDS("./4_Output/misspecify_SEIR_to_SIR.rds") %>%
  filter(model %in% model.list)

eff.truth <- readRDS("./4_Output/SEIR_RR.rds") %>%
  filter(trans_prob.base1=="0.1265", trans_prob.base2=="0.115") %>%
  group_by(trans_prob.base1, trans_prob.base2, eff.multi, model) %>%
  summarise(eff.true = mean(eff.true))
##############################################################################################################################
# Power / type I error rate
power.df2 <- p_out %>%
  group_by(model, eff.multi) %>%
  summarise(p = mean(p < 0.05)*100, nsim = n())
##############################################################################################################################
bias.df <- p_out %>%
  data.frame() %>%
  group_by(eff.multi) %>%
  mutate(AME.true = Y.trt - Y.untrt.true,
         eff.true = ifelse(model=="inc", AME.true, eff.multi))
##############################################################################################################################
bias.AME2 <- bias.df %>%
  group_by(model, eff.multi) %>%
  summarise(Y.untrt.true = mean(Y.untrt.true),
            AME.true = mean(AME.true), AME.fit = mean(AME)) %>%
  mutate(bias.fit = AME.fit - AME.true,
         bias.fit.pct = bias.fit / Y.untrt.true)
##############################################################################################################################
bias.original2 <- bias.df %>% 
  group_by(eff.multi, model) %>%
  summarise(eff = mean(effect)) %>%
  merge(eff.truth, by = c("eff.multi", "model"))
bias.original2$eff.true[bias.original2$model=="inc"] <- bias.AME2$AME.true[bias.AME2$model=="inc"]
bias.original2$eff.bias <- bias.original2$eff - bias.original2$eff.true
bias.original2$eff.bias.pct <- bias.original2$eff.bias / bias.original2$eff.true
##############################################################################################################################
# Make kable
# Table 2: Misspecify SEIR as SIR
tbl2.misspecify.SEIR <- format.tbl(power.df2, bias.original2, bias.AME2) %>% kable(format = "latex")
