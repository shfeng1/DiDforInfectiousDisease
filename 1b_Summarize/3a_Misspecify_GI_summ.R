here::i_am("1b_Summarize/3a_Misspecify_GI_summ.R")
source("./global_options.R")
source("./1a_Scripts/0_Format_Table.R")

p_out <- readRDS("./4_Output/misspecify_GI.rds")

eff.truth <- readRDS("./4_Output/SIR_RR.rds") %>%
  filter(trans_prob.base1=="0.1265", trans_prob.base2=="0.115") %>%
  group_by(trans_prob.base1, trans_prob.base2, eff.multi, model) %>%
  summarise(eff.true = mean(eff.true))
##############################################################################################################################
# Power / type I error rate
power.df.mean <- p_out %>%
  filter(var_spe==1) %>%
  group_by(model, eff.multi, mean_spe) %>%
  summarise(p = mean(p < 0.05)*100)

power.df.var <- p_out %>%
  filter(mean_spe==1) %>%
  group_by(model, eff.multi, var_spe) %>%
  summarise(p = mean(p < 0.05)*100)
##############################################################################################################################
bias.AME.mean <- p_out %>% # keep variance at true value, vary mean
  filter(var_spe==1) %>%
  group_by(model, eff.multi, mean_spe) %>%
  summarise(Y.untrt.fit = mean(Y.trt-AME),
            Y.trt = mean(Y.trt), Y.untrt.true = mean(Y.untrt.true)) %>%
  mutate(bias.fit = Y.untrt.fit - Y.untrt.true,
         bias.fit.pct = bias.fit / Y.untrt.true)

bias.AME.var <- p_out %>% # keep mean at true value, vary variance
  filter(mean_spe==1) %>%
  group_by(model, eff.multi, var_spe) %>%
  summarise(Y.untrt.fit = mean(Y.trt-AME),
            Y.trt = mean(Y.trt), Y.untrt.true = mean(Y.untrt.true)) %>%
  mutate(bias.fit = Y.untrt.fit - Y.untrt.true,
         bias.fit.pct = bias.fit / Y.untrt.true)
##############################################################################################################################
bias.original.mean <- p_out %>% 
  filter(var_spe==1) %>%
  group_by(eff.multi, model) %>%
  summarise(eff = mean(effect)) %>%
  merge(eff.truth, by = c("eff.multi", "model")) %>%
  mutate(eff.bias = eff - eff.true, eff.bias.pct = eff.bias / eff.true)

bias.original.var <- p_out %>% 
  filter(mean_spe==1) %>%
  group_by(eff.multi, model) %>%
  summarise(eff = mean(effect)) %>%
  merge(eff.truth, by = c("eff.multi", "model")) %>%
  mutate(eff.bias = eff - eff.true, eff.bias.pct = eff.bias / eff.true)
##############################################################################################################################
# Make kable
# Table 2: Misspecify the mean of SIR generation interval
tbl2.misspecify.mean <- format.tbl(power.df.mean, bias.original.mean, bias.AME.mean) %>% kable(format = "latex")

# Table 2: Misspecify the variance of SIR generation interval
tbl2.misspecify.var <- format.tbl(power.df.var, bias.original.var, bias.AME.var) %>% kable(format = "latex")
