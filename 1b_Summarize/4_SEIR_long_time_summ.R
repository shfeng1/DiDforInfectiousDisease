rm(list=ls())
here::i_am("1b_Summarize/4_SEIR_long_time_summ.R")
source("./global_options.R")
source("./1a_Scripts/0_Format_Table.R")

p_out <- readRDS("4_Output/SEIR_long_time.rds") %>% filter(model=="beta_exposure")

eff.truth <- readRDS("./4_Output/SEIR_RR.rds") %>%
  filter(trans_prob.base1=="0.1265", trans_prob.base2=="0.115") %>%
  group_by(trans_prob.base1, trans_prob.base2, eff.multi, model) %>%
  summarise(eff.true = mean(eff.true))
##############################################################################################################################
bias.df <- p_out %>%
  filter(eff.multi %in% c(0.9, 0.95)) %>% data.frame() %>%
  group_by(pop.size, trans_prob.base1, trans_prob.base2, eff.multi) %>%
  mutate(AME.true = Y.trt - Y.untrt.true, nsim = n())
##############################################################################################################################
bias.AME2 <- bias.df %>%
  group_by(pop.size, trans_prob.base1, trans_prob.base2, model, eff.multi) %>%
  summarise(AME.true = mean(AME.true), AME.fit = mean(AME),
            AME.adj = mean(AME.adj2)) %>%
  mutate(bias.fit = AME.fit - AME.true,
         bias.adj = AME.adj - AME.true)
##############################################################################################################################
bias.original2 <- bias.df %>% 
  group_by(trans_prob.base1, trans_prob.base2, eff.multi, model) %>%
  summarise(eff = mean(effect)) %>%
  merge(eff.truth, by = c("eff.multi", "model"))
bias.original2$eff.bias <- bias.original2$eff - bias.original2$eff.true
##############################################################################################################################
tbl.A3 <- format.tbl.A3(bias.original2, bias.AME2) %>%
  mutate(model="Log $\\beta_t$") %>%
  kable(format="latex", booktabs=TRUE, escape=FALSE,
        col.names=c(
          "Outcome specification",
          linebreak("Bias on the original scale: mean absolute value (range)"),
          linebreak("Bias on the AME scale before correction: mean absolute value (range)"),
          linebreak("Bias on the AME scale after correction: mean absolute value (range)")
        ))
tbl.A3
