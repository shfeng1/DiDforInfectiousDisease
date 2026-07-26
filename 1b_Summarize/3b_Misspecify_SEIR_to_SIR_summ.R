# rm(list=ls())
here::i_am("1b_Summarize/3b_Misspecify_SEIR_to_SIR_summ.R")
source("./global_options.R")
source("./1a_Scripts/0_Format_Table.R")
source("./1b_Summarize/0_AME_Summary_Helpers.R")

p_out <- read_required_rds("./4_Output/misspecify_SEIR_to_SIR.rds") %>%
  mutate(model.lab = case_when(model=="inc" ~ "incidence",
                               model=="loginc" ~ "log incidence",
                               model=="growth" ~ "log growth",
                               model=="Rt_wt" ~ "Rt (Wallinga Teunis)",
                               model=="Rt_est" ~ "Rt (Prevalence Estimation)",
                               model=="beta" ~ "\u03B2t (Prevalence Estimation)"),
         model.lab = factor(model.lab, levels = c("incidence", "log incidence", "log growth",
                                                  "Rt (Wallinga Teunis)", "Rt (Prevalence Estimation)",
                                                  "\u03B2t (Prevalence Estimation)"))) %>%
  add_direct_ame_estimands()

eff.truth <- read_required_rds("./4_Output/SEIR_RR.rds") %>%
  filter(model != "Rt_wt", !is.na(eff.true),
         near(trans_prob.base1, 0.1265), near(trans_prob.base2, 0.115)) %>%
  group_by(trans_prob.base1, trans_prob.base2, eff.multi, model) %>%
  summarise(eff.true = mean(eff.true), .groups = "drop") %>%
  mutate(model = ifelse(model=="Rt_cohort", "Rt_wt", model))
##############################################################################################################################
# Power / type I error rate
power.df2 <- p_out %>%
  filter(!is.na(p)) %>%
  group_by(model, model.lab, eff.multi) %>%
  summarise(p = mean(p < 0.05)*100, nsim = n(), .groups = "drop")
##############################################################################################################################
bias.AME2 <- p_out %>%
  group_by(model, model.lab, eff.multi) %>%
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
bias.original2 <- p_out %>%
  group_by(eff.multi, model, model.lab) %>%
  summarise(nsim = n(), eff = mean(effect), .groups = "drop") %>%
  merge(eff.truth %>% dplyr::select(eff.multi, model, eff.true),
        by = c("eff.multi", "model")) %>%
  mutate(eff.bias = eff - eff.true,
         eff.bias.pct = eff.bias / eff.true)
##############################################################################################################################
# Make kable
# Table 2: Misspecify SEIR as SIR
tbl2.misspecify.SEIR <- format.tbl(power.df2, bias.original2, bias.AME2) %>% kable(format = "latex")
