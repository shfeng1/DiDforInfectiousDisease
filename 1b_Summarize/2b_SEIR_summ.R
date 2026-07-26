# rm(list=ls())
here::i_am("1b_Summarize/2b_SEIR_summ.R")
source("./global_options.R")
source("./1a_Scripts/0_Format_Table.R")
source("./1b_Summarize/0_AME_Summary_Helpers.R")

p_out <- read_required_rds("./4_Output/SEIR_base_case.rds") %>%
  mutate(model.lab = case_when(model=="inc" ~ "incidence",
                               model=="loginc" ~ "log incidence",
                               model=="growth" ~ "log growth",
                               model=="Rt_wt" ~ "Rt (Wallinga Teunis)",
                               model=="Rt_est" ~ "Rt (Instantaneous Estimation)",
                               model=="beta" ~ "\u03B2t (Instantaneous Estimation)"),
         model.lab = factor(model.lab, levels = c("incidence", "log incidence", "log growth",
                                                  "Rt (Wallinga Teunis)", "Rt (Instantaneous Estimation)",
                                                  "\u03B2t (Instantaneous Estimation)")))

eff.truth <- read_required_rds("./4_Output/SEIR_RR.rds") %>%
  filter(model != "Rt_wt", !is.na(eff.true)) %>%
  group_by(trans_prob.base1, trans_prob.base2, eff.multi, model) %>%
  summarise(eff.true = mean(eff.true), .groups = "drop") %>%
  mutate(model = ifelse(model=="Rt_cohort", "Rt_wt", model))
##############################################################################################################################
# Power / type I error rate.  Keep transmission parameters numeric throughout;
# character comparisons can fail after RDS serialization or rounding changes.
power.df2 <- p_out %>%
  filter(!is.na(p), near(trans_prob.base1, 0.1265), near(trans_prob.base2, 0.115)) %>%
  group_by(pop.size, trans_prob.base1, trans_prob.base2, model, model.lab, eff.multi) %>%
  summarise(p = mean(p < 0.05)*100, nsim = n(), .groups = "drop")
##############################################################################################################################
bias.df <- p_out %>%
  add_direct_ame_estimands() %>%
  mutate(eff.true = ifelse(model=="inc", AME.true, eff.multi))
##############################################################################################################################
bias.AME2 <- bias.df %>%
  filter(near(trans_prob.base1, 0.1265), near(trans_prob.base2, 0.115)) %>%
  group_by(pop.size, trans_prob.base1, trans_prob.base2, model, model.lab, eff.multi) %>%
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
bias.original2 <- bias.df %>%
  filter(near(trans_prob.base1, 0.1265), near(trans_prob.base2, 0.115)) %>%
  group_by(trans_prob.base1, trans_prob.base2, eff.multi, model, model.lab) %>%
  summarise(nsim = n(), eff = mean(effect), .groups = "drop") %>%
  merge(eff.truth, by = c("trans_prob.base1", "trans_prob.base2", "eff.multi", "model")) %>%
  mutate(eff.bias = eff - eff.true,
         eff.bias.pct = eff.bias / eff.true)
##############################################################################################################################
p.power2 <- ggplot(power.df2 %>% filter(! model %in% c("true"),
                                        (model =="beta") | (eff.multi==1)),
                   aes(eff.multi, p, col=model.lab)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  geom_hline(yintercept = 5, lty = "dashed") +
  scale_y_continuous(limit=c(0, 100)) +
  scale_color_manual(name = "", values = pal) +
  labs(x = "", y = expression("P(reject"~H[0]*")")) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

p.bias.eff2 <- ggplot(bias.original2 %>% filter(model=="beta"), aes(eff.multi, eff.bias.pct*100)) +
  geom_point(col=pal[length(pal)]) +
  geom_line(col=pal[length(pal)]) +
  theme_bw() +
  labs(x = "", y = "Bias in RR (%)") +
  scale_y_continuous(limit=c(-40, 40), breaks = seq(-40, 40, 20)) +
  scale_color_manual(name = "", values = pal) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

p.bias.AME2 <- ggplot(bias.AME2 %>% filter(model=="beta"), aes(eff.multi, bias.adj2.pct*100)) +
  geom_point(col=pal[length(pal)]) +
  geom_line(col=pal[length(pal)]) +
  theme_bw() +
  labs(x = "", y = "Bias in AME as % of post-\nintervention case burden") +
  scale_y_continuous(limit=c(-40, 40), breaks = seq(-40, 40, 20)) +
  scale_color_manual(name = "", values = pal) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

p2 <- ggarrange(p.power2, p.bias.eff2, p.bias.AME2, ncol = 3, nrow = 1, common.legend = TRUE, legend = "bottom")
p2 <- annotate_figure(p2, top = text_grob(expression((b)~beta["0,t"]~"="~0.115*","~beta["1,t"]~"="~0.1265), size = 14, hjust = 0, x = 0.02))
p2
##############################################################################################################################
# Make kable
# Table 2: SEIR with unequal transmission between groups
tbl2.SEIR <- format.tbl(power.df2, bias.original2, bias.AME2) %>% kable(format = "latex")
