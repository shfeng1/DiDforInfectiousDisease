here::i_am("1b_Summarize/2a_SIR_summ.R")
source("./global_options.R")
source("./1a_Scripts/0_Format_Table.R")

p_out <- readRDS("4_Output/SIR_base_case.rds") %>%
  filter(model %in% model.list) %>%
  mutate(model.lab=factor(unname(model.labels[model]), levels=unname(model.labels)))

eff.truth <- readRDS("./4_Output/SIR_RR.rds") %>%
  group_by(trans_prob.base1, trans_prob.base2, eff.multi, model) %>%
  summarise(eff.true = mean(eff.true), nsim = n())
##############################################################################################################################
# Power / type I error rate
power.df1 <- p_out %>%
  filter(!is.na(p), trans_prob.base1 == trans_prob.base2) %>%
  group_by(pop.size, trans_prob.base1, trans_prob.base2, model, model.lab, eff.multi) %>%
  summarise(p = mean(p < 0.05)*100, nsim = n())

power.df2 <- p_out %>%
  filter(!is.na(p), trans_prob.base1 > trans_prob.base2) %>%
  group_by(pop.size, trans_prob.base1, trans_prob.base2, model, model.lab, eff.multi) %>%
  summarise(p = mean(p < 0.05)*100, nsim = n())
##############################################################################################################################
bias.df <- p_out %>%
  data.frame() %>%
  mutate(trans_prob.base1=as.character(trans_prob.base1), trans_prob.base2=as.character(trans_prob.base2)) %>%
  group_by(pop.size, trans_prob.base1, trans_prob.base2, eff.multi) %>%
  mutate(AME.true = Y.trt - Y.untrt.true, AME.fit = AME,
         eff.true = ifelse(model=="inc", AME.true, eff.multi))
##############################################################################################################################
bias.AME1 <- bias.df %>%
  filter(trans_prob.base1 == trans_prob.base2) %>%
  group_by(pop.size, trans_prob.base1, trans_prob.base2, model, model.lab, eff.multi) %>%
  summarise(nsim = n(), Y.untrt.true = mean(Y.untrt.true), AME.true = mean(AME.true), 
            AME.fit = mean(AME.fit)) %>%
  mutate(bias.fit = AME.fit - AME.true,
         bias.fit.pct = bias.fit / Y.untrt.true)

bias.AME2 <- bias.df %>%
  filter(trans_prob.base1 > trans_prob.base2) %>%
  group_by(pop.size, trans_prob.base1, trans_prob.base2, model, model.lab, eff.multi) %>%
  summarise(nsim = n(), Y.untrt.true = mean(Y.untrt.true), AME.true = mean(AME.true), 
            AME.fit = mean(AME.fit)) %>%
  mutate(bias.fit = AME.fit - AME.true,
         bias.fit.pct = bias.fit / Y.untrt.true)
##############################################################################################################################
bias.original1 <- bias.df %>% 
  filter(trans_prob.base1 == trans_prob.base2) %>%
  group_by(trans_prob.base1, trans_prob.base2, eff.multi, model, model.lab) %>%
  summarise(nsim = n(), eff = mean(effect)) %>%
  merge(eff.truth, by = c("trans_prob.base1", "trans_prob.base2", "eff.multi", "model"))
bias.original1$eff.true[bias.original1$model=="inc"] <- bias.AME1$AME.true[bias.AME1$model=="inc"]
bias.original1$eff.bias <- bias.original1$eff - bias.original1$eff.true
bias.original1$eff.bias.pct <- bias.original1$eff.bias / bias.original1$eff.true

bias.original2 <- bias.df %>% 
  filter(trans_prob.base1 > trans_prob.base2) %>%
  group_by(trans_prob.base1, trans_prob.base2, eff.multi, model, model.lab) %>%
  summarise(nsim = n(), eff = mean(effect)) %>%
  merge(eff.truth, by = c("trans_prob.base1", "trans_prob.base2", "eff.multi", "model"))
bias.original2$eff.true[bias.original2$model=="inc"] <- bias.AME2$AME.true[bias.AME2$model=="inc"]
bias.original2$eff.bias <- bias.original2$eff - bias.original2$eff.true
bias.original2$eff.bias.pct <- bias.original2$eff.bias / bias.original2$eff.true
##############################################################################################################################
p.power1 <- ggplot(power.df1, aes(eff.multi, p, col=model.lab)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  geom_hline(yintercept = 5, lty = "dashed") +
  scale_y_continuous(limit=c(0, 100)) +
  scale_color_manual(name = "", values = model.colors) +
  labs(x = "", y = expression("P(reject"~H[0]*")")) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

p.power2 <- ggplot(power.df2 %>% filter(model=="beta_exposure" | eff.multi==1),
                   aes(eff.multi, p, col=model.lab)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  geom_hline(yintercept = 5, lty = "dashed") +
  scale_y_continuous(limit=c(0, 100)) +
  scale_color_manual(name = "", values = model.colors) +
  labs(x = "", y = expression("P(reject"~H[0]*")")) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

p.bias.eff1 <- ggplot(bias.original1 %>% filter(model!="inc"),
                      aes(eff.multi, eff.bias.pct*100, col=model.lab)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(x = "", y = "Bias in RR (%)") +
  scale_y_continuous(limit=c(-50, 50), breaks = seq(-50, 50, 25)) +
  scale_color_manual(name = "", values = model.colors) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

p.bias.eff2 <- ggplot(bias.original2 %>% filter(model=="beta_exposure"), aes(eff.multi, eff.bias.pct*100)) +
  geom_point(col=model.colors[["\u03B2t"]]) +
  geom_line(col=model.colors[["\u03B2t"]]) +
  theme_bw() +
  labs(x = "", y = "Bias in RR (%)") +
  scale_y_continuous(limit=c(-50, 50), breaks = seq(-50, 50, 25)) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

# benchmarking against what % of the total case (infection) burden
p.bias.AME1 <- ggplot(bias.AME1, aes(eff.multi, bias.fit.pct*100, col=model.lab)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(x = "", y = "Bias in AME as % of post-\nintervention case burden") +
  scale_y_continuous(limit=c(-50, 50), breaks = seq(-50, 50, 25)) +
  scale_color_manual(name = "", values = model.colors) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

p.bias.AME2 <- ggplot(bias.AME2 %>% filter(model=="beta_exposure"), aes(eff.multi, bias.fit.pct*100)) +
  geom_point(col=model.colors[["\u03B2t"]]) +
  geom_line(col=model.colors[["\u03B2t"]]) +
  theme_bw() +
  labs(x = "", y = "Bias in AME as % of post-\nintervention case burden") +
  scale_y_continuous(limit=c(-50, 50), breaks = seq(-50, 50, 25)) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

p1 <- ggarrange(p.power1, p.bias.eff1, p.bias.AME1, ncol = 3, nrow = 1, common.legend = TRUE, legend = "none")
p1 <- annotate_figure(p1, top = text_grob(expression((a)~beta["0,t"]~"="~beta["1,t"]~"="~0.115), size = 14, hjust = 0, x = 0.02))

p2 <- ggarrange(p.power2, p.bias.eff2, p.bias.AME2, ncol = 3, nrow = 1, common.legend = TRUE, legend = "bottom")
p2 <- annotate_figure(p2, top = text_grob(expression((b)~beta["0,t"]~"="~0.115*","~beta["1,t"]~"="~0.1265), size = 14, hjust = 0, x = 0.02))

Figure2 <- ggarrange(p1, p2, ncol = 1, nrow = 2, heights = c(0.88, 1))
Figure2
##############################################################################################################################
# Make kable
# Table 2: SIR with equal transmission between groups
tbl2.SIR1 <- format.tbl(power.df1, bias.original1, bias.AME1) %>% kable(format = "latex")

# Table 2: SIR with unequal transmission between groups
tbl2.SIR2 <- format.tbl(power.df2, bias.original2, bias.AME2) %>% kable(format = "latex")
