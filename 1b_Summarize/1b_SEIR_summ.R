rm(list=ls())
here::i_am("1b_Summarize/1b_SEIR_summ.R")
source("./global_options.R")
source("./1a_Scripts/0_Format_Table.R")

p_out <- readRDS("4_Output/SEIR_base_case.rds") %>%
  mutate(model.lab = case_when(model=="inc" ~ "incidence", 
                               model=="loginc" ~ "log incidence",
                               model=="growth" ~ "log growth",
                               model=="Rt" ~ "Rt (Wallinga Teunis)",
                               model=="beta" ~ "\u03B2t (Prevalence Estimation)"),
         model.lab = factor(model.lab, levels = c("incidence", "log incidence", "log growth",
                                                  "Rt (Wallinga Teunis)", "\u03B2t (Prevalence Estimation)")))

eff.truth <- readRDS("./4_Output/SEIR_RR.rds") %>%
  filter(model != "Rt_wt", !is.na(eff.true)) %>%
  group_by(trans_prob.base1, trans_prob.base2, eff.multi, model) %>%
  summarise(eff.true = mean(eff.true)) %>%
  mutate(model = ifelse(model=="Rt_cohort", "Rt_wt", model))
##############################################################################################################################
# Power / type I error rate
power.df1 <- p_out %>%
  filter(!is.na(p), trans_prob.base1=="0.115", trans_prob.base2=="0.115") %>%
  group_by(pop.size, trans_prob.base1, trans_prob.base2, model, model.lab, eff.multi) %>%
  summarise(p = mean(p < 0.05)*100, nsim = n())

power.df2 <- p_out %>%
  filter(!is.na(p), trans_prob.base1=="0.1265", trans_prob.base2=="0.115") %>%
  group_by(pop.size, trans_prob.base1, trans_prob.base2, model, model.lab, eff.multi) %>%
  summarise(p = mean(p < 0.05)*100, nsim = n())
##############################################################################################################################
bias.df <- p_out %>%
  data.frame() %>%
  mutate(trans_prob.base1=as.character(trans_prob.base1), trans_prob.base2=as.character(trans_prob.base2)) %>%
  group_by(pop.size, trans_prob.base1, trans_prob.base2, eff.multi) %>%
  mutate(AME.true = Y.trt - Y.untrt.true,
         AME.fit = Y.trt - Y.untrt,
         AME.adj = Y.trt - Y.untrt.adj2,
         eff.true = ifelse(model=="inc", AME.true, eff.multi))
##############################################################################################################################
bias.AME1 <- bias.df %>%
  filter(trans_prob.base1=="0.115", trans_prob.base2=="0.115") %>%
  # filter(eff.multi != 1) %>%
  group_by(pop.size, trans_prob.base1, trans_prob.base2, model, model.lab, eff.multi) %>%
  summarise(nsim = n(), Y.trt = mean(Y.trt), Y.untrt = mean(Y.untrt), Y.untrt.true = mean(Y.untrt.true), 
            AME.true = mean(AME.true), AME.fit = mean(AME.fit), AME.adj = mean(AME.adj, na.rm = T)) %>%
  mutate(bias.fit = AME.fit - AME.true,
         bias.adj = ifelse(model %in% c("inc", "loginc"), bias.fit, AME.adj - AME.true),
         bias.fit.pct = bias.fit / Y.untrt.true,
         bias.adj.pct = bias.adj / Y.untrt.true)

bias.AME2 <- bias.df %>%
  filter(trans_prob.base1=="0.1265", trans_prob.base2=="0.115") %>%
  # filter(eff.multi != 1) %>%
  group_by(pop.size, trans_prob.base1, trans_prob.base2, model, model.lab, eff.multi) %>%
  summarise(nsim = n(), Y.trt = mean(Y.trt), Y.untrt = mean(Y.untrt), Y.untrt.true = mean(Y.untrt.true), 
            AME.true = mean(AME.true), AME.fit = mean(AME.fit), AME.adj = mean(AME.adj, na.rm = T)) %>%
  mutate(bias.fit = AME.fit - AME.true,
         bias.adj = ifelse(model %in% c("inc", "loginc"), bias.fit, AME.adj - AME.true),
         bias.fit.pct = bias.fit / Y.untrt.true,
         bias.adj.pct = bias.adj / Y.untrt.true)
##############################################################################################################################
bias.original1 <- bias.df %>% 
  filter(trans_prob.base1=="0.115", trans_prob.base2=="0.115") %>%
  group_by(trans_prob.base1, trans_prob.base2, eff.multi, model, model.lab) %>%
  summarise(nsim = n(), eff = mean(effect)) %>%
  merge(eff.truth, by = c("trans_prob.base1", "trans_prob.base2", "eff.multi", "model"))
bias.original1$eff.true[bias.original1$model=="inc"] <- bias.AME1$AME.true[bias.AME1$model=="inc"]
bias.original1$eff.bias <- bias.original1$eff - bias.original1$eff.true
bias.original1$eff.bias.pct <- bias.original1$eff.bias / bias.original1$eff.true

bias.original2 <- bias.df %>% 
  filter(trans_prob.base1=="0.1265", trans_prob.base2=="0.115") %>%
  group_by(trans_prob.base1, trans_prob.base2, eff.multi, model, model.lab) %>%
  summarise(nsim = n(), eff = mean(effect)) %>%
  merge(eff.truth, by = c("trans_prob.base1", "trans_prob.base2", "eff.multi", "model"))
bias.original2$eff.true[bias.original2$model=="inc"] <- bias.AME2$AME.true[bias.AME2$model=="inc"]
bias.original2$eff.bias <- bias.original2$eff - bias.original2$eff.true
bias.original2$eff.bias.pct <- bias.original2$eff.bias / bias.original2$eff.true
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
  geom_point(col=pal[5]) +
  geom_line(col=pal[5]) +
  theme_bw() +
  labs(x = "", y = "Bias in RR (%)") +
  scale_y_continuous(limit=c(-40, 40), breaks = seq(-40, 40, 20)) +
  scale_color_manual(name = "", values = pal) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

p.bias.AME2 <- ggplot(bias.AME2 %>% filter(model=="beta"), aes(eff.multi, bias.adj.pct*100)) +
  geom_point(col=pal[5]) +
  geom_line(col=pal[5]) +
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
format.tbl(power.df2, bias.original2, bias.AME2) %>% kable(format = "latex")
