here::i_am("1b_Summarize/5_Small_N1_summ.R")
source("./global_options.R")
source("./1a_Scripts/0_Format_Table.R")

p_out <- readRDS("./4_Output/SIR_N1=5.rds") %>% filter(model %in% model.list)
##############################################################################################################################
# Power / type I error rate
power.df2 <- p_out %>%
  group_by(model, eff.multi) %>%
  summarise(wild = mean(wild < 0.05, na.rm = TRUE)*100,
            normal = mean(normal < 0.05, na.rm = TRUE)*100,
            nsim = n())
##############################################################################################################################
# Make kable
# Type I error rate
typeIerr <- power.df2 %>%
  filter(eff.multi==1) %>%
  group_by(model) %>%
  summarise(error_wild = round(mean(wild)), error_norm = round(mean(normal)), .groups = "drop") %>%
  dplyr::select(model, error_wild, error_norm)

# Powers at selected effect sizes
power90 <- power.df2 %>%
  filter(eff.multi==0.90) %>%
  group_by(model) %>%
  summarise(wild90 = round(mean(wild)), norm90 = round(mean(normal)), .groups = "drop") %>%
  dplyr::select(model, wild90, norm90)

power95 <- power.df2 %>%
  filter(eff.multi==0.95) %>%
  group_by(model) %>%
  summarise(wild95 = round(mean(wild)), norm95 = round(mean(normal)), .groups = "drop") %>%
  dplyr::select(model, wild95, norm95)

power105 <- power.df2 %>%
  filter(eff.multi==1.05) %>%
  group_by(model) %>%
  summarise(wild105 = round(mean(wild)), norm105 = round(mean(normal)), .groups = "drop") %>%
  dplyr::select(model, wild105, norm105)

power110 <- power.df2 %>%
  filter(eff.multi==1.1) %>%
  group_by(model) %>%
  summarise(wild110 = round(mean(wild)), norm110 = round(mean(normal)), .groups = "drop") %>%
  dplyr::select(model, wild110, norm110)

power.out <- typeIerr %>% merge(power90) %>% merge(power95) %>% merge(power105) %>% merge(power110)

# Put in Latex format for wild score bootstrap
tblA3_wild <- power.out %>%
  dplyr::select(model, error_wild, wild90, wild95, wild105, wild110) %>%
  arrange(match(model, model.list)) %>%
  kable(format = "latex")

# Put in Latex format for normal-based SE
tblA3_norm <- power.out %>%
  dplyr::select(model, error_norm, norm90, norm95, norm105, norm110) %>%
  arrange(match(model, model.list)) %>%
  kable(format = "latex")
