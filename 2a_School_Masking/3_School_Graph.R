rm(list=ls())
here::i_am("2a_School_Masking/3_School_Graph.R")
source("./global_options.R")
source("./1a_Scripts/0_boottest.glm.R")
source("./1a_Scripts/0_Run_Estimators.R")
############################################# MAKE DATA ##############################################
data_in <- readRDS("./0_Data/School_Cleaned.rds") %>%
  mutate(col = ifelse(OrgName == "Boston", "Boston", ifelse(OrgName=="Chelsea", "Chelsea", "Comparison districts")),
         col = factor(col, levels = c("Comparison districts", "Boston", "Chelsea")))
######################################################################################################
# Figure A1
FigureA1 <- ggplot(data_in, aes(x=week, y=PosPer1K, group=OrgName, col=col)) + 
  geom_line(linewidth=0.15) +
  geom_line(data = subset(data_in,  OrgName=="Boston"), color = "tan3") +
  geom_line(data = subset(data_in,  OrgName=="Chelsea"), color = "steelblue3") +
  ylab("Student Cases Per 1K") +
  scale_x_continuous(breaks = seq(1, 41, 2)) +
  scale_y_continuous(limits = c(0, 180), breaks = seq(0, 180, 30)) +
  scale_colour_manual(name = "School districts", values = c("gray", "tan3", "steelblue3")) + 
  theme_bw() +
  xlab("") + 
  theme(legend.position = "right",
        legend.key.size = unit(2, "line"),
        panel.grid = element_blank()) +
  guides(colour = guide_legend(override.aes = list(linewidth=0.7)))
FigureA1
######################################################################################################
# Figure A2
df.inc <- readRDS("./0_Data/School_Cleaned.rds") %>% filter(week <= 40) %>%
  mutate(OrgCode = as.numeric(factor(OrgCode)))
######################################### INCIDENCE #########################################
inc.fit <- att_gt(yname = "PosPer1K", tname = "week", idname = "OrgCode", gname = "treat.time", 
                 control_group="notyettreated", weightsname = "total", data=df.inc) # original analysis had 44.9 (32.6, 57.1)

df.inc <- readRDS("./0_Data/df.inc.smooth.gp.rds")
inc.fit <- att_gt(yname = "PosPer1K", tname = "biweek", idname = "OrgCode", 
                  gname = "gname", control_group = "notyettreated", weightsname = "wt.total", 
                  data = df.inc, bstrap = TRUE, cband = TRUE)
inc.fit.agg <- aggte(inc.fit, type = "dynamic", na.rm = T)
FigureA2 <- ggdid(inc.fit.agg) + 
  geom_hline(yintercept=0, color="darkgray", linetype="dashed") +
  geom_vline(xintercept=0, color="darkgray", linetype="dashed") +
  scale_x_continuous(breaks = seq(-25, 15, 5)) +
  scale_y_continuous(limits = c(-15, 20.5), breaks = seq(-15, 20, 5)) +
  theme_bw() + xlab("Time to treatment") + ylab("Estimated treatment effect") + 
  ggtitle("Incidence model: replication of original analysis") +
  theme(plot.title=element_text(size=16, face="bold", hjust=0.5), 
        legend.position="none", panel.grid = element_blank())
FigureA2
