rm(list=ls())
here::i_am("2b_Kansas_Masking/6_Kansas_Graph.R")
source("./global_options.R")

df.in <- readRDS("./0_Data/Kansas_Cleaned.rds") %>%
  mutate(col = ifelse(trt.unit, "treated", "control"))

FigureA3 <- ggplot(df.in, aes(x=as.Date(start_date, format = "%Y-%m-%d"), y=stnnewcases7davg, group=ncounty, col=col)) + 
  geom_line(linewidth=0.2, alpha=0.7) +
  geom_line(data = subset(df.in, col=="treated"), color="tan3", linewidth=0.2, alpha=0.7) +
  ylab("Cases per 100K") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  scale_y_continuous(limits = c(0, 650), breaks = seq(0, 650, 150), expand = c(0.02, 0.02)) +
  scale_colour_manual(name = "Kansas counties", values = c("gray", "tan3")) + 
  theme_bw() +
  xlab("") + 
  theme_bw() + geom_vline(col = "darkgray", lty = 2, xintercept = as.Date("2020-07-03")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), legend.position = "right",
        legend.key.size = unit(2, "line"), legend.margin = margin(0, 0, 0, 0),
        legend.spacing.x = unit(0, "mm"), legend.spacing.y = unit(0, "mm")) +
  guides(colour = guide_legend(override.aes = list(linewidth=0.7)))
FigureA3
