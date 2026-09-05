rm(list=ls())
here::i_am("1b_Summarize/4_SIR_long_time_summ.R")
source("./global_options.R")

p_out <- readRDS("./4_Output/SIR_long_time_T1.rds")
T1.values.weeks <- sort(unique(p_out$T1))/agg

bias.AME <- p_out %>%
  mutate(AME.true=Y.trt-Y.untrt.true) %>%
  group_by(T1, model, eff.multi) %>%
  summarise(AME.true=mean(AME.true), AME.fit=mean(AME),
            Y.untrt.true=mean(Y.untrt.true), .groups="drop") %>%
  mutate(T1.weeks=T1/agg,
         bias.fit=AME.fit-AME.true,
         bias.pct=100*bias.fit/Y.untrt.true)

tbl.long.time.data <- bias.AME %>%
  group_by(model, T1.weeks) %>%
  summarise(bias.raw=round(mean(abs(bias.pct)), 1), .groups="drop") %>%
  mutate(model=factor(model, levels=model.list, labels=unname(model.labels[model.list])),
         T1.weeks=factor(T1.weeks, levels=T1.values.weeks),
         bias.raw=paste0(format(bias.raw, nsmall=1), "\\%")) %>%
  dplyr::select(model, T1.weeks, bias.raw) %>%
  pivot_wider(names_from=T1.weeks, values_from=bias.raw) %>%
  arrange(model)

tbl.long.time <- tbl.long.time.data %>%
  kable(format="latex", booktabs=TRUE, escape=FALSE,
        col.names=c("Outcome specification", paste0("T1 = ", T1.values.weeks)))
