rm(list=ls())
here::i_am("1b_Summarize/4_SEIR_long_time_summ.R")
source("./global_options.R")

models <- c("growth", "Rt_exposure", "beta_exposure")
T1.values.weeks <- c(3, 5, 10, 15, 20)

p_out <- readRDS("./4_Output/SEIR_long_time_T1.rds") %>%
  filter(
    model %in% models,
    eff.multi %in% c(0.9, 0.95, 1),
    (T1/agg) %in% T1.values.weeks
  )

bias.AME <- p_out %>%
  mutate(AME.true=Y.trt-Y.untrt.true) %>%
  group_by(T1, model, eff.multi) %>%
  summarise(
    AME.true=mean(AME.true),
    AME.fit=mean(AME),
    .groups="drop"
  ) %>%
  mutate(
    T1.weeks=T1/agg,
    bias.fit=AME.fit-AME.true
  )

model.labels <- c(
  growth="Log growth",
  Rt_exposure="Log $R_t$",
  beta_exposure="Log $\\beta_t$"
)

tbl.long.time.data <- bias.AME %>%
  group_by(model, T1.weeks) %>%
  summarise(
    mean=format(round(mean(abs(bias.fit)), 1), nsmall=1),
    min=format(round(min(bias.fit), 1), nsmall=1),
    max=format(round(max(bias.fit), 1), nsmall=1),
    .groups="drop"
  ) %>%
  mutate(
    bias.raw=paste0(mean, " (", min, ", ", max, ")"),
    model=factor(model, levels=models, labels=unname(model.labels[models])),
    T1.weeks=factor(T1.weeks, levels=T1.values.weeks)
  ) %>%
  dplyr::select(model, T1.weeks, bias.raw) %>%
  pivot_wider(names_from=T1.weeks, values_from=bias.raw) %>%
  arrange(model)

tbl.long.time <- tbl.long.time.data %>%
  kable(
    format="latex", booktabs=TRUE, escape=FALSE,
    col.names=c("Outcome specification", paste0("T1 = ", T1.values.weeks))
  )
