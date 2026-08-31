here::i_am("2b_Kansas_Masking/6a_Kansas_Table.R")
source("./global_options.R")

run_table3_analysis <- function(path) {
  results <- new.env(parent=globalenv())
  source(path, local=results)
  results
}

# Run each model in a separate temporary environment because the standalone
# analysis scripts clear their own workspaces.
kansas_incidence_results <- run_table3_analysis(
  "./2b_Kansas_Masking/1_Kansas_Incidence.R"
)
kansas_growth_results <- run_table3_analysis(
  "./2b_Kansas_Masking/2_Kansas_Growth.R"
)
kansas_Rt_results <- run_table3_analysis(
  "./2b_Kansas_Masking/3_Kansas_Rt_Exposure.R"
)
kansas_beta_exposure_results <- run_table3_analysis(
  "./2b_Kansas_Masking/4_Kansas_Beta_Exposure.R"
)
kansas_beta_COVIDEstim_results <- run_table3_analysis(
  "./2b_Kansas_Masking/5_Kansas_Beta_COVIDEstim.R"
)

format_table3_ci <- function(x, digits, marker="") {
  required <- c("estimate", "lower", "upper")
  if (!all(required %in% names(x))) {
    stop("Table 3 results must contain estimate, lower, and upper values.")
  }
  values <- formatC(x[required], format="f", digits=digits)
  paste0(values[1], marker, " (", values[2], ", ", values[3], ")")
}

table3_significance_marker <- function(x, null) {
  if (x["lower"] > null || x["upper"] < null) "$^{*}$" else ""
}

kansas_effects <- list(
  kansas_incidence_results$inc_effect_20,
  kansas_incidence_results$loginc_effect_20,
  kansas_growth_results$growth_effect_20,
  kansas_Rt_results$Rt_exposure_effect_20,
  kansas_beta_exposure_results$beta_exposure_effect_20,
  kansas_beta_COVIDEstim_results$beta_COVIDEstim_effect_20
)
kansas_AMEs <- list(
  kansas_incidence_results$inc_AME_20,
  kansas_incidence_results$loginc_AME_20,
  kansas_growth_results$growth_AME_20,
  kansas_Rt_results$Rt_exposure_AME_20,
  kansas_beta_exposure_results$beta_exposure_AME_20,
  kansas_beta_COVIDEstim_results$beta_COVIDEstim_AME_20
)
kansas_effect_digits <- c(1, 2, 2, 2, 2, 2)
kansas_markers <- c(
  table3_significance_marker(kansas_effects[[1]], 0),
  table3_significance_marker(kansas_effects[[2]], 1),
  "",
  "$^{\\dagger}$",
  table3_significance_marker(kansas_effects[[5]], 1),
  "$^{\\dagger}$"
)

Table3_Kansas <- tibble(
  Study=rep(
    "Mask mandates in Kansas counties\n(per 100,000 residents)",
    6
  ),
  `Length of post-intervention period`=rep("20 weeks", 6),
  `Outcome specification`=c(
    "Incidence$^{3}$", "Log incidence", "Log growth",
    "Log $R_t$ (instantaneous)",
    "Log $\\beta_t$ (instantaneous)",
    "Log $\\beta_t$ (COVIDestim)"
  ),
  `Treatment effect (95% CI)`=map2_chr(
    seq_along(kansas_effects),
    kansas_effect_digits,
    ~format_table3_ci(
      kansas_effects[[.x]], .y, kansas_markers[.x]
    )
  ),
  `Average marginal effect (95% CI)`=map2_chr(
    seq_along(kansas_AMEs),
    kansas_markers,
    ~format_table3_ci(kansas_AMEs[[.x]], 1, .y)
  )
)

tbl3.kansas <- Table3_Kansas %>%
  mutate(Study=linebreak(Study, align="c")) %>%
  kable(
    format="latex", booktabs=TRUE, escape=FALSE, linesep="",
    align=c("c", "c", "c", "c", "c"),
    col.names=c(
      "Study",
      linebreak("Length of post-\nintervention period", align="c"),
      linebreak("Outcome\nspecification", align="c"),
      linebreak("Treatment effect$^{1}$\n(95\\% CI)", align="c"),
      linebreak("Average marginal\neffect$^{2}$ (95\\% CI)", align="c")
    )
  ) %>%
  collapse_rows(columns=c(1, 2), valign="middle", latex_hline="major")

tbl3.kansas
