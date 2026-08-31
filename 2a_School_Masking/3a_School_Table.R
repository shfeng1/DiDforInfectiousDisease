here::i_am("2a_School_Masking/3a_School_Table.R")
source("./global_options.R")

# Run each analysis in its own temporary environment so that the incidence
# script's workspace cleanup does not remove results needed by this table.
school_incidence_results <- new.env(parent=globalenv())
source(
  "./2a_School_Masking/1_School_Incidence.R",
  local=school_incidence_results
)

school_growth_results <- new.env(parent=globalenv())
source(
  "./2a_School_Masking/2_School_Growth.R",
  local=school_growth_results
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

school_markers <- c(
  table3_significance_marker(school_incidence_results$inc_effect_5, 0),
  table3_significance_marker(school_incidence_results$loginc_effect_5, 1),
  table3_significance_marker(school_growth_results$growth_effect_5, 1),
  table3_significance_marker(school_incidence_results$inc_effect_15, 0),
  table3_significance_marker(school_incidence_results$loginc_effect_15, 1),
  table3_significance_marker(school_growth_results$growth_effect_15, 1)
)

school_effects <- list(
  school_incidence_results$inc_effect_5,
  school_incidence_results$loginc_effect_5,
  school_growth_results$growth_effect_5,
  school_incidence_results$inc_effect_15,
  school_incidence_results$loginc_effect_15,
  school_growth_results$growth_effect_15
)
school_AMEs <- list(
  school_incidence_results$inc_AME_5,
  school_incidence_results$loginc_AME_5,
  school_growth_results$growth_AME_5,
  school_incidence_results$inc_AME_15,
  school_incidence_results$loginc_AME_15,
  school_growth_results$growth_AME_15
)
school_effect_digits <- c(1, 2, 2, 1, 2, 2)

Table3_School <- tibble(
  Study=rep(
    "Removing school mask mandates\nin Massachusetts\n(per 1,000 students and staff)",
    6
  ),
  `Length of post-intervention period`=rep(c("5 weeks", "15 weeks"), each=3),
  `Outcome specification`=c(
    "Incidence", "Log incidence", "Log growth",
    "Incidence$^{3}$", "Log incidence", "Log growth"
  ),
  `Treatment effect (95% CI)`=map2_chr(
    seq_along(school_effects),
    school_effect_digits,
    ~format_table3_ci(
      school_effects[[.x]], .y, school_markers[.x]
    )
  ),
  `Average marginal effect (95% CI)`=map2_chr(
    seq_along(school_AMEs),
    school_markers,
    ~format_table3_ci(school_AMEs[[.x]], 1, .y)
  )
)

tbl3.school <- Table3_School %>%
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
  collapse_rows(columns=c(1, 2), valign="middle", latex_hline="major") %>%
  row_spec(3, extra_latex_after="\\cmidrule{2-5}")

tbl3.school
