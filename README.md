# Parallel Trends in an Unparalleled Pandemic

Replication code for *Parallel Trends in an Unparalleled Pandemic: Difference-in-Differences for Infectious Disease Policy Evaluation*. The repository contains the simulation studies and the Massachusetts school-masking and Kansas county-masking reanalyses.

## Requirements

The analyses were developed with R 4.2.2 and Stata/SE 15.1. Stata must have `boottest` installed. Parallel simulation uses `doMC`, so the code is intended for macOS or Linux.

Install the R dependencies with:

```r
install.packages(c(
  "tidyverse", "data.table", "readxl", "doMC", "RColorBrewer", "sandwich", "fixest", "lmtest", 
  "did", "haven", "RStata", "ggpubr", "kableExtra", "here", "vroom", "foreign"
))
```

In [`global_options.R`](global_options.R), set `RStata.StataPath` and `RStata.StataVersion` for your installation. Change `n.cores` there if five workers are not appropriate for your machine.

## Reproduce everything

`0_Master_script.R` is the single entry point for the complete workflow. By default, it reproduces the tables and figures using the supplied simulation results in `4_Output/`, runs the empirical analyses and the Figure 1 simulations, saves generated figures to `4_Output/`, and prints the tables and empirical estimates. The code block for regenerating the stored simulation results is commented out by default.

1. Clone or download the repository and open `Parallel_Trends_Replication.Rproj` (or start R with the repository root as the working directory).
2. Configure Stata and the worker count in `global_options.R`.
3. To reproduce the tables and figures using the supplied simulation results, keep the existing `.rds` files in `4_Output/`. To regenerate all simulations from scratch, first archive the current `.rds` output files under `4_Output` folder by moving them to a different folder, then uncomment the code block under `0_Master_script.R` to enable the simulation regeneration scripts. However, please be mindful of the extensive computing time (>8 hours).
4. Start a clean R session and run:

```r
source("0_Master_script.R")
```

## Manuscript outputs and scripts

The table below maps manuscript outputs to the scripts called by `0_Master_script.R`. Simulation summary scripts use the stored results in `4_Output/`; the master script saves the figures.

| Output | Script |
| --- | --- |
| Figure 1 | [`1b_Summarize/1_Simulate_comparison_of_models.R`](1b_Summarize/1_Simulate_comparison_of_models.R) |
| Figure 2, Table 2 | [`1b_Summarize/2a_SIR_summ.R`](1b_Summarize/2a_SIR_summ.R) (Figure 2 and Table 2)<br>[`1b_Summarize/2b_SEIR_summ.R`](1b_Summarize/2b_SEIR_summ.R) (Table 2)<br>[`1b_Summarize/3a_Misspecify_GI_summ.R`](1b_Summarize/3a_Misspecify_GI_summ.R) (Table 2)<br>[`1b_Summarize/3b_Misspecify_SEIR_to_SIR_summ.R`](1b_Summarize/3b_Misspecify_SEIR_to_SIR_summ.R) (Table 2) |
| Table 3 (Massachusetts) | [`2a_School_Masking/3a_School_Table.R`](2a_School_Masking/3a_School_Table.R) |
| Table 3 (Kansas) | [`2b_Kansas_Masking/6a_Kansas_Table.R`](2b_Kansas_Masking/6a_Kansas_Table.R) |
| Appendix tables and figures | [`2a_School_Masking/3b_School_Graph.R`](2a_School_Masking/3b_School_Graph.R) (Figures A1 and A2)<br>[`2b_Kansas_Masking/6b_Kansas_Graph.R`](2b_Kansas_Masking/6b_Kansas_Graph.R) (Figure A3)<br>[`2a_School_Masking/4_School_Callaway_SantAnna.R`](2a_School_Masking/4_School_Callaway_SantAnna.R) (Table A2)<br>[`1b_Summarize/4_SIR_long_time_summ.R`](1b_Summarize/4_SIR_long_time_summ.R) (Table A3)<br>[`1b_Summarize/5_Small_N1_summ.R`](1b_Summarize/5_Small_N1_summ.R) (Table A4) |

## Data and repository layout

- `0_Data/`: analysis inputs. The processed `School_Cleaned.rds` and `Kansas_Cleaned.rds` files required by the master workflow are included.
- `1a_Scripts/`: simulation models, estimators, and simulation drivers.
- `1b_Summarize/`: scripts that create the simulation figures and manuscript tables.
- `2a_School_Masking/`: Massachusetts school-masking reanalysis.
- `2b_Kansas_Masking/`: Kansas county-masking reanalysis.
- `4_Output/`: generated figures and simulation results.

The optional raw-data cleaning scripts are `2a_School_Masking/0_School_Clean_Data.R` and `2b_Kansas_Masking/0_Kansas_Clean_Data.R`. The latter additionally expects `0_Data/covidestim-daily-fips.csv.xz`, which is not required when using the included processed Kansas data.
