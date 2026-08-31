# Parallel Trends in an Unparalleled Pandemic

Replication code for *Parallel Trends in an Unparalleled Pandemic: Difference-in-Differences for Infectious Disease Policy Evaluation*. The repository contains the simulation studies and the Massachusetts school-masking and Kansas county-masking reanalyses.

## Requirements

The analyses were developed with R 4.2.2 and Stata/SE 15.1. Stata must have `boottest` installed. Parallel simulation uses `doMC`, so the code is intended for macOS or Linux.

Install the R dependencies with:

```r
install.packages(c(
  "tidyverse", "data.table", "readxl", "doMC", "RColorBrewer",
  "sandwich", "fixest", "lmtest", "did", "haven", "RStata",
  "ggpubr", "kableExtra", "here", "vroom", "foreign"
))
```

In [`global_options.R`](global_options.R), set `RStata.StataPath` and `RStata.StataVersion` for your installation. Change `n.cores` there if five workers are not appropriate for your machine.

## Reproduce everything

`0_Master_script.R` is the single entry point for the complete workflow. It runs the simulation, summary, empirical, and appendix scripts needed to recreate every manuscript table and figure; writes generated figures and simulation files to `4_Output/`; and prints the tables and empirical estimates.

1. Clone or download the repository and open `Parallel_Trends_Replication.Rproj` (or start R with the repository root as the working directory).
2. Configure Stata and the worker count in `global_options.R`.
3. For a clean simulation run, start from the existing `.rds` files in `4_Output/` folder. Running the commented out simulation scripts called in `0_Master_script.R` would allow the user to reproduce the exact simulation results; however, please be mindful of extensive computing time (>8 hours).
4. Start a clean R session and run:

```r
source("0_Master_script.R")
```

A full regeneration is long-running (more than eight hours; hardware and Stata I/O materially affect runtime).

## Data and repository layout

- `0_Data/`: analysis inputs. The processed `School_Cleaned.rds` and `Kansas_Cleaned.rds` files required by the master workflow are included.
- `1a_Scripts/`: simulation models, estimators, and simulation drivers.
- `1b_Summarize/`: scripts that create the simulation figures and manuscript tables.
- `2a_School_Masking/`: Massachusetts school-masking reanalysis.
- `2b_Kansas_Masking/`: Kansas county-masking reanalysis.
- `4_Output/`: generated figures and simulation results.

The optional raw-data cleaning scripts are `2a_School_Masking/0_School_Clean_Data.R` and `2b_Kansas_Masking/0_Kansas_Clean_Data.R`. The latter additionally expects `0_Data/covidestim-daily-fips.csv.xz`, which is not required when using the included processed Kansas data.
