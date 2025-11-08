library(tidyverse)
library(data.table)
library(readxl)
library(doParallel)
library(EpiEstim)
library(RColorBrewer)
library(fwildclusterboot)
library(sandwich)
library(clubSandwich)
library(fixest)
library(lmtest)
library(haven)
library(RStata)
library(ggpubr)
library(kableExtra)

options(warn=-1)
options(dplyr.summarise.inform = FALSE)
options("RStata.StataPath"="/Applications/Stata/StataSE.app/Contents/MacOS/stata-se")
options("RStata.StataVersion"=15.1)

doMC::registerDoMC(cores = parallel::detectCores())

pal <- c("darkgrey", "#984EA3", "#FF5932", "olivedrab4", "tan3", "#377EB8")

inf_mean <- 5*2; delta <- 3*2; pop.size <- 1e4; seed1 <- seed2 <- 100
N <- 50; N1_ratio <- 0.5; N1 <- N * N1_ratio
T0 <- 4*7; T1 <- 3*7; burnin <- 5*7; agg <- 7 # aggregate to weekly