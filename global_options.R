library(tidyverse)
library(data.table)
library(readxl)
library(doMC)
library(sandwich)
library(fixest)
library(lmtest)
library(did)
library(RStata)
library(ggpubr)
library(kableExtra)

options(warn=-1)
options(dplyr.summarise.inform = FALSE)
options("RStata.StataPath"="/Applications/Stata/StataSE.app/Contents/MacOS/stata-se")
options("RStata.StataVersion"=15.1)

n.cores <- 5L
doMC::registerDoMC(cores = n.cores)
simulation.chunk.size <- n.cores

make_rng_streams <- function(seed, n) {
  n <- as.integer(n)
  set.seed(as.integer(seed), kind="L'Ecuyer-CMRG")
  streams <- vector("list", n)
  if (n>0L) {
    streams[[1L]] <- .Random.seed
    if (n>1L) {
      for (i in 2:n) streams[[i]] <- parallel::nextRNGStream(streams[[i-1L]])
    }
  }
  streams
}

set_simulation_stream <- function(stream) {
  assign(".Random.seed", stream, envir=.GlobalEnv)
}

is_allowed_write_error <- function(e) {
  msg <- conditionMessage(e)
  grepl("unable to open file for writing", msg, fixed=TRUE) &&
    (grepl("Stale NFS file handle", msg, fixed=TRUE) ||
       grepl("Operation timed out", msg, fixed=TRUE))
}

capture_simulation_error <- function(expr) {
  tryCatch(
    force(expr),
    epidemic_extinct=function(e) NULL,
    error=function(e) {
      if (is_allowed_write_error(e)) return(NULL)
      stop(e)
    }
  )
}

pal <- c("darkgrey", "#984EA3", "#FF5932", "olivedrab4", "tan3", "#377EB8")
model.list <- c("inc", "loginc", "growth", "Rt_exposure", "beta_exposure")

nsim <- 1000 # 1K run at a time.
inf_mean <- 5*2; delta <- 3*2; pop.size <- 1e4; seed1 <- seed2 <- 100
N <- 50; N1_ratio <- 0.5; N1 <- N * N1_ratio
T0 <- 4*7; T1 <- 3*7; burnin <- 5*7; end_buffer <- 3*7; agg <- 7
