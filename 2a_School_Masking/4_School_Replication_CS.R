setwd("~/Library/CloudStorage/Dropbox/Research/COVID Rt/School Masking")
rm(list=ls())
library(tidyverse)
library(lubridate)
library(zoo)
library(RStata)
library(EpiEstim)
options("RStata.StataPath"="/Applications/Stata/StataSE.app/Contents/MacOS/stata-se")
options("RStata.StataVersion"=15.1)
options("max.print"=2000)

source("attgt.glm.R")
##########################################################################################################################################
df.clean <- readRDS("School Data/df.clean.rds")
data_in <- df.clean %>% filter(week <= 40)
######################################### INCIDENCE #########################################
inc <- attgt.glm(yname = "PosPer1K", tname = "week", idname = "OrgCode", gname = "treat.time", 
                 control_group="notyettreated", weightsname = "total", data_in=data_in,
                 family = gaussian(link="identity")) # original analysis had 44.9 (32.6, 57.1)
inc
sum(inc)
plot(1:15, inc)
abline(h=0)
######################################### LOG INCIDENCE #########################################
loginc <- attgt.glm(yname = "PosPer1K", tname = "week", idname = "OrgCode", gname = "treat.time", 
                    control_group="notyettreated", weightsname = "total", data_in=data_in,
                    family = poisson(link="log"))
loginc
sum(loginc) # 1.956731
plot(1:length(loginc), loginc)
abline(h=0)

loginc.obs <- data_in %>% 
  group_by(OrgName) %>%
  arrange(week) %>%
  filter(week >= min(data_in$treat.time[data_in$treat.time>0])) %>%
  mutate(loginc.ctl.case = PosPer1K/exp(loginc), # fitted control potential outcome on case scale, cannot convert on the original log scale b/c of zeros
         diff = PosPer1K - loginc.ctl.case # difference to get marginal effect for each unit i to ME
  ) %>%
  filter(treat.time > 0, week >= treat.time) %>%
  group_by(week) %>%
  summarise(diff = mean(diff)) # average over all units to get the average marginal effect
sum(loginc.obs$diff) # 10.0311
######################################### LOG GROWTH #########################################
data_in_growth <- data_in %>%
  group_by(OrgCode, OrgName) %>%
  arrange(week) %>%
  mutate(growth = ifelse(lag(PosPer1K)==0, PosPer1K/1, PosPer1K/lag(PosPer1K, 1)))

loggrowth <- attgt.glm(yname = "growth", tname = "week", idname = "OrgCode", gname = "treat.time", 
                       control_group="notyettreated", weightsname = "total", data_in=data_in_growth,
                       family = poisson(link="log"))
loggrowth
sum(loggrowth)
plot(1:length(loggrowth), loggrowth)
abline(h=0)

growth.obs <- data_in_growth %>% 
  group_by(OrgName) %>%
  arrange(week) %>%
  filter(week >= min(data_in$treat.time[data_in$treat.time>0])-1) %>%
  mutate(growth.ctl = log(growth) - c(0, loggrowth), # append a 0 to keep the last period before treatment
         Pos.ctl = NA) %>%
  filter(treat.time > 0, week >= treat.time-1,
         !is.na(growth))

for (unit in unique(growth.obs$OrgName)) {
  for (time in sort(unique(growth.obs$week[growth.obs$OrgName==unit]))) {
    if ((time+1) == unique(growth.obs$treat.time[growth.obs$OrgName==unit])) {
      growth.obs$Pos.ctl[growth.obs$OrgName==unit & growth.obs$week==time] <- 
        growth.obs$PosPer1K[growth.obs$OrgName==unit & growth.obs$week==time]
    } else {
      time.ind <- ifelse(nrow(growth.obs[growth.obs$OrgName==unit & growth.obs$week==time-1,])==0, 2, 1)
      growth.obs$Pos.ctl[growth.obs$OrgName==unit & growth.obs$week==time] <- 
        growth.obs$Pos.ctl[growth.obs$OrgName==unit & growth.obs$week==time-time.ind] * 
        exp(growth.obs$growth.ctl[growth.obs$OrgName==unit & growth.obs$week==time])
    }
  }
}

growth.obs <- growth.obs %>%
  filter(week >= treat.time) %>%
  mutate(diff = PosPer1K - Pos.ctl) %>%
  group_by(OrgName) %>%
  summarise(diff = sum(diff)) # difference to get marginal effect for each unit
mean(growth.obs$diff) # average over all units to get the AME: -179.4577
######################################## INCIDENCE CI #########################################
CI <- readRDS("CI.rds") %>% rbind(readRDS("CI2.rds"))
CI[,1:6] <- sapply(CI[,1:6], as.numeric)
CI$SE <- (CI$upper-CI$lower) / (qnorm(0.975) * 2)
glist <- sort(unique(data_in$treat.time[data_in$treat.time>0])) # ever-treated groups
weight <- sapply(glist, function(g) mean(data_in$treat.time[data_in$treat.time>0]==g))

CI.inc <- CI %>% filter(model == "Inc")
CI.inc.out <- c()
for (i in 1:5000){
  beta <- sapply(1:nrow(CI.inc), function(gt) {rnorm(1, CI.inc$beta[gt], CI.inc$SE[gt])})
  
  ATT_gt <- matrix(nrow = length(glist), ncol = length(min(glist):max(data_in$week)))
  rownames(ATT_gt) <- glist
  colnames(ATT_gt) <- min(glist):max(data_in$week)
  
  ATT_gt["26",] <- beta[CI.inc$g==26]
  ATT_gt["27",] <- c(NA, beta[CI.inc$g==27])
  ATT_gt["28",] <- c(NA, NA, beta[CI.inc$g==28])
  
  ATT.tmp <- c()
  for (col in 1:ncol(ATT_gt)) {
    ind <- 1:length(na.omit(ATT_gt[,col]))
    ATT.tmp <- c(ATT.tmp, as.numeric(weight[ind] %*% ATT_gt[ind,col]))
  }
  CI.inc.out <- c(CI.inc.out, sum(ATT.tmp))
}

quantile(CI.inc.out, c(0.025, 0.975)) # 34.09092 60.40981
######################################## LOG INCIDENCE CI ########################################
CI.loginc <- CI %>% filter(model == "Log inc")

CI.loginc.ATT <- c() # 1.956731 is the point ATT
CI.loginc.out <- c()
for (i in 1:5000){
  beta <- sapply(1:nrow(CI.loginc), function(gt) {rnorm(1, CI.loginc$beta[gt], CI.loginc$SE[gt])})
  
  ATT_gt <- matrix(nrow = length(glist), ncol = length(min(glist):max(data_in$week)))
  rownames(ATT_gt) <- glist
  colnames(ATT_gt) <- min(glist):max(data_in$week)
  
  ATT_gt["26",] <- beta[CI.loginc$g==26]
  ATT_gt["27",] <- c(NA, beta[CI.loginc$g==27])
  ATT_gt["28",] <- c(NA, NA, beta[CI.loginc$g==28])
  
  ATT.tmp <- c()
  for (col in 1:ncol(ATT_gt)) {
    ind <- 1:length(na.omit(ATT_gt[,col]))
    ATT.tmp <- c(ATT.tmp, as.numeric(weight[ind] %*% ATT_gt[ind,col]))
  }
  
  loginc.obs <- data_in %>% 
    group_by(OrgName) %>%
    arrange(week) %>%
    filter(week >= min(data_in$treat.time[data_in$treat.time>0])) %>%
    mutate(loginc.ctl.case = PosPer1K/exp(ATT.tmp), # fitted control potential outcome on case scale, cannot convert on the original log scale b/c of zeros
           diff = PosPer1K - loginc.ctl.case # difference to get marginal effect for each unit i to ME
    ) %>%
    filter(treat.time > 0, week >= treat.time) %>%
    group_by(week) %>%
    summarise(diff = mean(diff)) # average over all units to get the average marginal effect
  sum(loginc.obs$diff)
  
  CI.loginc.ATT <- c(CI.loginc.ATT, sum(ATT.tmp))
  CI.loginc.out <- c(CI.loginc.out, sum(loginc.obs$diff))
}
saveRDS(CI.loginc.ATT, "Temp/CI.loginc.ATT.rds")
saveRDS(CI.loginc.out, "Temp/CI.loginc.out.rds")

quantile(CI.loginc.ATT, c(0.025, 0.975)) # -12.64823  16.62681
quantile(CI.loginc.out, c(0.025, 0.975)) # -14.57948      43.62281
######################################## LOG GROWTH CI ########################################
CI.growth <- CI %>% filter(model == "Log growth")

CI.growth.ATT <- c() # -2.245529
CI.growth.out <- c() # -179.4577
for (i in 1:2000){
  beta <- sapply(1:nrow(CI.growth), function(gt) {rnorm(1, CI.growth$beta[gt], CI.growth$SE[gt])})
  
  ATT_gt <- matrix(nrow = length(glist), ncol = length(min(glist):max(data_in$week)))
  rownames(ATT_gt) <- glist
  colnames(ATT_gt) <- min(glist):max(data_in$week)
  
  ATT_gt["26",] <- beta[CI.growth$g==26]
  ATT_gt["27",] <- c(NA, beta[CI.growth$g==27])
  ATT_gt["28",] <- c(NA, NA, beta[CI.growth$g==28])
  
  ATT.tmp <- c()
  for (col in 1:ncol(ATT_gt)) {
    ind <- 1:length(na.omit(ATT_gt[,col]))
    ATT.tmp <- c(ATT.tmp, as.numeric(weight[ind] %*% ATT_gt[ind,col]))
  }
  
  CI.growth.ATT <- c(CI.growth.ATT, sum(ATT.tmp))
  
  growth.obs <- data_in_growth %>% 
    group_by(OrgName) %>%
    arrange(week) %>%
    filter(week >= min(data_in$treat.time[data_in$treat.time>0])-1) %>%
    mutate(growth.ctl = log(growth) - c(0, ATT.tmp), # append a 0 to keep the last period before treatment
           Pos.ctl = NA) %>%
    filter(treat.time > 0, week >= treat.time-1,
           !is.na(growth))
  
  for (unit in unique(growth.obs$OrgName)) {
    for (time in sort(unique(growth.obs$week[growth.obs$OrgName==unit]))) {
      if ((time+1) == unique(growth.obs$treat.time[growth.obs$OrgName==unit])) {
        growth.obs$Pos.ctl[growth.obs$OrgName==unit & growth.obs$week==time] <- 
          growth.obs$PosPer1K[growth.obs$OrgName==unit & growth.obs$week==time]
      } else {
        time.ind <- ifelse(nrow(growth.obs[growth.obs$OrgName==unit & growth.obs$week==time-1,])==0, 2, 1)
        growth.obs$Pos.ctl[growth.obs$OrgName==unit & growth.obs$week==time] <- 
          growth.obs$Pos.ctl[growth.obs$OrgName==unit & growth.obs$week==time-time.ind] * 
          exp(growth.obs$growth.ctl[growth.obs$OrgName==unit & growth.obs$week==time])
      }
    }
  }
  
  growth.obs <- growth.obs %>%
    filter(week >= treat.time) %>%
    mutate(diff = PosPer1K - Pos.ctl) %>%
    group_by(OrgName) %>%
    summarise(diff = sum(diff)) # difference to get marginal effect for each unit
  CI.growth.out <- c(CI.growth.out, mean(growth.obs$diff))
}
saveRDS(CI.growth.ATT, "Temp/CI.growth.ATT.rds")
saveRDS(CI.growth.out, "Temp/CI.growth.out.rds")

quantile(CI.growth.ATT, c(0.025, 0.975)) # -6.063438  1.649841
quantile(CI.growth.out, c(0.025, 0.975)) # -8040.6861   112.8495
