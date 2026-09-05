rm(list=ls())
here::i_am("2a_School_Masking/1_School_Incidence.R")
source("./global_options.R")
source("./1a_Scripts/0_boottest.glm.R", local=TRUE)
source("./2a_School_Masking/0_School_AME_Helpers.R", local=TRUE)
############################################# MAKE DATA ##############################################
data_in <- readRDS("./0_Data/School_Cleaned.rds") %>% filter(week <= 40) %>%
  group_by(OrgCode, OrgName) %>% arrange(week) %>%
  mutate(growth = ifelse(lag(PosPer1K)==0, PosPer1K/1, PosPer1K/lag(PosPer1K, 1)))
tname = "week";idname = "OrgCode";gname = "treat.time";control_group="notyettreated";weightsname = "total"
######################################### For Inc and Log Inc #########################################
yname = "PosPer1K"
y <- data_in[[yname]]
time <- data_in[[tname]]
ID <- as.numeric(data_in[[idname]])
group <- data_in[[gname]]
glist <- sort(unique(group[group>0])) # ever-treated groups
if (is.null(weightsname)) { # check if weights are specified
  wt <- rep(1, length(y))
} else {
  wt <- data_in[[weightsname]]
}

df_sunab <- data.frame(cbind(ID, group, time, y, wt)) %>%
  mutate(group = ifelse(group==0, 10000, group),
         time_to_trt = ifelse(group==10000, -1, time - group))

gweight <- sapply(glist, function(g) mean(df_sunab$group[df_sunab$group %in% glist]==g))
################################################# INC ################################################
fit <- feols(y ~ sunab(group, time, ref.c = 10000) | ID + time, weights = df_sunab$wt, data = df_sunab)
time_to_trt <- as.numeric(sapply(names(coef(fit)), function(var) substr(var, 7, nchar(var))))
ATT_gt <- boottest.glm(fit, gweight=gweight, model="inc")
rownames(ATT_gt) <- time_to_trt
ATT_boot <- colSums(ATT_gt[rownames(ATT_gt) >= 0,])
ATT <- sum(coef(fit)[time_to_trt >= 0])
inc_effect_15 <- c(estimate=ATT,
                   lower=unname(quantile(ATT_boot, 0.025)),
                   upper=unname(quantile(ATT_boot, 0.975)))
inc_AME_15 <- inc_effect_15
##################################################   KEEP 5 WEEKS POST INTERVENTION
ATT_boot <- colSums(ATT_gt[rownames(ATT_gt) %in% (0:4),])
ATT <- sum(coef(fit)[time_to_trt %in% (0:4)])
inc_effect_5 <- c(estimate=ATT,
                  lower=unname(quantile(ATT_boot, 0.025)),
                  upper=unname(quantile(ATT_boot, 0.975)))
inc_AME_5 <- inc_effect_5
############################################### LOG INC ##############################################
fit <- fepois(y ~ sunab(group, time, ref.c = 10000) | ID + time, weights = df_sunab$wt, data = df_sunab)
time_to_trt <- as.numeric(sapply(names(coef(fit)), function(var) substr(var, 7, nchar(var))))
ATT_gt <- boottest.glm(fit, gweight=gweight, model="loginc")
rownames(ATT_gt) <- time_to_trt
ATT_boot <- colMeans(ATT_gt[rownames(ATT_gt) >= 0,])
ATT <- mean(coef(fit)[time_to_trt >= 0])
AMEs <- apply(ATT_gt, 2, function(coef) school_loginc_ame(coef, time_to_trt, df_sunab))
loginc_effect_15 <- c(estimate=exp(ATT),
                      lower=unname(exp(quantile(ATT_boot, 0.025))),
                      upper=unname(exp(quantile(ATT_boot, 0.975))))
loginc_AME_15 <- c(estimate=school_loginc_ame(coef(fit), time_to_trt, df_sunab),
                   lower=unname(quantile(AMEs, 0.025)),
                   upper=unname(quantile(AMEs, 0.975)))
##################################################   KEEP 5 WEEKS POST INTERVENTION
ATT_boot <- colMeans(ATT_gt[rownames(ATT_gt) %in% (0:4),])
ATT <- mean(coef(fit)[time_to_trt %in% (0:4)])
AMEs <- apply(ATT_gt, 2, function(coef) school_loginc_ame(coef, time_to_trt, df_sunab, 0:4))
loginc_effect_5 <- c(estimate=exp(ATT),
                     lower=unname(exp(quantile(ATT_boot, 0.025))),
                     upper=unname(exp(quantile(ATT_boot, 0.975))))
loginc_AME_5 <- c(estimate=school_loginc_ame(coef(fit), time_to_trt, df_sunab, 0:4),
                  lower=unname(quantile(AMEs, 0.025)),
                  upper=unname(quantile(AMEs, 0.975)))
