rm(list=ls())
here::i_am("2a_School_Masking/1_School_Incidence.R")
source("./global_options.R")
source("./1a_Scripts/0_boottest.glm.R")
source("./1a_Scripts/0_Run_Estimators.R")
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
ATT_gt <- boottest.glm(fit, M=1000, gweight=gweight, model="inc")
rownames(ATT_gt) <- time_to_trt
ATT <- colSums(ATT_gt[rownames(ATT_gt) >= 0,])

print("------------------ INCIDENCE MODEL ------------------")
print(paste0("15-week treatment effect and AME are the same: ", round(mean(ATT), 1), " with CI: (",  
             round(quantile(ATT, 0.025), 1), ", ", round(quantile(ATT, 0.975), 1), ")"))
# 15-week treatment effect and AME are the same: 47.8 with CI: (38.4, 57.4)
##################################################   KEEP 5 WEEKS POST INTERVENTION
rownames(ATT_gt) <- time_to_trt
ATT <- colSums(ATT_gt[rownames(ATT_gt) %in% (0:4),])

print(paste0("5-week treatment effect and AME are the same: ", round(mean(ATT), 1), " with CI: (",  
             round(quantile(ATT, 0.025), 1), ", ", round(quantile(ATT, 0.975), 1), ")"))
# 5-week treatment effect and AME are the same: 8.5 with CI: (5.6, 11.5)
############################################### LOG INC ##############################################
fit <- fepois(y ~ sunab(group, time, ref.c = 10000) | ID + time, weights = df_sunab$wt, data = df_sunab)
time_to_trt <- as.numeric(sapply(names(coef(fit)), function(var) substr(var, 7, nchar(var))))
ATT_gt <- boottest.glm(fit, M=1000, gweight=gweight, model="loginc")
rownames(ATT_gt) <- time_to_trt
ATT <- colMeans(ATT_gt[rownames(ATT_gt) >= 0,])

print("------------------ LOG INCIDENCE MODEL ------------------")
print(paste0("15-week treatment effect: ", round(exp(mean(ATT)), 2), " with CI: (",  
             round(exp(quantile(ATT, 0.025)), 2), ", ", round(exp(quantile(ATT, 0.975)), 2), ")"))
# 15-week treatment effect: 1.19 with CI: (0.89, 1.56)
2*pnorm(abs(mean(ATT)/sd(ATT)), lower.tail=FALSE) # p = 0.239113

AMEs <- apply(ATT_gt, 2, function(coef) loginc_AME(coef))
print(paste0("15-week AME: ", round(loginc_AME(rowMeans(ATT_gt)), 1), " with CI: (",  
             round(quantile(AMEs, 0.025), 1), ", ", round(quantile(AMEs, 0.975), 1), ")"))
# 15-week AME: 12.9 with CI: (-33.2, 47.0)
##################################################   KEEP 5 WEEKS POST INTERVENTION
rownames(ATT_gt) <- time_to_trt
ATT <- colMeans(ATT_gt[rownames(ATT_gt) %in% (0:4),])
print(paste0("5-week treatment effect: ", round(exp(mean(ATT)), 2), " with CI: (",  
             round(exp(quantile(ATT, 0.025)), 2), ", ", round(exp(quantile(ATT, 0.975)), 2), ")"))
# 5-week treatment effect: 1.62 with CI: (1.24, 2.09)
2*pnorm(abs(mean(ATT)/sd(ATT)), lower.tail=FALSE) # p = 0.0002542965

AMEs <- apply(ATT_gt, 2, function(coef) loginc_AME(coef, c(0:4)))
print(paste0("5-week AME: ", round(loginc_AME(coef=rowMeans(ATT_gt), subset=c(0:4)), 1), " with CI: (",  
             round(quantile(AMEs, 0.025), 1), ", ", round(quantile(AMEs, 0.975), 1), ")"))
# 5-week AME: 6.9 with CI: (1.6, 10.8)
