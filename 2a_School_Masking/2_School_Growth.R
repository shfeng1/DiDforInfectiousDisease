# rm(list=ls())
here::i_am("2a_School_Masking/2_School_Growth.R")
source("./global_options.R")
source("./1a_Scripts/0_boottest.glm.R")
source("./1a_Scripts/0_Run_Estimators.R")
############################################# MAKE DATA ##############################################
data_in <- readRDS("./0_Data/School_Cleaned.rds") %>% filter(week <= 40) %>%
  group_by(OrgCode, OrgName) %>% arrange(week) %>%
  mutate(growth = ifelse(lag(PosPer1K)==0, PosPer1K/1, PosPer1K/lag(PosPer1K, 1)))
tname = "week";idname = "OrgCode";gname = "treat.time";control_group="notyettreated";weightsname = "total"
############################################## LOG GROWTH #############################################
yname = "growth"
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

df_sunab <- data.frame(cbind(ID, group, time, y, wt, PosPer1K = data_in$PosPer1K)) %>% 
  mutate(group = ifelse(group==0, 10000, group),
         time_to_trt = ifelse(group==10000, -1, time - group)) %>%
  filter(!is.na(y))
gweight <- sapply(glist, function(g) mean(df_sunab$group[df_sunab$group %in% glist]==g))

fit <- fepois(y ~ sunab(group, time, ref.c = 10000) | ID + time, weights = df_sunab$wt, data = df_sunab)
time_to_trt <- as.numeric(sapply(names(coef(fit)), function(var) substr(var, 7, nchar(var))))
ATT_gt <- boottest.glm(fit, M=1000, gweight=gweight, model="growth")

rownames(ATT_gt) <- time_to_trt
ATT <- colMeans(ATT_gt[as.numeric(rownames(ATT_gt)) >= 0,])

print("------------------ LOG GROWTH MODEL ------------------")
print(paste0("15-week treatment effect: ", round(exp(mean(ATT)), 2), " with CI: (",  
             round(exp(quantile(ATT, 0.025)), 2), ", ", round(exp(quantile(ATT, 0.975)), 2), ")"))
# 15-week treatment effect: 0.89 with CI: (0.7, 1.12)
2*pnorm(abs(mean(ATT)/sd(ATT)), lower.tail=FALSE) # 0.3435519

AMEs <- apply(ATT_gt, 2, function(coef) growth_AME(coef) )
print(paste0("15-week AME: ", round(growth_AME(coef=rowMeans(ATT_gt)), 1), " with CI: (",  
             round(quantile(AMEs, 0.025), 1), ", ", round(quantile(AMEs, 0.975), 1), ")"))
# 15-week AME: -343.4 with CI: (-8402.1, 88.4)
##################################################   KEEP 5 WEEKS POST INTERVENTION
ATT <- colMeans(ATT_gt[rownames(ATT_gt) %in% (0:4),])

print(paste0("5-week treatment effect: ", round(exp(mean(ATT)), 2), " with CI: (",  
             round(exp(quantile(ATT, 0.025)), 2), ", ", round(exp(quantile(ATT, 0.975)), 2), ")"))
# 5-week treatment effect: 0.97 with CI: (0.77, 1.2)
2*pnorm(abs(mean(ATT)/sd(ATT)), lower.tail=FALSE) # 0.7682157

AMEs <- apply(ATT_gt, 2, function(coef) growth_AME(coef=coef, subset=c(0:4)) )
print(paste0("5-week AME: ", round(growth_AME(coef=rowMeans(ATT_gt), subset=c(0:4)) , 1), " with CI: (",  
             round(quantile(AMEs, 0.025), 1), ", ", round(quantile(AMEs, 0.975), 1), ")"))
# 5-week AME: -7.4 with CI: (-51.2, 8.9)
