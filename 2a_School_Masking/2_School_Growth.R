rm(list=ls())
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
ATT <- colMeans(ATT_gt[rownames(ATT_gt) >= 0,])
mean(ATT) # -0.1135901
exp(mean(ATT)) # RR 0.8926238
quantile(ATT, c(0.025, 0.975)) # CI (-0.3506821, 0.1270042)
exp(quantile(ATT, c(0.025, 0.975))) # RR (0.7042076, 1.1354218)
2*pnorm(abs(mean(ATT)/sd(ATT)), lower.tail=FALSE) # 0.3553724

growth_AME(coef=rowMeans(ATT_gt)) # -317.9912
AMEs <- apply(ATT_gt, 2, function(coef) growth_AME(coef) )
quantile(AMEs, c(0.025, 0.975)) # (-6844.8307, 93.3456)
##################################################   KEEP 5 WEEKS POST INTERVENTION
ATT <- colMeans(ATT_gt[rownames(ATT_gt) %in% (0:4),]) # -0.03513041 (-0.2552712, 0.1953011), RR 0.9654795 (0.7747064, 1.2156770)
mean(ATT) # -0.03064354
exp(mean(ATT)) # RR 0.9698212
quantile(ATT, c(0.025, 0.975)) # CI (-0.2492987, 0.1919316)
exp(quantile(ATT, c(0.025, 0.975))) # RR (0.7793471, 1.2115876)
2*pnorm(abs(mean(ATT)/sd(ATT)), lower.tail=FALSE) # 0.7893808

growth_AME(coef=rowMeans(ATT_gt), subset=c(0:4)) # -7.032802
AMEs <- apply(ATT_gt, 2, function(coef) growth_AME(coef=coef, subset=c(0:4)) )
quantile(AMEs, c(0.025, 0.975)) # CI (-47.376461, 9.134382)
