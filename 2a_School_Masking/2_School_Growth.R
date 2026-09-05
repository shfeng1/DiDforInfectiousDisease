here::i_am("2a_School_Masking/2_School_Growth.R")
source("./global_options.R")
source("./1a_Scripts/0_boottest.glm.R", local=TRUE)
source("./1a_Scripts/0_Run_Estimators.R", local=TRUE)
source("./2a_School_Masking/0_School_AME_Helpers.R", local=TRUE)
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

fit <- fepois(y ~ sunab(group, time, ref.c = 10000) | ID + time,
              weights=df_sunab$wt, data=df_sunab)
time_to_trt <- as.numeric(sapply(names(coef(fit)), function(var) substr(var, 7, nchar(var))))
ATT_gt <- boottest.glm(fit, gweight=gweight, model="growth")
rownames(ATT_gt) <- time_to_trt
variance.model <- school_growth_variance_model(df_sunab)

ATT_boot <- colMeans(ATT_gt[rownames(ATT_gt) >= 0,])
ATT <- mean(coef(fit)[time_to_trt >= 0])
AMEs <- school_growth_ame(ATT_gt, time_to_trt, df_sunab, variance.model)
growth_effect_15 <- c(estimate=exp(ATT),
                      lower=unname(exp(quantile(ATT_boot, 0.025))),
                      upper=unname(exp(quantile(ATT_boot, 0.975))))
growth_AME_15 <- c(estimate=mean(AMEs["AME",]),
                   lower=unname(quantile(AMEs["AME",], 0.025)),
                   upper=unname(quantile(AMEs["AME",], 0.975)))
##################################################   KEEP 5 WEEKS POST INTERVENTION
ATT_boot <- colMeans(ATT_gt[rownames(ATT_gt) %in% (0:4),])
ATT <- mean(coef(fit)[time_to_trt %in% (0:4)])
AMEs <- school_growth_ame(ATT_gt, time_to_trt, df_sunab, variance.model, subset=0:4)
growth_effect_5 <- c(estimate=exp(ATT),
                     lower=unname(exp(quantile(ATT_boot, 0.025))),
                     upper=unname(exp(quantile(ATT_boot, 0.975))))
growth_AME_5 <- c(estimate=mean(AMEs["AME",]),
                  lower=unname(quantile(AMEs["AME",], 0.025)),
                  upper=unname(quantile(AMEs["AME",], 0.975)))
