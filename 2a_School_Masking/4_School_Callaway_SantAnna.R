rm(list=ls())
here::i_am("2a_School_Masking/4_School_Callaway_SantAnna.R")
source("./global_options.R")
source("./1a_Scripts/0_attgt.glm.R")
###########################################################################################################################################
# Callaway & Sant'Anna estimates replication in Table A2
###########################################################################################################################################
df.clean <- readRDS("./0_Data/School_Cleaned.rds")
data_in <- df.clean %>% filter(week <= 40)
###########################################################################################################################################
# Incidence
inc <- attgt.glm(yname="PosPer1K", tname="week", idname="OrgCode", gname="treat.time",
                 weightsname="total", data_in=data_in, family=gaussian(link="identity"))
inc.ATT <- sum(inc)
inc.AME <- inc.ATT
###########################################################################################################################################
# Log incidence
loginc <- attgt.glm(yname="PosPer1K", tname="week", idname="OrgCode", gname="treat.time",
                    weightsname="total", data_in=data_in, family=poisson(link="log"))
loginc.ATT <- sum(loginc)

loginc.AME <- data_in %>%
  group_by(OrgName) %>%
  arrange(week) %>%
  filter(week >= min(data_in$treat.time[data_in$treat.time > 0])) %>%
  mutate(Pos.ctl = PosPer1K / exp(loginc),
         diff = PosPer1K - Pos.ctl) %>%
  filter(treat.time > 0, week >= treat.time) %>%
  group_by(week) %>%
  summarise(diff=mean(diff), .groups="drop") %>%
  summarise(AME=sum(diff)) %>%
  pull(AME)
###########################################################################################################################################
# Log growth
data_in_growth <- data_in %>%
  group_by(OrgCode, OrgName) %>%
  arrange(week) %>%
  mutate(growth=ifelse(lag(PosPer1K) == 0, PosPer1K, PosPer1K / lag(PosPer1K, 1)))

loggrowth <- attgt.glm(yname="growth", tname="week", idname="OrgCode", gname="treat.time",
                       weightsname="total", data_in=data_in_growth, family=poisson(link="log"))
loggrowth.ATT <- sum(loggrowth)

growth.obs <- data_in_growth %>%
  group_by(OrgName) %>%
  arrange(week) %>%
  filter(week >= min(data_in$treat.time[data_in$treat.time > 0]) - 1) %>%
  mutate(growth.ctl = log(growth) - c(0, loggrowth),
         Pos.ctl = NA_real_) %>%
  filter(treat.time > 0, week >= treat.time - 1, !is.na(growth))

for (unit in unique(growth.obs$OrgName)) {
  for (time in sort(unique(growth.obs$week[growth.obs$OrgName == unit]))) {
    if ((time + 1) == unique(growth.obs$treat.time[growth.obs$OrgName == unit])) {
      growth.obs$Pos.ctl[growth.obs$OrgName == unit & growth.obs$week == time] <-
        growth.obs$PosPer1K[growth.obs$OrgName == unit & growth.obs$week == time]
    } else {
      time.ind <- ifelse(nrow(growth.obs[growth.obs$OrgName == unit & growth.obs$week == time - 1,]) == 0, 2, 1)
      growth.obs$Pos.ctl[growth.obs$OrgName == unit & growth.obs$week == time] <-
        growth.obs$Pos.ctl[growth.obs$OrgName == unit & growth.obs$week == time - time.ind] *
        exp(growth.obs$growth.ctl[growth.obs$OrgName == unit & growth.obs$week == time])
    }
  }
}

loggrowth.AME <- growth.obs %>%
  filter(week >= treat.time) %>%
  mutate(diff=PosPer1K - Pos.ctl) %>%
  group_by(OrgName) %>%
  summarise(diff=sum(diff), .groups="drop") %>%
  summarise(AME=mean(diff)) %>%
  pull(AME)

###########################################################################################################################################
# Table A2
###########################################################################################################################################
Table_A2 <- tibble(`Outcome specification` = c("Incidence", "Log incidence", "Log growth"),
  `Treatment effect` = c(inc.ATT, exp(loginc.ATT), exp(loggrowth.ATT)),
  `Average marginal effect` = c(inc.AME, loginc.AME, loggrowth.AME)) %>%
  mutate(across(where(is.numeric), ~format(round(.x, 1), nsmall=1)))
