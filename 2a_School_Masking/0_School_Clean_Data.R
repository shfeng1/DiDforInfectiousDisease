rm(list=ls())
here::i_am("2a_School_Masking/0_School_Clean_Data.R")
source("./global_options.R")
source("./1a_Scripts/0_Estimate_Rt.R")

# 72 Districts included in the study 
# pulled from xxx pdf
district.list <- c("Andover", "Belmont", "Boston", "Braintree", "Brookline", "Burlington", "Cambridge", "Canton", 
                   "Carlisle", "Chelsea", "Cohasset", "Concord", "Dedham", "Dover", "Everett", "Foxborough", "Franklin", 
                   "Hanover", "Hingham", "Holbrook", "Hull", "Lexington", "Lincoln", "Lynnfield", "Malden", "Mansfield", 
                   "Marshfield", "Maynard", "Medfield", "Medford", "Medway", "Melrose", "Millis", "Milton", "Needham", 
                   "Newton", "Norfolk", "North Reading", "Norwell", "Norwood", "Quincy", "Randolph", "Reading", "Revere", 
                   "Rockland", "Saugus", "Sharon", "Somerville", "Stoneham", "Stoughton", "Sudbury", "Wakefield", "Walpole", 
                   "Waltham", "Wayland", "Wellesley", "Westwood", "Weymouth", "Wilmington", "Winchester", "Winthrop",
                   "Woburn", "Wrentham", "Acton-Boxborough", "Ayer Shirley School District", "Concord-Carlisle", 
                   "Dover-Sherborn", "Groton-Dunstable", "King Philip", "Lincoln-Sudbury", "Nashoba", "Whitman-Hanson")
one_wk <- c("Andover", "Braintree", "Burlington", "Canton", "Cohasset", "Dedham", "Dover", "Everett", "Foxborough", "Franklin", 
            "Hanover", "Hingham", "Holbrook", "Hull", "Lynnfield", "Mansfield", "Marshfield", "Maynard", "Medfield", "Medway", 
            "Millis", "Nashoba", "Norfolk", "North Reading", "Norwell", "Norwood", "Quincy", "Reading", "Rockland", "Saugus", 
            "Stoughton", "Walpole", "Waltham", "Wellesley", "Westwood", "Weymouth", "Wilmington", "Winchester", "Winthrop", 
            "Woburn", "Wrentham", "Acton-Boxborough", "Dover-Sherborn", "Groton-Dunstable", "King Philip", "Whitman-Hanson")
two_wk <- c("Belmont", "Carlisle", "Concord", "Lexington", "Malden", "Melrose", "Milton", "Needham", "Newton", "Sharon", 
            "Stoneham", "Sudbury", "Wakefield", "Wayland", "Ayer Shirley School District", "Concord-Carlisle", "Lincoln-Sudbury")
three_wk <- c("Brookline", "Cambridge", "Lincoln", "Medford", "Randolph", "Revere", "Somerville")
never <- c("Boston", "Chelsea")
###########################################################################################################################################
# merge school registrar with reported case data
df <- readRDS("./0_Data/School.rds") # source: xxx
pop.student <- as.data.frame(read_excel("./0_Data/enrollmentbygrade.xlsx", skip=1, col_names=T)) %>%
  dplyr::select("District Code", "Total") %>% `colnames<-`(c("OrgCode", "students.total"))
pop.staff <- as.data.frame(read_excel("./0_Data/teacherdata.xlsx", skip=1, col_names=T)) %>% 
  dplyr::select("District Code", "Total # of Teachers (FTE)") %>% `colnames<-`(c("OrgCode", "staffs.total"))
N <- length(unique(df$OrgName))
df.clean <- df %>% dplyr::select("OrgName", "OrgCode", "StudentsPos", "StaffPos", "RoutinePos", "RoutineTests", "end_date")
df.clean$RoutineTests <- as.numeric(gsub(",", "", df.clean$RoutineTests))
###########################################################################################################################################
##### In the following chunks, we replicated the same data pre-processing procedure as in the original paper
# https://www.nejm.org/doi/10.1056/NEJMoa2211029
###########################################################################################################################################
# Correct for 4 weeks when the schools were closed:
# The cases (two-week aggregated) are reported in the next following the closure, assuming equal incidence across both weeks
wk_missing <- c("2021-11-25", # Thanksgiving
                "2021-12-30", # Winter Recess
                "2022-02-24", # February Recess
                "2022-04-21") # Spring Recess
wk_after <- c("2021-12-02", "2022-01-06", "2022-03-03", "2022-04-28")
missing.T <- length(wk_missing)
df.add <- df %>% group_by(OrgName, OrgCode) %>% summarise(end_date = NA)
df.add <- do.call("rbind", replicate(missing.T, df.add, simplify = FALSE))
df.add$end_date <- rep(wk_missing, each=N)
df.add$wk_after <- rep(wk_after, each=N)
df.add <- merge(df.add, df.clean, by.x = c("OrgName", "OrgCode", "wk_after"), by.y = c("OrgName", "OrgCode", "end_date"))
df.add$StudentsPos <- df.add$StudentsPos/2 # assuming equal rates across both weeks
df.add$RoutinePos <- df.add$RoutinePos/2
df.add$wk_after <- NULL
df.complete <- rbind(df.clean, df.add)
df.complete$StudentsPos[df.complete$end_date %in% wk_after] <- df.complete$StudentsPos[df.complete$end_date %in% wk_after]/2
df.complete$RoutinePos[df.complete$end_date %in% wk_after] <- df.complete$RoutinePos[df.complete$end_date %in% wk_after]/2
###########################################################################################################################################
# Correct for zero-reporting
df.0report <- df.complete %>%
  group_by(OrgName, OrgCode) %>%
  arrange(OrgName, OrgCode, end_date) %>%
  mutate(week = 1:n(), school.total = StudentsPos + StaffPos,
         lag1 = lag(RoutinePos, 1), lag2 = lag(RoutinePos, 2))
df.0report$past4 <- rowMeans(df.0report[,(ncol(df.0report)-1):ncol(df.0report)], na.rm = T)
df.0report$past4[is.nan(df.0report$past4)] <- -1

# As long as there are at least a mean of 2 reported cases over the past 2 weeks but 0 cases reported this week, we flag it and set it to missing.
# Missing cases are corrected in the next step to share half of the cases from next week's reporting
# --- changes very little of the data, especially when we group to biweekly level later
df.0report$flag <- ifelse( (df.0report$school.total==0) & (df.0report$RoutinePos>0 | df.0report$past4>=2), 1, 0 )
missing.report <- df.0report %>% 
  filter(flag == 1, week != max(df.0report$week)) %>%
  dplyr::select("OrgName", "OrgCode", "week", "end_date", "flag")

df.correct <- df.0report %>% dplyr::select("OrgName", "OrgCode", "StudentsPos", "StaffPos", "RoutinePos", "RoutineTests", "end_date", "week")
for (ind in 1:nrow(missing.report)) {
  Org <- missing.report$OrgName[ind]
  wk <- missing.report$week[ind]
  df.correct$StudentsPos[df.correct$OrgName==Org & df.correct$week==wk] <- df.correct$StudentsPos[df.correct$OrgName==Org & df.correct$week==(wk+1)]/2
  df.correct$StudentsPos[df.correct$OrgName==Org & df.correct$week==(wk+1)] <- df.correct$StudentsPos[df.correct$OrgName==Org & df.correct$week==(wk+1)]/2
  
  df.correct$StaffPos[df.correct$OrgName==Org & df.correct$week==wk] <- df.correct$StaffPos[df.correct$OrgName==Org & df.correct$week==(wk+1)]/2
  df.correct$StaffPos[df.correct$OrgName==Org & df.correct$week==(wk+1)] <- df.correct$StaffPos[df.correct$OrgName==Org & df.correct$week==(wk+1)]/2
}
###########################################################################################################################################
# merge in weights
pop.student2 <- pop.student %>%
  filter(OrgCode %in% df.correct$OrgCode) %>%
  mutate(students.total = as.numeric(gsub(",","",students.total)), wt.student = students.total/sum(students.total))
pop.staff2 <- pop.staff %>%
  filter(OrgCode %in% df.correct$OrgCode) %>%
  mutate(staffs.total = as.numeric(gsub(",","",staffs.total)), wt.staff = staffs.total/sum(staffs.total))
pop.total <- merge(pop.student2, pop.staff2) %>%
  mutate(total = students.total + staffs.total, wt.total = total/sum(total))
df.clean <- merge(df.correct, pop.student2) %>% merge(pop.staff2) %>% merge(pop.total)

int.t <- min(df.clean$week[df.clean$end_date>="2022-02-28"])
df.clean <- df.clean %>%
  mutate(PosPer1K = (StudentsPos + StaffPos) / ((students.total + staffs.total)/1000),
         treat.time=ifelse(OrgName %in% one_wk, int.t+1,
                           ifelse(OrgName %in% two_wk, int.t+2,
                                  ifelse(OrgName %in% three_wk, int.t+3, 0))))
saveRDS(df.clean, "./0_Data/School_Cleaned.rds")
###########################################################################################################################################
# Bi-weekly smoothing <-- only required for generating results in the Appendix Table xxx
smooth.df <- df.clean %>%
  # filter(week > 1) %>% # change the biweek structure
  filter(week <= 40) %>% # most of school had 0 reporting in the last week
  arrange(OrgName, OrgCode, end_date) %>%
  group_by(OrgName) %>%
  mutate(biweek = rep(1:20, each = 2)) %>%
  group_by(OrgName, OrgCode, biweek) %>%
  reframe(end_date = max(end_date),
          StudentsPos = sum(StudentsPos),
          StaffPos = sum(StaffPos),
          RoutineTests = sum(RoutineTests),
          RoutinePos = sum(RoutinePos),
          students.total = mean(students.total),
          staffs.total = mean(staffs.total),
          wt.student = mean(wt.student),
          wt.staff = mean(wt.staff),
          wt.total = mean(wt.total)) %>%  
  filter(biweek > 2, biweek < 20) %>% # most zeros took place in the very first or last week
  group_by(OrgName, OrgCode) %>%
  arrange(end_date) %>%
  mutate(biweek = 1:n())

smooth.df <- smooth.df %>%
  group_by(OrgName, OrgCode, biweek, end_date) %>%
  summarise(StudentsPos = sum(StudentsPos),
            StaffPos = sum(StaffPos),
            RoutineTests = sum(RoutineTests),
            RoutinePos = sum(RoutinePos),
            students.total = sum(students.total),
            staffs.total = sum(staffs.total),
            wt.student = sum(wt.student),
            wt.staff = sum(wt.staff),
            wt.total = sum(wt.total)) %>%
  group_by(OrgName, OrgCode) %>%
  mutate(OrgCode = as.numeric(OrgCode),
         StudentsPosPer1K = StudentsPos/students.total*1000,
         StaffPosPer1K = StaffPos/staffs.total*1000)
###########################################################################################################################################
# Group small schools
# summary(smooth.df$students.total) # 1st Quartile at 2142 -- so group small schools with < 2000 students 
# Min.  1st Qu.  Median   Mean   3rd Qu.    Max. 
# 502    2142     3262    4103    4461     46001 

df.smooth.gp <- smooth.df
df.smooth.gp$OrgName[df.smooth.gp$students.total < 2000 & df.smooth.gp$OrgName %in% one_wk] <- "SmallSchools_1wk"
df.smooth.gp$OrgCode[df.smooth.gp$OrgName == "SmallSchools_1wk"] <- 99999991
df.smooth.gp$OrgName[df.smooth.gp$students.total < 2000 & df.smooth.gp$OrgName %in% c(two_wk, three_wk)] <- "SmallSchools_2or3wk"
df.smooth.gp$OrgCode[df.smooth.gp$OrgName == "SmallSchools_2or3wk"] <- 99999992

df.smooth.gp <- df.smooth.gp %>%
  group_by(OrgName, OrgCode, biweek, end_date) %>%
  summarise(StudentsPos = sum(StudentsPos),
            StaffPos = sum(StaffPos),
            RoutineTests = sum(RoutineTests),
            RoutinePos = sum(RoutinePos),
            students.total = sum(students.total),
            staffs.total = sum(staffs.total),
            wt.student = sum(wt.student),
            wt.staff = sum(wt.staff),
            wt.total = sum(wt.total)) %>%
  group_by(OrgName, OrgCode) %>%
  mutate(OrgCode = as.numeric(OrgCode),
         StudentsPosPer1K = StudentsPos/students.total*1000,
         StaffPosPer1K = StaffPos/staffs.total*1000,
         PosPer1K = (StudentsPos + StaffPos) / ((students.total + staffs.total)/1000))

unique(df.smooth.gp$OrgName[df.smooth.gp$StudentsPos==0]) # make sure there is now no school with 0 cases
unique(df.smooth.gp$OrgName[df.smooth.gp$PosPer1K==0])

saveRDS(df.smooth.gp, "School Data/df.smooth.gp.rds")
