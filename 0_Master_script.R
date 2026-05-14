rm(list=ls())
here::i_am("0_Master_script.R")
#################################################################################
# Figure 1
source("1b_Summarize/1_Simulate_comparison_of_models.R")
Figure1
# ggsave("4_Output/Figure 1. Summary of Models.png", width = 16, height = 10)
#################################################################################
# Generation of Figure 2
source("1b_Summarize/2a_SIR_summ.R")
Figure2
# ggsave("./4_Output/Figure 2. Simulation Results.png", width = 10, height = 6)
#################################################################################
# Generation of Table 2
# BEGIN OF TABLE 2
source("1b_Summarize/3a_Misspecify_GI_summ.R")
source("1b_Summarize/3b_Misspecify_SEIR_to_SIR_summ.R")
source("1b_Summarize/2b_SEIR_summ.R")

print("Table 2: SIR with equal transmission between groups")
tbl2.SIR1

print("Table 2: SIR with unequal transmission between groups")
tbl2.SIR2
#################################################################################
print("Table 2: Misspecify the mean of SIR generation interval")
tbl2.misspecify.mean

print("Table 2: Misspecify the variance of SIR generation interval")
tbl2.misspecify.var
#################################################################################
print("Table 2: Misspecify SEIR as SIR")
tbl2.misspecify.SEIR
#################################################################################
print("Table 2: SEIR with unequal transmission between groups")
tbl2.SEIR
# END OF TABLE 2
#################################################################################
# Generation of Table 3
# BEGIN OF TABLE 3
# Removing school mask mandates in Massachusetts
# Note: the two scripts below take about 10 minutes to run
source("2a_School_Masking/1_School_Incidence.R")
source("2a_School_Masking/2_School_Growth.R")

# Mask mandates in in Kansas counties
# Note: the scripts take about 10 minutes with parallel computing
source("2b_Kansas_Masking/1_Kansas_Incidence.R")
source("2b_Kansas_Masking/2_Kansas_Growth.R")
source("2b_Kansas_Masking/3_Kansas_Rt.R")
source("2b_Kansas_Masking/3_Kansas_Rt_Exposure.R")
source("2b_Kansas_Masking/4_Kansas_Beta.R")
source("2b_Kansas_Masking/4_Kansas_Beta_Exposure.R")
source("2b_Kansas_Masking/5_Kansas_Beta_COVIDEstim.R")
# END OF TABLE 3
#################################################################################
################################### APPENDIX ####################################
# Generation of Figures A1, A2, and A3
source("2a_School_Masking/3_School_Graph.R")
FigureA1
ggsave("./4_Output/Figure A1. Incidence By School.png", width = 10, height = 5)

FigureA2
ggsave("./4_Output/Figure A2. Incidence Model Replication.png", width = 10, height = 5)

source("2b_Kansas_Masking/6_Kansas_Graph.R")
FigureA3
ggsave("./4_Output/Figure A3. Incidence By Kansas County.png", width = 10, height = 5)

