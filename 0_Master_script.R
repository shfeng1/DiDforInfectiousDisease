rm(list=ls())
here::i_am("0_Master_script.R")
#################################################################################
####################### Generate simulations from scratch #######################
### They are commented out intentionally, because it may take days on a local
###   machine to reproduce all simulations.
### source("1a_Scripts/1a_Simulate_SIR_Base_Case.R")
### source("1a_Scripts/1b_Simulate_SEIR_Base_Case.R")
### source("1a_Scripts/2a_Calculate_SIR_True_RR.R")
### source("1a_Scripts/2b_Calculate_SEIR_True_RR.R")
### source("1a_Scripts/3a_Misspecify_GI.R")
### source("1a_Scripts/3b_Misspecify_SEIR_to_SIR.R")
### source("1a_Scripts/4_Simulate_Long_Time.R")
### source("1a_Scripts/5_Simulate_Inference.R")
#######################################################################################
# Figure 1
source("1b_Summarize/1_Simulate_comparison_of_models.R")
Figure1
ggsave("4_Output/Figure 1. Summary of Models.png", width = 16, height = 10)
#######################################################################################
# Generation of Figure 2
source("1b_Summarize/2a_SIR_summ.R")
Figure2
ggsave("./4_Output/Figure 2. Simulation Results.png", width = 10, height = 6)
#######################################################################################
# Generation of Table 2
# BEGIN OF TABLE 2
source("1b_Summarize/2a_SIR_summ.R")
source("1b_Summarize/2b_SEIR_summ.R")
source("1b_Summarize/3a_Misspecify_GI_summ.R")
source("1b_Summarize/3b_Misspecify_SEIR_to_SIR_summ.R")

print("Table 2: SIR with equal transmission between groups")
tbl2.SIR1

print("Table 2: SIR with unequal transmission between groups")
tbl2.SIR2
#######################################################################################
print("Table 2: Misspecify the mean of SIR generation interval")
tbl2.misspecify.mean

print("Table 2: Misspecify the variance of SIR generation interval")
tbl2.misspecify.var
#######################################################################################
print("Table 2: Misspecify SEIR as SIR")
tbl2.misspecify.SEIR
#######################################################################################
print("Table 2: SEIR with unequal transmission between groups")
tbl2.SEIR
# END OF TABLE 2
#######################################################################################
# Generation of Table 3
# BEGIN OF TABLE 3
# Removing school mask mandates in Massachusetts
source("2a_School_Masking/3a_School_Table.R")
print("Table 3: Removing school mask mandates in Massachusetts")
print(tbl3.school)

# Mask mandates in Kansas counties
source("2b_Kansas_Masking/6a_Kansas_Table.R")
print("Table 3: Mask mandates in Kansas counties")
print(tbl3.kansas)
# END OF TABLE 3
#######################################################################################
###################################### APPENDIX #######################################
# Generation of Figures A1, A2, and A3
source("2a_School_Masking/3b_School_Graph.R")
FigureA1
ggsave("./4_Output/Figure A1. Incidence By School.png", width = 10, height = 5)

FigureA2
ggsave("./4_Output/Figure A2. Incidence Model Replication.png", width = 10, height = 5)

source("2b_Kansas_Masking/6b_Kansas_Graph.R")
FigureA3
ggsave("./4_Output/Figure A3. Incidence By Kansas County.png", width = 10, height = 5)

# Generation of Table A2
source("2a_School_Masking/4_School_Callaway_SantAnna.R")
print(Table_A2)

# Generation of Table A3
source("1b_Summarize/4_SEIR_long_time_summ.R")
print(tbl.long.time)

# Generation of Table A4
source("1b_Summarize/5_Small_N1_summ.R")
print(tbl_inference_wild)
print(tbl_inference_norm)
# END OF REPRODUCING EXHIBITS
#######################################################################################
