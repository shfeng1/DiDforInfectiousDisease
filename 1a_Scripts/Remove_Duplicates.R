here::i_am("1a_Scripts/Remove_Duplicates.R")

output.file <- "./4_Output/SIR_N1=5.rds"
clean.file  <- output.file

sim.out <- readRDS(output.file)

# Each row in sim.out is one p_out; require equality across every p_out column.
p_out.names <- names(sim.out)

dup <- duplicated(sim.out[, p_out.names, drop=FALSE])

cat("Original rows:   ", nrow(sim.out), "\n")
cat("Duplicate rows:  ", sum(dup), "\n")
cat("Rows after clean:", sum(!dup), "\n")

sim.out.clean <- sim.out[!dup, , drop=FALSE]
rownames(sim.out.clean) <- NULL
# Verify that no exact duplicates remain across all p_out columns.
stopifnot(!any(duplicated(sim.out.clean[, p_out.names, drop=FALSE])))

saveRDS(sim.out.clean, clean.file)
