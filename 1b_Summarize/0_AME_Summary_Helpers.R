# Shared helpers for the simulation summary scripts.
#
# REQUIRED OUTPUT STRUCTURE
# -------------------------
# Each estimator row produced by 0_Run_Estimators.R must already contain:
#
#   Y.obs
#   Y.untrt
#   Y.untrt.adj1
#   Y.untrt.adj2
#   AME
#   AME.adj1
#   AME.adj2
#
# The true simulation AME is calculated in the same direction:
#
#   AME.true = Y.obs - Y.untrt.true

read_required_rds <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Required corrected output does not exist: ", path,
      "\nRun the corresponding corrected simulation script first."
    )
  }
  readRDS(path)
}


filter_simulation_start <- function(data, simulate_from_trt = FALSE) {
  if (!is.logical(simulate_from_trt) || length(simulate_from_trt) != 1L ||
      is.na(simulate_from_trt)) {
    stop("simulate_from_trt must be exactly TRUE or FALSE.")
  }
  if (!("simulate_from_trt" %in% names(data))) {
    return(data)
  }

  keep <- data$simulate_from_trt %in% simulate_from_trt
  if (!any(keep)) {
    stop("No simulation rows were found for simulate_from_trt=", simulate_from_trt, ".")
  }
  data[keep, , drop = FALSE]
}

assert_required_columns <- function(data, columns,
                                    object_name = deparse(substitute(data))) {
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      object_name, " is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  invisible(TRUE)
}

same_with_na <- function(x, y, tolerance = 1e-8) {
  if (length(x) != length(y)) {
    return(FALSE)
  }
  both_na <- is.na(x) & is.na(y)
  comparable <- !is.na(x) & !is.na(y)
  all(both_na | (comparable &
                   abs(x - y) <= tolerance * pmax(1, abs(x), abs(y))))
}

# Validate the stored simulation output and add only the true AME.  Fitted and
# adjusted untreated outcomes are pulled directly from the RDS file; they are
# never reverse-engineered from AME columns.
add_direct_ame_estimands <- function(data) {
  assert_required_columns(
    data,
    c(
      "Y.obs", "Y.untrt", "Y.untrt.adj1", "Y.untrt.adj2",
      "Y.untrt.true", "AME", "AME.adj1", "AME.adj2"
    ),
    "simulation output"
  )

  expected_ame <- data$Y.obs - data$Y.untrt
  expected_ame_adj1 <- data$Y.obs - data$Y.untrt.adj1
  expected_ame_adj2 <- data$Y.obs - data$Y.untrt.adj2

  if (!same_with_na(data$AME, expected_ame)) {
    stop("Stored AME is not equal to Y.obs - Y.untrt.")
  }
  if (!same_with_na(data$AME.adj1, expected_ame_adj1)) {
    stop("Stored AME.adj1 is not equal to Y.obs - Y.untrt.adj1.")
  }
  if (!same_with_na(data$AME.adj2, expected_ame_adj2)) {
    stop("Stored AME.adj2 is not equal to Y.obs - Y.untrt.adj2.")
  }

  data %>% mutate(AME.true = Y.obs - Y.untrt.true)
}
