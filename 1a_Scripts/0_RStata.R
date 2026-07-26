# Run Stata from R using a unique worker-local temporary directory.
#
# This implementation is safe for foreach parallel execution and for direct
# Kansas-analysis calls:
#   * every invocation uses unique input, output, do-file, and log paths;
#   * Stata is launched in the operating system's documented batch/exit mode;
#   * the Stata command return code is written to and parsed from an explicit
#     sentinel in the Stata log;
#   * the wrapper does not mistake a launcher exit status for a model failure;
#   * p-values are returned as data rather than scraped from printed output.

rstata_stop <- function(stage, message, log = character()) {
  log <- log[nzchar(trimws(log))]
  log.text <- if (length(log)) {
    paste0(
      "\n--- Stata log (last 120 lines) ---\n",
      paste(tail(log, 120L), collapse = "\n")
    )
  } else {
    ""
  }

  stop(
    sprintf("my_RStata failed during %s: %s%s", stage, message, log.text),
    call. = FALSE
  )
}

resolve_stata_executable <- function(stata.path) {
  if (length(stata.path) != 1L || is.na(stata.path) || !nzchar(stata.path)) {
    rstata_stop(
      "Stata configuration",
      paste0(
        "RStata.StataPath is not set. Set it in global_options.R before ",
        "running the analysis or registering the foreach backend."
      )
    )
  }

  stata.path <- path.expand(stata.path)
  # Tolerate an option value stored with surrounding quotes.
  stata.path <- sub('^"(.*)"$', '\\1', stata.path)

  if (file.exists(stata.path)) {
    return(normalizePath(stata.path, winslash = "/", mustWork = TRUE))
  }

  resolved <- Sys.which(stata.path)
  if (nzchar(resolved)) {
    return(normalizePath(resolved, winslash = "/", mustWork = TRUE))
  }

  rstata_stop(
    "Stata configuration",
    paste0("Stata executable does not exist and is not on PATH: ", stata.path)
  )
}

stata_quote_path <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  paste0('"', gsub('"', '""', path, fixed = TRUE), '"')
}

format_stata_number <- function(x) {
  if (length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop("The Stata null-hypothesis value must be one finite number.", call. = FALSE)
  }

  # Keep enough significant digits for the fine CI-inversion grids while
  # preventing scientific notation from entering the Stata hypothesis string.
  trimws(format(x, scientific = FALSE, digits = 17L))
}

is_allowed_transient_write_error <- function(message) {
  write.failure <-
    grepl("unable to open file for writing", message, fixed = TRUE) ||
    grepl("cannot open the connection", message, fixed = TRUE)

  transient.nfs <-
    grepl("Stale NFS file handle", message, fixed = TRUE) ||
    grepl("Operation timed out", message, fixed = TRUE)

  write.failure && transient.nfs
}

handle_simulation_error <- function(e, parameter_row, replicate) {
  message <- conditionMessage(e)

  # Preserve only the two explicitly allowed transient file-system failures.
  if (is_allowed_transient_write_error(message)) {
    return(NULL)
  }

  stop(
    sprintf(
      "Simulation parameter row %s, replicate %s failed: %s",
      as.character(parameter_row), as.character(replicate), message
    ),
    call. = FALSE
  )
}

my_RStata <- local({
  # Capture these values when the file is sourced on the master. The closure
  # carries them to PSOCK workers whose process-local options may be empty.
  default.stata.path <- getOption("RStata.StataPath", NULL)
  default.stata.version <- getOption("RStata.StataVersion", NULL)
  default.stata.echo <- getOption("RStata.StataEcho", TRUE)

  function(
      src = stop("At least 'src' must be specified"),
      data.in = NULL,
      data.out = FALSE,
      stata.path = default.stata.path,
      stata.version = default.stata.version,
      stata.echo = default.stata.echo,
      id = "temp",
      ...
  ) {
    if (!is.character(src) || length(src) < 1L) {
      stop("src must be a non-empty character vector", call. = FALSE)
    }
    if (!(is.null(data.in) || is.data.frame(data.in))) {
      stop("data.in must be NULL or a data.frame", call. = FALSE)
    }
    if (!is.logical(data.out) || length(data.out) < 1L || is.na(data.out[1L])) {
      stop("data.out must be logical", call. = FALSE)
    }
    if (!is.numeric(stata.version) || length(stata.version) < 1L || is.na(stata.version[1L])) {
      rstata_stop(
        "Stata configuration",
        paste0(
          "RStata.StataVersion is not set. Set it in global_options.R before ",
          "running the analysis or registering the foreach backend."
        )
      )
    }
    if (!is.logical(stata.echo) || length(stata.echo) < 1L || is.na(stata.echo[1L])) {
      stop("stata.echo must be logical", call. = FALSE)
    }

    stata.exec <- resolve_stata_executable(stata.path)
    stata.version <- stata.version[1L]
    stata.echo <- stata.echo[1L]
    dataIn <- is.data.frame(data.in)
    dataOut <- isTRUE(data.out[1L])
    OS <- unname(Sys.info()["sysname"])

    safe.id <- gsub("[^A-Za-z0-9_.-]", "_", as.character(id)[1L])
    work.dir <- tempfile(
      pattern = paste0("rstata_", safe.id, "_"),
      tmpdir = tempdir()
    )

    created <- dir.create(work.dir, recursive = TRUE, showWarnings = FALSE)
    if (!created || !dir.exists(work.dir)) {
      rstata_stop("temporary-directory creation", work.dir)
    }

    keep.on.error <- isTRUE(getOption("RStata.KeepTempOnError", FALSE))
    success <- FALSE
    on.exit({
      if (success || !keep.on.error) {
        unlink(work.dir, recursive = TRUE, force = TRUE)
      }
    }, add = TRUE)

    doFile <- file.path(work.dir, "run.do")
    dtaInFile <- file.path(work.dir, "data_in.dta")
    dtaOutFile <- file.path(work.dir, "data_out.dta")
    commandLog <- file.path(work.dir, "stata_command.log")
    launcherLog <- file.path(work.dir, "stata_launcher.log")
    automaticBatchLog <- file.path(work.dir, "run.log")

    if (dataIn) {
      tryCatch(
        foreign::write.dta(
          data.in,
          file = dtaInFile,
          version = if (stata.version >= 7) 7L else 6L,
          ...
        ),
        error = function(e) {
          rstata_stop(
            "writing the worker-local Stata input file",
            conditionMessage(e)
          )
        }
      )

      if (!file.exists(dtaInFile)) {
        rstata_stop("writing the worker-local Stata input file", dtaInFile)
      }
    }

    SRC <- unlist(
      lapply(src, strsplit, split = "\n", fixed = TRUE),
      use.names = FALSE
    )

    if (length(SRC) == 1L && file.exists(SRC[1L])) {
      SRC <- tryCatch(
        readLines(SRC[1L], warn = FALSE),
        error = function(e) {
          rstata_stop(
            "reading the supplied Stata source file",
            conditionMessage(e)
          )
        }
      )
    }

    cut.me.here <- "RSTATA: cut me here"
    cut.me.comment <- paste0("/*", cut.me.here, "*/")
    rc.sentinel <- "RSTATA_RETURN_CODE="

    commands <- c(
      "set more off",
      "capture log close _all",
      paste0(
        "log using ", stata_quote_path(commandLog),
        ", text replace name(rstata)"
      )
    )

    if (dataIn) {
      commands <- c(
        commands,
        paste0("use ", stata_quote_path(dtaInFile), ", clear")
      )
    }

    commands <- c(
      commands,
      "capture noisily {",
      cut.me.comment,
      SRC,
      cut.me.comment,
      "}",
      "local rstata_rc = _rc"
    )

    if (dataOut) {
      save.command <- if (stata.version >= 14) {
        paste0(
          "saveold ", stata_quote_path(dtaOutFile),
          ", version(12) replace"
        )
      } else if (stata.version >= 13) {
        paste0("saveold ", stata_quote_path(dtaOutFile), ", replace")
      } else {
        paste0("save ", stata_quote_path(dtaOutFile), ", replace")
      }

      # Preserve the model/boottest return code. Only attempt to save when the
      # requested Stata commands succeeded; if saving fails, return that code.
      commands <- c(
        commands,
        "if `rstata_rc' == 0 {",
        paste0("    capture noisily ", save.command),
        "    local rstata_rc = _rc",
        "}"
      )
    }

    commands <- c(
      commands,
      paste0("display as text \"", rc.sentinel, "`rstata_rc'\""),
      "log close rstata"
    )

    tryCatch(
      writeLines(commands, con = doFile, useBytes = TRUE),
      error = function(e) {
        rstata_stop(
          "writing the worker-local Stata do-file",
          conditionMessage(e)
        )
      }
    )

    if (!file.exists(doFile)) {
      rstata_stop("writing the worker-local Stata do-file", doFile)
    }

    old.wd <- getwd()
    on.exit(setwd(old.wd), add = TRUE)
    tryCatch(
      setwd(work.dir),
      error = function(e) {
        rstata_stop(
          "entering the worker-local Stata directory",
          conditionMessage(e)
        )
      }
    )

    # NECESSARY FIX [STATA-BATCH]: v5 passed only "do <file>" on macOS/Linux.
    # Stata's documented shell interface requires a batch option. On macOS,
    # -e also closes the GUI automatically; on Unix/Linux, -b is the documented
    # batch mode; Windows uses /e.
    batch.option <- if (identical(OS, "Windows")) {
      "/e"
    } else if (identical(OS, "Darwin")) {
      "-e"
    } else {
      "-b"
    }

    args <- c(batch.option, "do", shQuote(doFile))
    launcher.status <- tryCatch(
      system2(
        command = stata.exec,
        args = args,
        stdout = launcherLog,
        stderr = launcherLog,
        wait = TRUE
      ),
      error = function(e) {
        rstata_stop("starting the Stata process", conditionMessage(e))
      }
    )

    read_log <- function(path) {
      if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0L) {
        return(character())
      }
      tryCatch(readLines(path, warn = FALSE), error = function(e) character())
    }

    command.log <- read_log(commandLog)
    launcher.log <- read_log(launcherLog)
    batch.log <- read_log(automaticBatchLog)
    combined.log <- c(command.log, batch.log, launcher.log)

    # The explicit sentinel is the authoritative Stata command result. This is
    # intentionally separate from the operating-system launcher status, which
    # can be nonzero for GUI/batch-launch reasons even after the do-file ran.
    rc.lines <- grep(
      paste0(rc.sentinel, "[0-9]+"),
      combined.log,
      value = TRUE
    )

    if (!length(rc.lines)) {
      rstata_stop(
        "Stata execution",
        paste0(
          "Stata did not write its return-code sentinel. Launcher status was ",
          as.character(launcher.status),
          ". This usually means Stata could not start the do-file or could not ",
          "open its explicit command log."
        ),
        combined.log
      )
    }

    stata.rc <- suppressWarnings(as.integer(sub(
      paste0(".*", rc.sentinel, "([0-9]+).*"),
      "\\1",
      tail(rc.lines, 1L)
    )))

    if (is.na(stata.rc)) {
      rstata_stop(
        "parsing the Stata return code",
        tail(rc.lines, 1L),
        combined.log
      )
    }

    if (stata.echo && length(command.log)) {
      cutpoints <- grep(cut.me.here, command.log, fixed = TRUE)
      if (length(cutpoints) >= 2L) {
        cat(
          command.log[seq.int(cutpoints[1L] + 1L, cutpoints[2L] - 1L)],
          sep = "\n"
        )
      } else {
        cat(command.log, sep = "\n")
      }
    }

    if (stata.rc != 0L) {
      rstata_stop(
        "Stata commands",
        paste0("Stata returned r(", stata.rc, ")."),
        combined.log
      )
    }

    if (dataOut) {
      if (!file.exists(dtaOutFile)) {
        rstata_stop(
          "reading Stata output",
          paste0(
            "Stata returned r(0) but did not create the requested ",
            "data_out.dta file. Launcher status was ",
            as.character(launcher.status), "."
          ),
          combined.log
        )
      }

      result <- tryCatch(
        foreign::read.dta(dtaOutFile, ...),
        error = function(e) {
          rstata_stop(
            "reading the worker-local Stata output file",
            conditionMessage(e),
            combined.log
          )
        }
      )

      success <- TRUE
      return(invisible(result))
    }

    success <- TRUE
    invisible(NULL)
  }
})

# Run one wild-cluster bootstrap test and return r(p) directly as a scalar.
# This replaces all fragile capture.output()/last-line parsing in the Kansas
# scripts. Any Stata error is handled by my_RStata before this function returns.
run_stata_boottest_p <- function(
    model_command,
    hypothesis,
    cluster,
    data.in,
    reps = 10000L,
    id = "boottest"
) {
  if (!is.character(model_command) || length(model_command) != 1L || !nzchar(model_command)) {
    stop("model_command must be one non-empty Stata command string.", call. = FALSE)
  }
  if (!is.character(hypothesis) || length(hypothesis) != 1L || !nzchar(hypothesis)) {
    stop("hypothesis must be one non-empty boottest hypothesis.", call. = FALSE)
  }
  if (!is.character(cluster) || length(cluster) != 1L || !nzchar(cluster)) {
    stop("cluster must name one Stata cluster variable.", call. = FALSE)
  }
  if (!is.numeric(reps) || length(reps) != 1L || is.na(reps) || reps < 1) {
    stop("reps must be one positive integer.", call. = FALSE)
  }

  command <- paste(
    model_command,
    paste0(
      "boottest ", hypothesis,
      ", cluster(", cluster, ") reps(", as.integer(reps), ") quietly"
    ),
    "gen double p = r(p) in 1",
    "keep p",
    "keep if _n == 1",
    sep = "\n"
  )

  out <- my_RStata(
    src = command,
    data.in = data.in,
    data.out = TRUE,
    stata.echo = FALSE,
    id = id
  )

  if (!is.data.frame(out) || !"p" %in% names(out) || nrow(out) < 1L) {
    stop("Stata completed but did not return the expected p variable.", call. = FALSE)
  }

  p <- suppressWarnings(as.numeric(out$p[1L]))
  if (length(p) != 1L || is.na(p) || !is.finite(p) || p < 0 || p > 1) {
    stop("Stata returned an invalid wild-bootstrap p-value.", call. = FALSE)
  }

  p
}
