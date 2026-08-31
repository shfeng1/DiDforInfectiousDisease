format_stata_log <- function(stata.log, max.lines=40L) {
  stata.log <- as.character(stata.log)
  stata.log <- stata.log[nzchar(trimws(stata.log))]
  if (!length(stata.log)) return("")
  stata.log <- tail(stata.log, as.integer(max.lines))
  paste0("\nStata log (last ", length(stata.log), " non-empty lines):\n",
         paste(stata.log, collapse="\n"))
}

stop_stata_condition <- function(message, stata.log=character(), class="stata_execution_error") {
  condition <- structure(
    list(message=paste0(message, format_stata_log(stata.log)),
         call=NULL, stata.log=stata.log),
    class=c(class, "error", "condition")
  )
  stop(condition)
}

my_RStata = function (src = stop("At least 'src' must be specified"), data.in = NULL,
                   data.out = FALSE, stata.path = getOption("RStata.StataPath",
                                                            stop("You need to set up a Stata path; ?chooseStataBin")),
                   stata.version = getOption("RStata.StataVersion", stop("You need to specify your Stata version")),
                   stata.echo = getOption("RStata.StataEcho", TRUE), id = "temp", ...)
{
  if (!is.character(src))
    stop("src must be a character")
  if (!(is.null(data.in) | is.data.frame(data.in)))
    stop("data.in must be NULL or a data.frame")
  if (!is.logical(data.out))
    stop("data.out must be logical")
  if (!is.numeric(stata.version))
    stop("stata.version must be numeric")
  if (!is.logical(stata.echo))
    stop("stata.echo must be logical")
  OS <- Sys.info()["sysname"]
  OS.type <- .Platform$OS.type
  SRC <- unlist(lapply(src, strsplit, "\n"))
  dataIn <- is.data.frame(data.in)
  dataOut <- data.out[1L]
  stataVersion <- stata.version[1L]
  stataEcho <- stata.echo[1L]
  if (OS %in% "Windows") {
    winRStataLog <- "RStata.log"
    on.exit(unlink(winRStataLog))
  }
  temp.root <- file.path(tempdir(), paste0("parallel_trends_stata_", Sys.getpid()))
  dir.create(temp.root, recursive=TRUE, showWarnings=FALSE)
  id <- gsub("[^A-Za-z0-9_-]", "_", as.character(id))
  doFile <- tempfile(paste0("RStata_", id, "_"), tmpdir=temp.root, fileext=".do")
  on.exit(unlink(doFile), add = TRUE)
  if (dataIn) {
    dtaInFile <- tempfile(paste0("RStataDataIn_", id, "_"), tmpdir=temp.root, fileext=".dta")
    on.exit(unlink(dtaInFile), add = TRUE)
    foreign::write.dta(data.in, file = dtaInFile, version = if (stataVersion >=
                                                                7)
      7L
      else 6L, ...)
  }
  if (dataOut) {
    dtaOutFile <- tempfile(paste0("RStataDataOut_", id, "_"), tmpdir=temp.root, fileext=".dta")
    on.exit(unlink(dtaOutFile), add = TRUE)
  }
  if (file.exists(SRC[1L]))
    SRC <- readLines(SRC[1L])
  cut_me_here <- "RSTATA: cut me here"
  cut_me_comment <- paste0("/*", cut_me_here, "*/")
  return_code_marker <- "RSTATA: return code "
  SRC <- c({
    if (dataIn) sprintf("use \"%s\"", tools::file_path_sans_ext(dtaInFile)) else ""
  }, "capture noisily {", cut_me_comment, SRC, cut_me_comment,
  "} /* end capture noisily */",
  "local RStata_rc = _rc",
  paste0("display \"", return_code_marker, "\" `RStata_rc'"),
  "if `RStata_rc' != 0 {",
  "exit `RStata_rc'",
  "}")
  SRC <- c("set more off", SRC)
  if (dataOut) {
    save_cmd <- sprintf("%s \"%s\"%s", if (stataVersion >= 13)
      "saveold"
      else "save", tools::file_path_sans_ext(dtaOutFile), if (stataVersion >=
                                                              14)
        ", version(12)"
      else "")
    SRC <- c(SRC, save_cmd)
  }
  SRC <- c(SRC, "exit, clear STATA")
  stataCmd <- paste(shQuote(stata.path), if (OS %in% "Windows")
    "/e"
    else "", "do", shQuote(doFile))
  con <- file(doFile, "w")
  writeLines(SRC, con)
  close(con)
  rdl <- pipe(stataCmd, "r")
  stataLog <- readLines(rdl, warn=FALSE)
  close(rdl)
  if (OS %in% "Windows") stataLog <- readLines(winRStataLog, warn=FALSE)
  trimmedLog <- trimws(stataLog)
  returnCodeLine <- grep(paste0("^", return_code_marker, "[0-9]+$"),
                         trimmedLog, value=TRUE)
  cutpoints <- grep(cut_me_here, stataLog)
  commandLog <- if (length(cutpoints)>=2L) {
    stataLog[seq.int(cutpoints[1L]+1L, cutpoints[2L]-1L)]
  } else {
    stataLog
  }
  if (stataEcho) {
    cat(commandLog, sep = "\n")
  }
  if (!length(returnCodeLine)) {
    stop_stata_condition(
      "Stata did not report a return code; it may have failed to start or exited unexpectedly.",
      stataLog
    )
  }
  returnCode <- as.integer(sub(return_code_marker, "", tail(returnCodeLine, 1L),
                               fixed=TRUE))
  if (!identical(returnCode, 0L)) {
    stop_stata_condition(paste0("Stata exited with return code ", returnCode, "."),
                         commandLog)
  }
  if (dataOut) {
    if (!file.exists(dtaOutFile)) {
      stop_stata_condition("Stata completed without creating the requested output data.",
                           commandLog)
    }
    res <- tryCatch(
      foreign::read.dta(dtaOutFile, ...),
      error=function(e) {
        stop_stata_condition(
          paste0("Unable to read Stata output: ", conditionMessage(e)),
          commandLog
        )
      }
    )
    attr(res, "stata.log") <- commandLog
    invisible(res)
  }
}

run_stata_bootstrap_batch <- function(data.in, models, id) {
  models <- as.character(models)
  supported <- c("inc", "loginc", "growth", "Rt_exposure", "beta_exposure")
  if (!all(models %in% supported)) stop("Unsupported bootstrap model.")
  if ("Rt_exposure" %in% models) data.in$Rt <- data.in$Rt_exposure

  model.commands <- vapply(models, function(model) {
    command <- if (model=="inc") {
      "quietly glm inc i.unit i.week i.trt_post, family(gaussian) link(identity)"
    } else if (model=="loginc") {
      "quietly glm inc i.unit i.week 1.trt_post, family(poisson) link(log)"
    } else if (model=="Rt_exposure") {
      "quietly glm Rt i.unit i.week 1.trt_post, family(poisson) link(log)"
    } else {
      paste0("quietly glm ", model,
             " i.unit i.week 1.trt_post, family(poisson) link(log)")
    }
    paste(
      "set seed 123456789",
      command,
      "quietly boottest 1.trt_post, cluster(unit) reps(10000) quietly",
      paste0("post `results' (`s') (\"", model, "\") (r(p))"),
      sep="\n"
    )
  }, character(1))

  command <- paste(
    "set matsize 1000",
    "tempname results",
    "tempfile resultdata",
    "postfile `results' long replication_id str20 model float p using `resultdata'",
    "levelsof replication_id, local(simulations)",
    "foreach s of local simulations {",
    "preserve",
    "keep if replication_id==`s'",
    "drop replication_id",
    paste(model.commands, collapse="\n"),
    "restore",
    "}",
    "postclose `results'",
    "use `resultdata', clear",
    sep="\n"
  )

  out <- my_RStata(command, data.in=data.in, data.out=TRUE,
                   stata.echo=FALSE, id=id)
  stata.log <- attr(out, "stata.log", exact=TRUE)
  required <- c("replication_id", "model", "p")
  if (!all(required %in% names(out))) {
    stop_stata_condition(
      paste0("Stata bootstrap batch ", id, " did not return columns ",
             paste(setdiff(required, names(out)), collapse=", "), "."),
      stata.log, class="stata_bootstrap_error"
    )
  }
  out$replication_id <- as.integer(out$replication_id)
  out$model <- as.character(out$model)
  expected <- expand.grid(
    replication_id=sort(unique(data.in$replication_id)),
    model=models, stringsAsFactors=FALSE
  )
  key <- paste(out$replication_id, out$model, sep="::")
  expected.key <- paste(expected$replication_id, expected$model, sep="::")
  if (anyDuplicated(key) || !setequal(key, expected.key)) {
    missing.key <- setdiff(expected.key, key)
    duplicate.key <- unique(key[duplicated(key)])
    detail <- c(
      if (length(missing.key)) paste0("missing: ", paste(head(missing.key, 10L), collapse=", ")),
      if (length(duplicate.key)) paste0("duplicates: ", paste(head(duplicate.key, 10L), collapse=", "))
    )
    stop_stata_condition(
      paste0("Stata bootstrap batch ", id,
             " returned incomplete or duplicate results",
             if (length(detail)) paste0(" (", paste(detail, collapse="; "), ")") else "",
             "."),
      stata.log, class="stata_bootstrap_error"
    )
  }
  out
}

is_retryable_stata_bootstrap_error <- function(e) {
  known.write.error <- exists("is_allowed_write_error", mode="function") &&
    is_allowed_write_error(e)
  known.write.error || inherits(e, "stata_execution_error") ||
    inherits(e, "stata_bootstrap_error")
}

run_stata_bootstrap_attempts <- function(data.in, models, id, attempts=2L) {
  attempts <- max(1L, as.integer(attempts))
  last.error <- NULL
  for (attempt in seq_len(attempts)) {
    value <- tryCatch(
      run_stata_bootstrap_batch(data.in, models, id),
      error=function(e) e
    )
    if (!inherits(value, "error")) return(value)
    last.error <- value
    if (!is_retryable_stata_bootstrap_error(value) || attempt==attempts) break
    first.line <- strsplit(conditionMessage(value), "\n", fixed=TRUE)[[1L]][1L]
    message("Retrying Stata bootstrap batch ", id, " after: ", first.line)
    Sys.sleep(0.5*attempt + (Sys.getpid() %% 5L)/10)
  }
  stop(last.error)
}

run_stata_bootstrap_with_fallback <- function(data.in, models, id, attempts=2L) {
  batch <- tryCatch(
    run_stata_bootstrap_attempts(data.in, models, id, attempts),
    error=function(e) e
  )
  if (!inherits(batch, "error")) return(batch)
  if (!is_retryable_stata_bootstrap_error(batch)) stop(batch)

  message("Retrying Stata bootstrap batch ", id,
          " one replication at a time after the batch retry failed.")
  pieces <- lapply(sort(unique(data.in$replication_id)), function(replication.id) {
    one <- data.in[data.in$replication_id==replication.id, , drop=FALSE]
    single <- tryCatch(
      run_stata_bootstrap_attempts(
        one, models, paste0(id, "_", replication.id), attempts
      ),
      error=function(e) e
    )
    if (!inherits(single, "error")) return(single)
    if (exists("is_allowed_write_error", mode="function") &&
        is_allowed_write_error(single)) return(NULL)
    stop_stata_condition(
      paste0("Stata bootstrap failed for replication ", replication.id,
             " after retrying the batch and the individual replication:\n",
             conditionMessage(single)),
      class="stata_bootstrap_error"
    )
  })
  data.table::rbindlist(pieces, use.names=TRUE, fill=TRUE)
}

run_simulation_chunk <- function(indices, rng.streams, simulate.one,
                                 models, id, p.column="p") {
  pieces <- lapply(indices, function(s) {
    capture_simulation_error({
      set_simulation_stream(rng.streams[[s]])
      value <- simulate.one(s)
      value$replication_id <- as.integer(s)
      value
    })
  })
  pieces <- Filter(Negate(is.null), pieces)
  if (!length(pieces)) return(NULL)

  bootstrap.data <- data.table::rbindlist(lapply(pieces, function(piece) {
    transform(piece$bootstrap_data, replication_id=piece$replication_id)
  }), use.names=TRUE, fill=TRUE)
  bootstrap <- run_stata_bootstrap_with_fallback(bootstrap.data, models, id)
  if (!nrow(bootstrap)) return(NULL)

  complete.ids <- bootstrap %>%
    dplyr::count(replication_id) %>%
    dplyr::filter(n==length(models)) %>%
    dplyr::pull(replication_id)
  results <- data.table::rbindlist(lapply(pieces, function(piece) {
    transform(piece$result, replication_id=piece$replication_id)
  }), use.names=TRUE, fill=TRUE)
  results <- results[results$replication_id %in% complete.ids, , drop=FALSE]
  key <- paste(results$replication_id, results$model, sep="::")
  bootstrap.key <- paste(bootstrap$replication_id, bootstrap$model, sep="::")
  results[[p.column]] <- bootstrap$p[match(key, bootstrap.key)]
  results$replication_id <- NULL
  results
}

run_pure_simulation_chunk <- function(indices, rng.streams, simulate.one) {
  pieces <- lapply(indices, function(s) {
    capture_simulation_error({
      set_simulation_stream(rng.streams[[s]])
      simulate.one(s)
    })
  })
  pieces <- Filter(Negate(is.null), pieces)
  if (!length(pieces)) return(NULL)
  data.table::rbindlist(pieces, use.names=TRUE, fill=TRUE)
}
