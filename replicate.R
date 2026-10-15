#!/usr/bin/env Rscript
# Reproduce the simulation, application, and all manuscript tables and figures.
# See README.md for installation, run times, and the output map.

args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
script <- if (length(source_files)) tail(source_files, 1L)[[1L]] else
  if (length(file_arg)) sub("^--file=", "", file_arg[1L]) else "replicate.R"
# Rscript can encode spaces in --file; keep literal filenames when they exist.
if (!file.exists(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
package_root <- dirname(normalizePath(script, mustWork = TRUE))

if ("--help" %in% args) {
  cat(paste0(
    "Usage: Rscript replicate.R [options]\n\n",
    "  --stage=all          Run all analyses (default).\n",
    "  --stage=simulation   Baseline, few-treated, and storage experiments.\n",
    "  --stage=application  CHC application only.\n",
    "  --report-only        Use the included saved draws and estimates.\n",
    "  --quick              Eight replications; writes to results_quick/.\n",
    "  --reps=2000          Replications per simulation design.\n",
    "  --workers=4          Monte Carlo workers; default is one on Windows.\n",
    "  --output=PATH        Alternative output directory.\n",
    "  --help               Print this message.\n"))
  quit(save = "no", status = 0)
}

allowed <- "^--(report-only$|quick$|stage=(all|simulation|application)$|reps=[0-9]+$|workers=[0-9]+$|output=.+$)"
if (any(!grepl(allowed, args))) stop("Unknown option. Run with --help.", call. = FALSE)
option <- function(name, default) {
  values <- args[startsWith(args, paste0("--", name, "="))]
  if (length(values) > 1L) stop("Repeated option: ", name, call. = FALSE)
  if (length(values)) sub("^[^=]+=", "", values) else default
}

quick <- "--quick" %in% args
report_only <- "--report-only" %in% args
if (quick && report_only) stop("Choose --quick or --report-only, not both.")
stage <- option("stage", "all")
repetitions <- as.integer(option("reps", if (quick) 8L else 2000L))
workers <- as.integer(option("workers", if (.Platform$OS.type == "windows") 1L else 4L))
if (repetitions < 2L || workers < 1L || workers > 8L) stop("Invalid replication or worker count.")
if (.Platform$OS.type == "windows" && workers > 1L) stop("Use --workers=1 on Windows.")
default_output <- file.path(package_root, if (quick) "results_quick" else "results")
output <- option("output", default_output)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output, mustWork = TRUE)
config <- list(root = package_root, output = output, report_only = report_only,
  repetitions = repetitions, workers = workers,
  figures_dir = file.path(output, "figures"), tables_dir = file.path(output, "tables"))
for (path in c(config$figures_dir, config$tables_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

required <- c("ggplot2", "patchwork", "data.table")
if (!report_only) {
  required <- c(required, "Rcpp", "RcppArmadillo", "plm", "clubSandwich")
  if (stage != "simulation") {
    required <- c(required, "haven", "fixest", "etwfe", "digest", "tidyselect")
  }
}

required <- unique(required)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install these R packages before running: ", paste(missing, collapse = ", "),
    ". See README.md for the installation command.", call. = FALSE)
}

Sys.setenv(VECLIB_MAXIMUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
if (!report_only) {
  Rcpp::sourceCpp(file.path(package_root, "code", "numerics.cpp"),
    cacheDir = file.path(tempdir(), "streaming_panel_rcpp"), env = globalenv())
}

started <- proc.time()[["elapsed"]]
simulation <- application <- NULL
if (stage != "application") {
  simulation_env <- new.env(parent = globalenv())
  sys.source(file.path(package_root, "code", "simulation.R"), envir = simulation_env)
  simulation <- simulation_env$run_simulation(config)
}

if (stage != "simulation") {
  application_env <- new.env(parent = globalenv())
  sys.source(file.path(package_root, "code", "application.R"), envir = application_env)
  application <- application_env$run_application(config)
}

# These small reference tables are frozen from the manuscript's validated run.
# All calculations above are independent of them; they are used only for checks.
reference <- readRDS(file.path(package_root, "reference", "paper_results.rds"))
checks <- list()
compare <- function(label, actual, expected, tolerance = 1e-9) {
  result <- all.equal(actual, expected, tolerance = tolerance, check.attributes = FALSE)
  checks[[length(checks) + 1L]] <<- data.frame(
    check = label, passed = isTRUE(result),
    detail = if (isTRUE(result)) "Matched" else paste(result, collapse = "; "))
  if (!isTRUE(result)) stop("Manuscript benchmark failed: ", label, "\n", paste(result, collapse = "; "))
}

if (!is.null(simulation)) {
  if (simulation$baseline$configuration$repetitions == 2000L &&
    simulation$few_treated$configuration$repetitions == 2000L) {
    for (name in names(reference$baseline)) {
      compare(paste("Baseline", name), simulation$baseline$summary[[name]], reference$baseline[[name]])
    }
    for (name in names(reference$few_treated)) {
      compare(paste("Few treated", name), simulation$few_treated$summary[[name]], reference$few_treated[[name]])
    }
  } else message("Short Monte Carlo run: production benchmarks are not applied.")
}

if (!is.null(application)) {
  for (name in names(reference$application)) {
    path <- file.path(output, "application", paste0(name, ".csv"))
    actual <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    compare(paste("CHC", name), actual, reference$application[[name]])
  }
}

check_table <- if (length(checks)) do.call(rbind, checks) else
  data.frame(check = character(), passed = logical(), detail = character())
write.csv(check_table, file.path(output, "replication_checks.csv"), row.names = FALSE)
software <- data.frame(package = required,
  version = vapply(required, function(name) as.character(packageVersion(name)), character(1)))
write.csv(software, file.path(output, "software_versions.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output, "sessionInfo.txt"))
run_info <- list(stage = stage, report_only = report_only, quick = quick,
  repetitions = repetitions, workers = workers,
  elapsed_seconds = proc.time()[["elapsed"]] - started,
  R_version = R.version.string, platform = R.version$platform,
  source_md5 = setNames(unname(tools::md5sum(file.path(package_root,
  c("replicate.R", "code/simulation.R", "code/application.R", "code/numerics.cpp")))),
  c("replicate.R", "code/simulation.R", "code/application.R", "code/numerics.cpp")))
capture.output(dput(run_info), file = file.path(output, "run_info.txt"))
cat(sprintf("\nCompleted %s in %.1f seconds. %d manuscript benchmarks passed.\n",
  stage, run_info$elapsed_seconds, nrow(check_table)))
