#!/usr/bin/env Rscript

# Thin process entrypoint. Statistical functions remain under R/, while the
# Electron protocol, analysis commands, and visualization commands live in
# focused worker modules.

script_arguments <- enc2utf8(commandArgs(trailingOnly = FALSE))
script_flag <- grep("^--file=", script_arguments, value = TRUE)
script_path <- if (length(script_flag) > 0L) {
    sub("^--file=", "", script_flag[[1]])
} else {
    file.path(getwd(), "inst", "worker", "analysis_worker.R")
}
worker_root <- dirname(normalizePath(script_path, mustWork = TRUE))

for (module in c("runtime.R", "mageck.R", "analysis_commands.R", "visualization_commands.R", "dispatcher.R")) {
    sys.source(file.path(worker_root, "modules", module), envir = .GlobalEnv)
}

run_worker_loop()
