#!/usr/bin/env Rscript

# Shared helpers for CKM deliverable pipelines.

log_step <- function(msg) {
  cat(sprintf("\n[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}

assert_file <- function(path) {
  if (!file.exists(path)) stop(sprintf("Required script/file not found: %s", path), call. = FALSE)
}

run_script <- function(path) {
  assert_file(path)
  log_step(sprintf("Running: %s", path))
  source(path, chdir = TRUE, encoding = "UTF-8")
  log_step(sprintf("Done: %s", path))
}
