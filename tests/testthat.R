library(testthat)
repo <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), ".."), mustWork = FALSE)
if (!dir.exists(file.path(repo, "R"))) repo <- "."
`%||%` <- function(a,b) if (is.null(a)) b else a
Sys.setenv(PD_REPO_ROOT = repo)
test_dir("tests/testthat", reporter = "summary")
