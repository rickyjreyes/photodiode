# ==============================================================================
# check_photodiode_dependencies.R
# Report required / optional package availability and environment facts.
# ==============================================================================

pd_check_dependencies <- function(verbose = TRUE) {
  required <- c("jsonlite","ggplot2","scales")
  optional <- c("patchwork","viridis","ggridges","dplyr","tidyr","data.table",
                "knitr","rmarkdown","multitaper","RcppCNPy","reticulate",
                "testthat","cowplot")
  status <- function(p) requireNamespace(p, quietly = TRUE)
  req <- data.frame(package = required, required = TRUE,
                    available = vapply(required, status, logical(1)))
  opt <- data.frame(package = optional, required = FALSE,
                    available = vapply(optional, status, logical(1)))
  tab <- rbind(req, opt)
  tab$version <- vapply(tab$package, function(p)
    if (status(p)) as.character(utils::packageVersion(p)) else NA_character_,
    character(1))
  missing_req <- tab$package[tab$required & !tab$available]
  if (verbose) {
    cat("== photodiode R audit dependency check ==\n")
    for (i in seq_len(nrow(tab))) {
      cat(sprintf("  [%s] %-12s %s%s\n",
                  ifelse(tab$available[i], "ok", if (tab$required[i]) "MISS" else "--"),
                  tab$package[i],
                  ifelse(tab$available[i], tab$version[i], ""),
                  ifelse(tab$required[i], " (required)", "")))
    }
    cat(sprintf("  pandoc: %s\n", ifelse(rmarkdown::pandoc_available(), "available", "absent")))
  }
  list(table = tab, ok = length(missing_req) == 0, missing_required = missing_req)
}

if (identical(environment(), globalenv()) &&
    !is.null(sys.calls()) && length(commandArgs(trailingOnly = TRUE)) >= 0 &&
    sys.nframe() == 0) {
  # allow `Rscript R/check_photodiode_dependencies.R`
}
