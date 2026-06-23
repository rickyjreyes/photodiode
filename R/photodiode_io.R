# ==============================================================================
# photodiode_io.R
# Input adapters, checksums, and shared helpers for the photodiode R audit suite.
#
# All numeric .npy waveforms are read with a pure-R parser (no RcppCNPy / reticulate
# required for ordinary <f8 / <f4 / <i. arrays). reticulate is used ONLY as a
# fallback adapter for object-style ("pickled") .npy files and is never used to
# call the Python analysis pipeline.
# ==============================================================================

suppressWarnings(suppressMessages({
  library(jsonlite)
}))

# ---- Colorblind-safe palette (Okabe-Ito) -------------------------------------
PD_PALETTE <- c(
  baseline       = "#0072B2",  # blue
  post           = "#D55E00",  # vermillion
  control        = "#009E73",  # bluish green
  sham           = "#CC79A7",  # reddish purple
  disconnected   = "#56B4E9",  # sky blue
  terminated     = "#E69F00",  # orange
  dark           = "#000000",  # black
  channel_swap   = "#F0E442",  # yellow
  null           = "#999999",  # grey
  observed       = "#D55E00",
  unknown        = "#999999"
)

pd_condition_color <- function(cond) {
  cond <- tolower(as.character(cond))
  out <- PD_PALETTE[cond]
  out[is.na(out)] <- PD_PALETTE[["unknown"]]
  unname(out)
}

# ---- Output directory constants ----------------------------------------------
pd_dirs <- function(repo_root = ".") {
  list(
    repo            = normalizePath(repo_root, mustWork = FALSE),
    parity          = file.path(repo_root, "outputs_r", "parity"),
    audit_out       = file.path(repo_root, "outputs_r", "statistical_audit"),
    tables          = file.path(repo_root, "tables_r", "statistical_audit"),
    figures         = file.path(repo_root, "figures_r", "statistical_audit"),
    rendered        = file.path(repo_root, "reports", "rendered")
  )
}

pd_ensure_dirs <- function(repo_root = ".") {
  d <- pd_dirs(repo_root)
  for (p in c(d$parity, d$audit_out, d$tables, d$figures, d$rendered)) {
    if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(d)
}

# ---- Checksums ---------------------------------------------------------------
pd_checksum <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}

# ---- Pure-R .npy reader ------------------------------------------------------
# Supports 1-D and 2-D numeric arrays with little-endian dtypes <f8, <f4,
# <i8, <i4, <i2, <i1, <u1. Returns a numeric vector (1-D) or matrix (2-D).
pd_read_npy <- function(path) {
  con <- file(path, "rb"); on.exit(close(con))
  magic <- readBin(con, "raw", 6L)
  if (!identical(magic, as.raw(c(0x93, 0x4e, 0x55, 0x4d, 0x50, 0x59)))) {
    stop("Not a .npy file: ", path)
  }
  ver <- readBin(con, "raw", 2L)
  if (as.integer(ver[1]) == 1L) {
    hlen <- readBin(con, "integer", n = 1L, size = 2L, signed = FALSE, endian = "little")
  } else {
    hlen <- readBin(con, "integer", n = 1L, size = 4L, endian = "little")
  }
  hdr <- rawToChar(readBin(con, "raw", hlen))

  descr <- sub(".*'descr':\\s*'([^']+)'.*", "\\1", hdr)
  fortran <- grepl("'fortran_order':\\s*True", hdr)
  shape_str <- sub(".*'shape':\\s*\\(([^)]*)\\).*", "\\1", hdr)
  shape <- suppressWarnings(as.integer(strsplit(gsub("\\s", "", shape_str), ",")[[1]]))
  shape <- shape[!is.na(shape)]

  # object/pickled arrays cannot be parsed natively
  if (grepl("^\\|?O", descr)) {
    close(con); on.exit()
    return(pd_read_npy_reticulate(path))
  }

  endian <- if (substr(descr, 1, 1) %in% c("<", "|")) "little" else
            if (substr(descr, 1, 1) == ">") "big" else "little"
  kind <- gsub("[<>|=]", "", descr)
  typ <- substr(kind, 1, 1); sz <- as.integer(substr(kind, 2, nchar(kind)))
  n <- if (length(shape) == 0) 1L else prod(shape)

  vals <- switch(typ,
    f = readBin(con, "double",  n = n, size = sz, endian = endian),
    i = readBin(con, "integer", n = n, size = sz, signed = TRUE,  endian = endian),
    u = readBin(con, "integer", n = n, size = sz, signed = FALSE, endian = endian),
    stop("Unsupported .npy dtype: ", descr)
  )
  vals <- as.numeric(vals)

  if (length(shape) <= 1L) return(vals)
  if (length(shape) == 2L) {
    m <- if (fortran) matrix(vals, nrow = shape[1], ncol = shape[2])
         else matrix(vals, nrow = shape[1], ncol = shape[2], byrow = TRUE)
    return(m)
  }
  stop(">2-D .npy not supported: shape ", shape_str)
}

# reticulate fallback strictly as an INPUT adapter for object-style .npy
pd_read_npy_reticulate <- function(path) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Object-style .npy requires reticulate as an input adapter, which is not installed: ", path)
  }
  np <- reticulate::import("numpy", convert = TRUE)
  obj <- np$load(path, allow_pickle = TRUE)
  as.numeric(obj)
}

# ---- Waveform loader (npy or CSV) --------------------------------------------
# Returns list(t, y_raw, y, sample_rate, n, source). DC is removed in `y`.
pd_load_waveform <- function(path, sample_rate = NULL) {
  ext <- tolower(tools::file_ext(path))
  t <- NULL; y <- NULL
  if (ext == "npy") {
    obj <- pd_read_npy(path)
    if (is.matrix(obj)) {
      if (ncol(obj) >= 2) { t <- obj[, 1]; y <- obj[, 2] }
      else                 { y <- obj[, 1] }
    } else {
      y <- as.numeric(obj)
    }
  } else if (ext == "csv") {
    dat <- pd_read_waveform_csv(path)
    t <- dat$t; y <- dat$y
    if (is.null(sample_rate) && !is.null(dat$sample_rate)) sample_rate <- dat$sample_rate
  } else {
    stop("Unsupported waveform extension: ", ext)
  }

  good <- is.finite(y)
  y <- y[good]
  if (!is.null(t)) {
    t <- t[good]
    ord <- order(t); t <- t[ord]; y <- y[ord]
    dt <- stats::median(diff(t), na.rm = TRUE)
    if (is.finite(dt) && dt > 0 && is.null(sample_rate)) sample_rate <- 1 / dt
  }
  if (is.null(sample_rate)) stop("No sample rate available for ", path)
  if (is.null(t)) t <- seq(0, by = 1 / sample_rate, length.out = length(y))

  y_dc <- y - mean(y)
  list(t = t, y_raw = y, y = y_dc, sample_rate = as.numeric(sample_rate),
       n = length(y_dc), source = path)
}

# ---- Waveform CSV reader (handles time/voltage and Tektronix-style exports) --
pd_read_waveform_csv <- function(path) {
  first <- readLines(path, n = 1L, warn = FALSE)
  if (grepl("Record Length|Sample Interval|Vertical Units", first, ignore.case = TRUE)) {
    return(pd_read_tek_csv(path))
  }
  d <- utils::read.csv(path, check.names = FALSE)
  nm <- tolower(names(d))
  tcol <- which(grepl("time", nm))[1]
  vcol <- which(grepl("volt|signal|amplitude", nm))[1]
  if (is.na(vcol)) vcol <- ncol(d)
  t <- if (!is.na(tcol)) as.numeric(d[[tcol]]) else NULL
  y <- as.numeric(d[[vcol]])
  list(t = t, y = y, sample_rate = NULL)
}

# Tektronix "usr_wf_data" style: key/value header rows then a data column.
pd_read_tek_csv <- function(path) {
  raw <- readLines(path, warn = FALSE)
  si <- NA_real_
  for (ln in raw) {
    if (grepl("^Sample Interval", ln, ignore.case = TRUE)) {
      si <- suppressWarnings(as.numeric(strsplit(ln, ",")[[1]][2]))
    }
  }
  parts <- strsplit(raw, ",")
  # data rows: last numeric field is the voltage; time may be in a column too
  vals <- suppressWarnings(vapply(parts, function(p) {
    p <- p[p != ""]
    if (length(p) == 0) return(NA_real_)
    as.numeric(p[length(p)])
  }, numeric(1)))
  y <- vals[is.finite(vals)]
  sr <- if (is.finite(si) && si > 0) 1 / si else NULL
  list(t = NULL, y = y, sample_rate = sr)
}

# ---- Siglent metadata.txt parser ---------------------------------------------
pd_read_metadata_txt <- function(path) {
  if (!file.exists(path)) return(list())
  ln <- readLines(path, warn = FALSE)
  kv <- list()
  for (l in ln) {
    if (grepl("=", l)) {
      k <- sub("=.*", "", l); v <- sub("^[^=]*=", "", l)
      kv[[trimws(k)]] <- trimws(v)
    }
  }
  kv
}

# ---- Siglent raw-byte -> voltage calibration ---------------------------------
# voltage = raw_byte * vscale / 25 - offset  (Siglent SDS waveform convention)
pd_siglent_volts <- function(raw_byte, vscale, offset) {
  as.numeric(raw_byte) * vscale / 25 - offset
}

# ---- small util --------------------------------------------------------------
pd_write_json <- function(x, path) {
  writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE, null = "null",
                              digits = NA), path)
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a
