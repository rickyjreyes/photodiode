# ==============================================================================
# ratio_geometry.R
# Q = 2/3 ratio-pair and Koide-style triplet geometry, plus the crucial
# arithmetic-consequence / triplet-overlap bookkeeping.
#
#   Q_low  = f1 / f2
#   Q_high = f3 / (2 f2)
#   Q_mean = (Q_low + Q_high) / 2
#   koide_error = sqrt((Q_low - 2/3)^2 + (Q_high - 2/3)^2)
# ==============================================================================

PD_KOIDE_Q <- 2/3

pd_ratio_pairs <- function(peaks) {
  fs <- peaks$frequency_Hz; amps <- peaks$amplitude
  n <- length(fs); if (n < 2) return(data.frame())
  cb <- utils::combn(n, 2)
  rows <- lapply(seq_len(ncol(cb)), function(c) {
    i <- cb[1, c]; j <- cb[2, c]
    if (fs[i] <= fs[j]) { flo <- fs[i]; fhi <- fs[j]; alo <- amps[i]; ahi <- amps[j] }
    else                { flo <- fs[j]; fhi <- fs[i]; alo <- amps[j]; ahi <- amps[i] }
    r <- flo / fhi
    data.frame(f_low = flo, f_high = fhi, ratio_low_over_high = r,
               amp_low = alo, amp_high = ahi,
               amp_ratio_low_over_high = if (ahi != 0) alo/ahi else NA,
               delta_from_2_over_3 = abs(r - 2/3),
               delta_from_1_over_2 = abs(r - 0.5),
               delta_from_3_over_4 = abs(r - 0.75))
  })
  out <- do.call(rbind, rows)
  out[order(out$delta_from_2_over_3), ]
}

pd_koide_triplets <- function(peaks) {
  d <- peaks[order(peaks$frequency_Hz), ]
  fs <- d$frequency_Hz; amps <- d$amplitude
  n <- length(fs); if (n < 3) return(data.frame())
  cb <- utils::combn(n, 3)
  rows <- lapply(seq_len(ncol(cb)), function(c) {
    i <- cb[1, c]; j <- cb[2, c]; k <- cb[3, c]
    f1 <- fs[i]; f2 <- fs[j]; f3 <- fs[k]
    q_low <- f1 / f2; q_high <- f3 / (2 * f2)
    data.frame(f1 = f1, f2 = f2, f3 = f3,
               A1 = amps[i], A2 = amps[j], A3 = amps[k],
               Q_low = q_low, Q_high = q_high, Q_mean = 0.5*(q_low+q_high),
               koide_error = sqrt((q_low - PD_KOIDE_Q)^2 + (q_high - PD_KOIDE_Q)^2),
               # exact 2:3:4 means f1:f2:f3 = 2:3:4
               ratio_234_error = sqrt((f1/f2 - 2/3)^2 + (f2/f3 - 3/4)^2))
  })
  out <- do.call(rbind, rows)
  out[order(out$koide_error), ]
}

# Triplet overlap: flag triplets that share frequencies (NOT independent).
pd_triplet_overlap <- function(triplets, koide_tol = 1e-6) {
  hit <- triplets[triplets$koide_error <= koide_tol, , drop = FALSE]
  if (!nrow(hit)) return(data.frame())
  combos <- if (nrow(hit) >= 2) utils::combn(nrow(hit), 2) else matrix(nrow = 2, ncol = 0)
  rows <- lapply(seq_len(ncol(combos)), function(c) {
    a <- hit[combos[1,c], ]; b <- hit[combos[2,c], ]
    sa <- c(a$f1, a$f2, a$f3); sb <- c(b$f1, b$f2, b$f3)
    shared <- intersect(sa, sb)
    data.frame(
      triplet_a = sprintf("(%g,%g,%g)", a$f1, a$f2, a$f3),
      triplet_b = sprintf("(%g,%g,%g)", b$f1, b$f2, b$f3),
      shared_frequencies = paste(shared, collapse = ","),
      n_shared = length(shared),
      overlapping = length(shared) > 0
    )
  })
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

# Note on arithmetic consequence: returned as a structured annotation.
pd_q_arithmetic_note <- function(triplets, koide_tol = 1e-6) {
  hit <- triplets[triplets$koide_error <= koide_tol, , drop = FALSE]
  on_ladder <- if (nrow(hit))
    all(abs(c(hit$f1, hit$f2, hit$f3) / 60 - round(c(hit$f1, hit$f2, hit$f3) / 60)) < 1e-6)
    else NA
  list(
    n_exact_q23_triplets = nrow(hit),
    all_on_60Hz_ladder = on_ladder,
    arithmetic_consequence = isTRUE(on_ladder),
    note = paste("Any integer-harmonic triplet proportional to 2:3:4 gives",
                 "Q_low = Q_high = 2/3 by construction; once a 60 Hz ladder is",
                 "established these triplets are arithmetic consequences, not",
                 "independent evidence.")
  )
}
