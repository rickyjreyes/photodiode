# ==============================================================================
# persistence_analysis.R
# Across-capture persistence vs elapsed time since excitation. Strictly gated:
# long-duration claims require multiple sequential timestamped post captures.
# ==============================================================================

# Build an elapsed-time persistence curve from a registry of post captures that
# carry a per-capture summary (e.g. 60 Hz amplitude, ladder statistic).
pd_persistence_curve <- function(capture_summaries, excitation_epoch = NULL) {
  df <- capture_summaries
  if (!"epoch" %in% names(df) || all(is.na(df$epoch))) {
    return(list(testable = FALSE, curve = data.frame(),
                reason = "No timestamped sequential post captures with epochs."))
  }
  post <- df[df$phase %in% c("post", "post-excitation", "post_excitation"), , drop = FALSE]
  if (nrow(post) < 3) {
    return(list(testable = FALSE, curve = post,
                reason = sprintf("Only %d post-excitation capture(s); >=3 needed for a curve.",
                                 nrow(post))))
  }
  e0 <- if (!is.null(excitation_epoch)) excitation_epoch else min(post$epoch, na.rm = TRUE)
  post$elapsed_s <- post$epoch - e0
  post <- post[order(post$elapsed_s), ]
  list(testable = TRUE, curve = post,
       reason = "Sequential post captures available.")
}
