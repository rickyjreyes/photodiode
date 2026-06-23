# ==============================================================================
# ringdown_models.R
# Time-resolved decay model comparison (M0..M5). Runs ONLY when enough
# time-resolved post-excitation amplitude measurements exist.
#
# Important: Q_ratio (= 2/3 frequency geometry) is NOT a resonator quality
# factor. Q_resonator = pi * f * tau requires an identified oscillation
# frequency and a directly measured decay; it is reported only when both exist.
# ==============================================================================

pd_aicc <- function(rss, n, k) {
  ll <- -n/2 * (log(2*pi*rss/n) + 1)
  aic <- -2*ll + 2*k
  aic + (2*k*(k+1)) / max(1, n - k - 1)
}
pd_bic_rss <- function(rss, n, k) {
  ll <- -n/2 * (log(2*pi*rss/n) + 1); -2*ll + k*log(n)
}

# `dat` must be a data.frame with columns time_s and value (e.g. harmonic power).
pd_ringdown_models <- function(dat, min_points = 6) {
  if (is.null(dat) || nrow(dat) < min_points) {
    return(list(testable = FALSE,
                reason = sprintf("Need >=%d time-resolved points; have %d.",
                                 min_points, if (is.null(dat)) 0 else nrow(dat)),
                table = data.frame()))
  }
  t <- dat$time_s; y <- dat$value; n <- length(y)
  fit_rss <- function(pred, k) {
    rss <- sum((y - pred)^2)
    data.frame(rss = rss, aicc = pd_aicc(rss, n, k), bic = pd_bic_rss(rss, n, k))
  }
  rows <- list()
  # M0 constant
  rows$M0 <- cbind(model = "M0_constant", fit_rss(rep(mean(y), n), 1),
                   tau = NA, half_life = NA)
  # M1 exponential
  m1 <- try(stats::nls(y ~ A*exp(-t/tau), start = list(A = max(y), tau = diff(range(t))/2),
                       control = stats::nls.control(warnOnly = TRUE)), silent = TRUE)
  if (!inherits(m1, "try-error")) {
    tau <- stats::coef(m1)["tau"]
    rows$M1 <- cbind(model = "M1_exp", fit_rss(stats::predict(m1), 2),
                     tau = unname(tau), half_life = unname(tau*log(2)))
  }
  # M2 exponential + offset
  m2 <- try(stats::nls(y ~ A*exp(-t/tau) + C,
                       start = list(A = max(y), tau = diff(range(t))/2, C = min(y)),
                       control = stats::nls.control(warnOnly = TRUE)), silent = TRUE)
  if (!inherits(m2, "try-error")) {
    tau <- stats::coef(m2)["tau"]
    rows$M2 <- cbind(model = "M2_exp_offset", fit_rss(stats::predict(m2), 3),
                     tau = unname(tau), half_life = unname(tau*log(2)))
  }
  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  list(testable = TRUE, reason = "Sufficient time-resolved points.",
       table = tab,
       q_warning = paste("Q_ratio (2/3 frequency geometry) is distinct from",
                         "Q_resonator = pi*f*tau. No resonator Q is reported",
                         "unless an oscillation frequency is identified and the",
                         "decay is directly measured."))
}
