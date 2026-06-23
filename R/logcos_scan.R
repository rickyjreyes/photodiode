# ==============================================================================
# logcos_scan.R
# Canonical log-frequency cosine-modulation scan (PARITY_MODE) plus a
# parameterized variant used by the window/detrending sensitivity audit.
#
#   ell   = log(f)
#   y_raw = log1p(amplitude)
#   trend = quadratic poly in ell (default), removed by least squares
#   y     = standardized residual (population sd, ddof = 0)
#   M0: y = C
#   M1: y = C + a cos(k ell) + b sin(k ell)
#   delta_chi2 = RSS(M0) - RSS(M1)
# ==============================================================================

# Build detrended, standardized log-domain residual.
pd_build_log_domain <- function(spec, detrend = c("quadratic","linear","cubic")) {
  detrend <- match.arg(detrend)
  f <- spec$frequency_Hz; amp <- spec$amplitude
  ell <- log(f)
  y_raw <- log1p(amp)
  X <- switch(detrend,
    linear    = cbind(1, ell),
    quadratic = cbind(1, ell, ell^2),
    cubic     = cbind(1, ell, ell^2, ell^3))
  beta <- qr.solve(X, y_raw)
  trend <- as.numeric(X %*% beta)
  yd <- y_raw - trend
  sdpop <- sqrt(mean((yd - mean(yd))^2))      # numpy std, ddof = 0
  y <- if (sdpop > 0) yd / sdpop else yd
  data.frame(frequency_Hz = f, ell_ln_f = ell, log1p_amplitude = y_raw,
             trend = trend, y_detrended = y)
}

# Single-k constant vs constant+cos/sin least-squares fit.
pd_fit_logcos <- function(ell, y, k) {
  chi0 <- sum((y - mean(y))^2)
  c_ <- cos(k * ell); s_ <- sin(k * ell)
  X1 <- cbind(1, c_, s_)
  b1 <- qr.solve(crossprod(X1), crossprod(X1, y))
  r1 <- y - X1 %*% b1
  chi1 <- sum(r1 * r1)
  C <- b1[1]; a <- b1[2]; b <- b1[3]
  list(k = k, delta_chi2 = chi0 - chi1, chi2_null = chi0, chi2_logcos = chi1,
       C = C, a = a, b = b, A = sqrt(a^2 + b^2), phi = atan2(-b, a))
}

# Vectorized engine: closed-form constant-vs-cos/sin least squares across ALL k.
# Returns per-k C, a, b, delta_chi2 using 3x3 normal equations solved analytically.
pd_logcos_engine <- function(ell, y, kgrid) {
  m <- length(ell); sy <- sum(y); syy <- sum(y * y)
  P  <- outer(kgrid, ell)            # K x m phase matrix
  C  <- cos(P); S <- sin(P)
  A01 <- rowSums(C);      A02 <- rowSums(S)
  A11 <- rowSums(C * C);  A22 <- rowSums(S * S); A12 <- rowSums(C * S)
  r1  <- C %*% y;         r2  <- S %*% y; r0 <- sy
  # Solve [[m,A01,A02],[A01,A11,A12],[A02,A12,A22]] beta = [r0,r1,r2] per row.
  a11 <- m;  a12 <- A01; a13 <- A02
  a22 <- A11; a23 <- A12; a33 <- A22
  det <- a11*(a22*a33 - a23*a23) - a12*(a12*a33 - a23*a13) + a13*(a12*a23 - a22*a13)
  # inverse-times-rhs via Cramer's rule
  b0 <- ( (a22*a33 - a23*a23)*r0 + (a13*a23 - a12*a33)*r1 + (a12*a23 - a13*a22)*r2 ) / det
  ba <- ( (a23*a13 - a12*a33)*r0 + (a11*a33 - a13*a13)*r1 + (a13*a12 - a11*a23)*r2 ) / det
  bb <- ( (a12*a23 - a13*a22)*r0 + (a12*a13 - a11*a23)*r1 + (a11*a22 - a12*a12)*r2 ) / det
  chi0 <- syy - r0*r0/m
  chi1 <- syy - (b0*r0 + ba*as.numeric(r1) + bb*as.numeric(r2))
  list(k = kgrid, C = as.numeric(b0), a = as.numeric(ba), b = as.numeric(bb),
       chi0 = chi0, chi1 = chi1, delta = chi0 - chi1)
}

# Full scan over a k grid; returns scan data.frame and the best row.
pd_scan_logcos <- function(ell, y, kmin = 0.5, kmax = 80, nk = 4000) {
  kgrid <- seq(kmin, kmax, length.out = nk)
  e <- pd_logcos_engine(ell, y, kgrid)
  A <- sqrt(e$a^2 + e$b^2); phi <- atan2(-e$b, e$a)
  scan <- data.frame(k = kgrid, delta_chi2 = e$delta, chi2_null = e$chi0,
                     chi2_logcos = e$chi1, C = e$C, a = e$a, b = e$b, A = A, phi = phi)
  best <- scan[which.max(scan$delta_chi2), , drop = FALSE]
  list(scan = scan, best = as.list(best))
}

# Fast scan-max only (used in nulls / bootstrap).
pd_scanmax_logcos <- function(ell, y, kgrid) {
  e <- pd_logcos_engine(ell, y, kgrid)
  i <- which.max(e$delta)
  list(delta_chi2 = e$delta[i], k = kgrid[i])
}

pd_active_winding <- function(k, fmin, fmax) k * log(fmax / fmin) / (2 * pi)
pd_active_delta_ell <- function(fmin, fmax) log(fmax / fmin)

# p-value with (r+1)/(B+1) convention; never returns 0.
pd_pvalue_geq <- function(real, null_vals) {
  (1 + sum(null_vals >= real)) / (1 + length(null_vals))
}
