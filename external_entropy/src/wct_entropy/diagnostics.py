from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

Array = NDArray[np.float64]


def normalized_density(field: Array) -> Array:
    """Normalize |psi|^2 to a discrete probability distribution."""

    density = np.abs(field) ** 2
    total = float(np.sum(density))
    if not np.isfinite(total) or total <= 0.0:
        raise ValueError("field has zero or non-finite intensity")
    return density / total


def shannon_entropy(field: Array) -> float:
    """Discrete spatial Shannon entropy of normalized field intensity."""

    probability = normalized_density(field)
    nonzero = probability > 0.0
    return float(-np.sum(probability[nonzero] * np.log(probability[nonzero])))


def inverse_participation_ratio(field: Array) -> float:
    """Return sum(p^2); larger values indicate stronger concentration."""

    probability = normalized_density(field)
    return float(np.sum(probability**2))


def centered_autocorrelation(field: Array) -> Array:
    """Return a real, fft-shifted autocorrelation normalized at zero lag."""

    spectrum = np.fft.fft2(field)
    correlation = np.fft.ifft2(np.abs(spectrum) ** 2).real
    correlation = np.fft.fftshift(correlation)
    center = (field.shape[0] // 2, field.shape[1] // 2)
    peak = float(correlation[center])
    if not np.isfinite(peak) or peak <= 0.0:
        raise ValueError("autocorrelation has a non-positive zero-lag value")
    return correlation / peak


def coherence_profile(field: Array) -> Array:
    """Average the positive-lag center row and column autocorrelations."""

    correlation = centered_autocorrelation(field)
    cy, cx = correlation.shape[0] // 2, correlation.shape[1] // 2
    max_radius = min(cx, cy)
    row = correlation[cy, cx : cx + max_radius]
    column = correlation[cy : cy + max_radius, cx]
    return 0.5 * (row + column)


def coherence_length(
    field: Array,
    dx: float,
    dy: float,
    threshold: float = np.exp(-1.0),
) -> float:
    """First positive lag where the axis-averaged correlation crosses threshold.

    Linear interpolation is used between samples. If no crossing occurs, the
    largest sampled lag is returned.
    """

    if not 0.0 < threshold < 1.0:
        raise ValueError("threshold must lie strictly between zero and one")
    profile = coherence_profile(field)
    spacing = 0.5 * (dx + dy)
    crossings = np.flatnonzero(profile <= threshold)
    if crossings.size == 0:
        return float((len(profile) - 1) * spacing)
    index = int(crossings[0])
    if index == 0:
        return 0.0
    y0, y1 = float(profile[index - 1]), float(profile[index])
    if y0 == y1:
        fraction = 1.0
    else:
        fraction = (y0 - threshold) / (y0 - y1)
    return float((index - 1 + np.clip(fraction, 0.0, 1.0)) * spacing)


def l2_norm(field: Array, dx: float, dy: float) -> float:
    return float(np.sqrt(np.sum(field**2) * dx * dy))


def gradient_energy(field: Array, dx: float, dy: float) -> float:
    grad_x = (np.roll(field, -1, axis=1) - np.roll(field, 1, axis=1)) / (2.0 * dx)
    grad_y = (np.roll(field, -1, axis=0) - np.roll(field, 1, axis=0)) / (2.0 * dy)
    return float(0.5 * np.sum(grad_x**2 + grad_y**2) * dx * dy)
