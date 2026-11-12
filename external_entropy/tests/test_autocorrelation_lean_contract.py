"""Exact 2x2 DFT and autocorrelation source correspondence fixtures."""
import numpy as np
import pytest
from wct_entropy.diagnostics import centered_autocorrelation, coherence_profile


@pytest.mark.parametrize("field", [
    np.array([[1., 2.], [3., 4.]]),
    np.array([[1., -2.], [-3., 4.]]),
    np.array([[1., 0.], [0., 0.]]),
])
def test_fft2_power_inverse_matches_periodic_2x2_formula(field):
    a, b, c, d = field.ravel()
    F = np.fft.fft2(field)
    expected_F = np.array([[a+b+c+d, a-b+c-d],
                           [a+b-c-d, a-b-c+d]])
    np.testing.assert_allclose(F, expected_F)
    power = np.abs(F)**2
    energy = a*a + b*b + c*c + d*d
    np.testing.assert_allclose(np.sum(power)/4, energy)
    corr = np.fft.ifft2(power).real
    np.testing.assert_allclose(corr, [[energy, 2*(a*b+c*d)],
                                      [2*(a*c+b*d), 2*(a*d+b*c)]])
    shifted = centered_autocorrelation(field)
    np.testing.assert_allclose(shifted[1, 1], 1.0, atol=1e-14)
    np.testing.assert_allclose(shifted[0, 0], 2*(a*d+b*c)/energy, atol=1e-14)
    np.testing.assert_allclose(coherence_profile(field)[0], 1.0, atol=1e-14)


@pytest.mark.parametrize("shape", [(2,2),(3,5),(4,6),(7,4)])
def test_fftshift_center_and_normalization_any_parity(shape):
    field = np.random.default_rng(314).normal(size=shape)
    ac = centered_autocorrelation(field)
    np.testing.assert_allclose(ac[shape[0]//2,shape[1]//2],1.0,atol=1e-13)
    np.testing.assert_allclose(ac, np.fft.fftshift(np.fft.ifft2(
        np.abs(np.fft.fft2(field))**2).real) / np.sum(field**2),atol=1e-13)


def test_zero_field_rejected():
    with pytest.raises(ValueError, match="non-positive zero-lag"):
        centered_autocorrelation(np.zeros((2,2)))
