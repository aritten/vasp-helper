import numpy as np

from vasphelper import math_functions as mf


def test_gaussian_convolution_returns_expected_shape_and_values():
    E_pts = np.array([-1.0, 1.0])
    sigma = 1.0

    E_grid, I_smooth = mf.gaussian_convolution(E_pts, sigma)

    expected_grid = np.linspace(E_pts.min(), E_pts.max(), 600)
    expected_kernel = np.exp(-((expected_grid[:, None] - E_pts[None, :]) ** 2) / (2 * sigma**2))
    expected_I_smooth = expected_kernel @ np.ones_like(E_pts)

    assert E_grid.shape == (600,)
    assert I_smooth.shape == (600,)
    assert np.allclose(E_grid, expected_grid)
    assert np.allclose(I_smooth, expected_I_smooth)
    assert np.all(I_smooth > 0)