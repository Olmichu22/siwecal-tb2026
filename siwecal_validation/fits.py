"""
Fitting helpers for the validation plots.
"""

import numpy as np
from scipy.optimize import curve_fit


def gauss(x, amplitude, mean, sigma):
    """Unnormalised Gaussian used for the energy-peak fit."""
    return amplitude * np.exp(-(x - mean) ** 2 / (2 * sigma ** 2))


def fit_gaussian(values, percentiles=(5, 95), bins=50, maxfev=10000):
    """Fit a Gaussian to the core of a distribution.

    The fit uses only the ``percentiles`` window of ``values`` (to ignore tails)
    and returns
    ``(amplitude, mean, sigma, mean_err, sigma_err, bin_centers, counts)``.
    ``sigma`` is always returned positive: it only enters :func:`gauss` as
    ``sigma**2``, so the fit is bounded to ``sigma >= 0`` and the magnitude is
    taken as a guard.

    Statistical (Poisson) errors
    ----------------------------
    Each histogram bin carries a Poisson error ``sqrt(N)``. Feeding those errors
    to ``curve_fit`` together with ``absolute_sigma=True`` keeps the returned
    covariance in real statistical units (otherwise ``curve_fit`` rescales it by
    the reduced chi-square, which reflects fit quality, not the statistics). The
    1-sigma parameter uncertainties are then the square root of the covariance
    diagonal::

        perr = sqrt(diag(pcov))  ->  (amplitude_err, mean_err, sigma_err)

    Raises ``RuntimeError`` if the fit does not converge (the caller decides
    whether to skip the fit and keep the histogram).
    """
    low, high = np.percentile(values, percentiles)
    core = values[(values >= low) & (values <= high)]
    counts, edges = np.histogram(core, bins=bins)
    centers = 0.5 * (edges[:-1] + edges[1:])

    # Per-bin Poisson error sqrt(N). Empty bins carry no Poisson information, so
    # floor their weight to 1 to avoid division by zero while leaving them
    # effectively unconstrained in the fit.
    bin_errors = np.where(counts > 0, np.sqrt(counts), 1.0)

    p0 = [counts.max(), np.mean(core), np.std(core)]
    bounds = ([0, -np.inf, 0], [np.inf, np.inf, np.inf])
    popt, pcov = curve_fit(gauss, centers, counts, p0=p0, sigma=bin_errors,
                           absolute_sigma=True, bounds=bounds, maxfev=maxfev)
    amplitude, mean, sigma = popt
    # 1-sigma statistical uncertainties of (amplitude, mean, sigma).
    perr = np.sqrt(np.diag(pcov))
    mean_err, sigma_err = perr[1], perr[2]
    return amplitude, mean, abs(sigma), mean_err, sigma_err, centers, counts
