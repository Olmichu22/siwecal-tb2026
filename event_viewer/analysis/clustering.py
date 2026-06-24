"""
Unsupervised clustering of events on their per-event variables.

``ClusteringService`` wraps a handful of scikit-learn algorithms behind one
``fit`` method. Features are standardised first (so distance-based algorithms are
not dominated by the largest-scale variable). Rows with non-finite values (e.g.
NaN shower variables) are excluded from the fit and returned with label ``-1``,
the same sentinel DBSCAN uses for noise.
"""

from __future__ import annotations

from typing import List

import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN, KMeans, SpectralClustering
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

# Algorithm key -> human label (drives the UI dropdown).
ALGORITHMS = {
    "kmeans": "K-Means",
    "dbscan": "DBSCAN",
    "gmm": "Gaussian Mixture",
    "spectral": "Spectral",
}

UNCLUSTERED = -1


class ClusteringService:
    """Fit a clustering model on selected event-level features."""

    def _build(self, algo: str, n_clusters: int, eps: float, min_samples: int):
        if algo == "kmeans":
            return KMeans(n_clusters=n_clusters, n_init=10, random_state=0)
        if algo == "dbscan":
            return DBSCAN(eps=eps, min_samples=min_samples)
        if algo == "gmm":
            return GaussianMixture(n_components=n_clusters, random_state=0)
        if algo == "spectral":
            return SpectralClustering(n_clusters=n_clusters,
                                      assign_labels="kmeans", random_state=0)
        raise ValueError(f"Unknown clustering algorithm: {algo!r}")

    def fit(self, df: pd.DataFrame, features: List[str], algo: str,
            n_clusters: int = 3, eps: float = 0.5,
            min_samples: int = 5) -> np.ndarray:
        """Return a cluster label per row of ``df`` (``-1`` for excluded rows).

        Parameters
        ----------
        df : pandas.DataFrame
            The (already cut-filtered) per-event table.
        features : list of str
            Columns to cluster on. Must contain at least one column.
        algo : str
            One of :data:`ALGORITHMS`.
        n_clusters, eps, min_samples : numbers
            Algorithm parameters (``n_clusters`` for K-Means/GMM/Spectral,
            ``eps``/``min_samples`` for DBSCAN).
        """
        if not features:
            raise ValueError("Select at least one feature to cluster on.")

        labels = np.full(len(df), UNCLUSTERED, dtype=int)
        if len(df) == 0:
            return labels

        matrix = df[features].to_numpy(dtype=float)
        valid = np.all(np.isfinite(matrix), axis=1)
        if valid.sum() < 2:
            return labels

        scaled = StandardScaler().fit_transform(matrix[valid])
        model = self._build(algo, n_clusters, eps, min_samples)
        labels[valid] = model.fit_predict(scaled)
        return labels
