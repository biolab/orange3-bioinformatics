import numpy as np
from scipy.stats import zscore, rankdata
from sklearn.preprocessing import quantile_transform

from Orange.data.table import Table
from Orange.preprocess.preprocess import Preprocess


class LogarithmicScale(Preprocess):
    def __call__(self, data) -> Table:
        _data = data.copy()
        _data.X = np.log2(data.X + 1)
        return _data


class ZScore(Preprocess):
    """
    Compute the z score.

    Detailed description: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.zscore.html
    """

    def __init__(self, axis=0):
        self.axis = axis

    def __call__(self, data) -> Table:
        _data = data.copy()
        _data.X = zscore(data.X, axis=self.axis)
        _data.X[np.isnan(_data.X)] = 0
        return _data


class QuantileTransform(Preprocess):
    """
    Transform features to follow a uniform or a normal distribution.

    Detailed description: https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.quantile_transform.html
    """

    def __init__(self, axis=0, n_quantiles=1000, output_distribution='uniform'):
        self.axis = axis
        self.n_quantiles = n_quantiles
        self.output_distribution = output_distribution

    def __call__(self, data) -> Table:
        _data = data.copy()
        _data.X = quantile_transform(
            _data.X,
            n_quantiles=self.n_quantiles,
            output_distribution=self.output_distribution,
            copy=True,
            axis=self.axis,
        )
        return _data


class QuantileNormalization(Preprocess):
    """
    Quantile normalize a test distribution to a reference distribution
    of the same length by taking the average of each quantile across samples.

    Detailed description: https://en.wikipedia.org/wiki/Quantile_normalization
    """

    def __call__(self, data) -> Table:
        _data = data.copy()

        mean = np.mean(np.sort(_data.X, axis=1), axis=0)
        rank = rankdata(_data.X, method='average', axis=1) - 1

        rank_floor = rank.astype(int)
        rank_ceil = np.ceil(rank).astype(int)
        _data.X = (mean.take(rank_floor) + mean.take(rank_ceil)) / 2

        return _data
