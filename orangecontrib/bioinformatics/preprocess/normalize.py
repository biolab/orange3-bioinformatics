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
    def __init__(self, axis=0):
        self.axis = axis

    def __call__(self, data) -> Table:
        _data = data.copy()
        _data.X = zscore(data.X, axis=self.axis)
        _data.X[np.isnan(_data.X)] = 0
        return _data


class QuantileTransform(Preprocess):
    def __init__(self, axis=0, n_quantiles=100, output_distribution='uniform'):
        self.axis = axis
        self.n_quantiles = n_quantiles
        self.output_distribution = output_distribution

    def __call__(self, data) -> Table:
        _data = data.copy()
        _data.X = quantile_transform(
            _data.X, n_quantiles=self.n_quantiles, output_distribution=self.output_distribution, copy=True
        )
        return _data


class QuantileNormalization(Preprocess):
    def __call__(self, data) -> Table:
        _data = data.copy()
        sorted_indices = np.argsort(_data.X, axis=1)
        mean = np.mean(np.take_along_axis(_data.X, sorted_indices, axis=1), axis=0)
        rank = rankdata(_data.X, method='min', axis=1) - 1
        _data.X = mean.take(rank)
        return _data
