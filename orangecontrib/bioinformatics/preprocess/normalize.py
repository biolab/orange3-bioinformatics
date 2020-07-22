import numpy as np
from scipy.stats import zscore

from Orange.preprocess.preprocess import Preprocess


class LogarithmicScale(Preprocess):
    def __call__(self, data):
        _data = data.copy()
        _data.X = np.log2(data.X + 1)
        return _data


class ZScore(Preprocess):
    def __init__(self, axis=0):
        self.axis = axis

    def __call__(self, data):
        _data = data.copy()
        _data.X = zscore(data.X, axis=self.axis)
        _data.X[np.isnan(_data.X)] = 0
        return _data
