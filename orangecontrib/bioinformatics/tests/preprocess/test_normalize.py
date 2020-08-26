import unittest

import numpy as np

from Orange.data.table import Table, Domain, ContinuousVariable

from orangecontrib.bioinformatics.preprocess import ZScore, LogarithmicScale, QuantileNormalization


class TestQuantileNormalization(unittest.TestCase):
    def test_quantile_normalization(self):
        domain = Domain(
            [ContinuousVariable('A'), ContinuousVariable('B'), ContinuousVariable('C'), ContinuousVariable('D')]
        )
        table = Table.from_list(domain, [[5, 2, 3, 4], [4, 1, 4, 2], [3, 4, 6, 8]])
        _table = QuantileNormalization()(table)

        # expected result
        result = np.array([[5.66666667, 2, 3, 4.66666667], [5.166667, 2, 5.166667, 3], [2, 3, 4.66666667, 5.66666667]])
        np.testing.assert_array_almost_equal(result, _table.X)


class TestZScore(unittest.TestCase):
    def test_z_score(self):
        domain = Domain([ContinuousVariable('A'), ContinuousVariable('B')])
        table = Table.from_list(domain, [[1, 2], [3, 4]])

        _table = ZScore(axis=0)(table)
        np.testing.assert_array_almost_equal([[-1, -1], [1, 1]], _table.X)

        _table = ZScore(axis=1)(table)
        np.testing.assert_array_almost_equal([[-1, 1], [-1, 1]], _table.X)


class TestLogarithmicScale(unittest.TestCase):
    def test_log2_normalization(self):
        domain = Domain([ContinuousVariable('A'), ContinuousVariable('B')])
        table = Table.from_list(domain, [[0, 1], [3, 7], [15, 31]])

        _table = LogarithmicScale()(table)
        np.testing.assert_array_almost_equal([[0, 1], [2, 3], [4, 5]], _table.X)


if __name__ == '__main__':
    unittest.main()
