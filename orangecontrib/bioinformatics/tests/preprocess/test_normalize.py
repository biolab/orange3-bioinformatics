import unittest

import numpy as np

from Orange.data.table import Table, Domain, ContinuousVariable

from orangecontrib.bioinformatics.preprocess import QuantileNormalization


class TestQuantileNormalization(unittest.TestCase):
    def test_quantile_normalization(self):
        domain = Domain(
            [ContinuousVariable('A'), ContinuousVariable('B'), ContinuousVariable('C'), ContinuousVariable('D')]
        )
        table = Table.from_list(domain, [[5, 2, 3, 4], [4, 1, 4, 2], [3, 4, 6, 8]])
        _table = QuantileNormalization()(table)

        # expected result
        result = np.array(
            [[5.66666667, 2, 3, 4.66666667], [4.66666667, 2, 4.66666667, 3], [2, 3, 4.66666667, 5.66666667]]
        )
        np.testing.assert_array_almost_equal(result, _table.X)


if __name__ == '__main__':
    unittest.main()
