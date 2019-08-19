import unittest

import numpy as np
from scipy.stats import multivariate_normal as mvn

from orangecontrib.bioinformatics.utils import statistics


class TestStatistics(unittest.TestCase):
    def test_hypergeometric(self):
        threshold = 1

        # Expression matrix (cells x genes)
        x = np.array(
            [
                [1.2, 0.38, 0.14, 0.5, 0.94, 0.55, 3.88, 0.49],
                [0.26, 1.05, 1.83, 0.35, 0.19, 2.62, 0.69, 0.79],
                [0.28, 0.91, 1.78, 0.38, 0.22, 3.25, 0.66, 0.61],
                [0.26, 1.02, 1.61, 0.35, 0.2, 2.67, 0.56, 0.62],
                [0.81, 1.45, 0.48, 0.85, 3.57, 0.63, 0.16, 0.52],
                [0.78, 2.04, 0.53, 0.83, 3.05, 0.71, 0.14, 0.63],
                [0.65, 1.84, 0.54, 0.98, 3.11, 0.58, 0.16, 0.64],
                [1.34, 0.43, 0.16, 0.57, 0.73, 0.55, 4.11, 0.4],
                [1.49, 0.41, 0.15, 0.6, 0.73, 0.44, 3.49, 0.48],
                [0.29, 0.87, 1.57, 0.36, 0.18, 2.67, 0.65, 0.65],
            ]
        )

        cluster = np.zeros((x.shape[0],), dtype=bool)
        cluster[:4] = True
        scores, pvalues = statistics.score_hypergeometric_test(
            a=x[cluster], b=x[np.logical_not(cluster),], threshold=threshold
        )
        self.assertIsNotNone(scores)

    def test_alternatives(self):
        """ Test implemented alternative hypotheses. """
        np.random.seed(42)
        n = 100
        a = mvn.rvs(mean=0, cov=1, size=n).reshape((n, 1))
        b = mvn.rvs(mean=1, cov=1, size=n).reshape((n, 1))
        c = mvn.rvs(mean=-1, cov=1, size=n).reshape((n, 1))

        for method in (statistics.score_t_test, statistics.score_hypergeometric_test, statistics.score_mann_whitney):

            _, pval_less = method(a, b, alternative=statistics.ALT_LESS)
            _, pval_two = method(a, b, alternative=statistics.ALT_TWO)
            _, pval_greater = method(a, b, alternative=statistics.ALT_GREATER)
            assert np.all(pval_less < pval_two < pval_greater)
            assert np.linalg.norm(pval_less + pval_greater - np.array([1.0])) < 1e-8

            _, pval_less = method(b, c, alternative=statistics.ALT_LESS)
            _, pval_two = method(b, c, alternative=statistics.ALT_TWO)
            _, pval_greater = method(b, c, alternative=statistics.ALT_GREATER)
            assert np.all(pval_greater < pval_two < pval_less)
            assert np.linalg.norm(pval_less + pval_greater - np.array([1.0])) < 1e-8

    def test_alternatives_2d(self):
        """ Test implemented alternative hypotheses. """
        np.random.seed(42)
        n = 100
        a1 = mvn.rvs(mean=0, cov=1, size=n).reshape((n, 1))
        a2 = mvn.rvs(mean=1, cov=1, size=n).reshape((n, 1))
        b1 = mvn.rvs(mean=2, cov=1, size=n).reshape((n, 1))
        b2 = mvn.rvs(mean=3, cov=1, size=n).reshape((n, 1))
        _a = np.hstack((a1, a2))
        _b = np.hstack((b1, b2))

        for method in (statistics.score_t_test, statistics.score_hypergeometric_test, statistics.score_mann_whitney):
            _, pval_less = method(_a, _b, alternative=statistics.ALT_LESS)
            _, pval_two = method(_a, _b, alternative=statistics.ALT_TWO)
            _, pval_greater = method(_a, _b, alternative=statistics.ALT_GREATER)
            assert np.all(np.logical_and(pval_less < pval_two, pval_two < pval_greater))
            assert np.linalg.norm(pval_less + pval_greater - np.array([1.0])) < 1e-8

    def test_identical(self):
        """ Test with identical values. """
        n = 100
        a = np.random.rand(n, 1)
        b = np.random.rand(n, 1)
        _a = np.hstack([a, a])
        _b = np.hstack([a, b])
        for alt in statistics.ALTERNATIVES:
            for method in (
                statistics.score_t_test,
                statistics.score_hypergeometric_test,
                statistics.score_mann_whitney,
            ):
                p, r = method(_a, _b, alternative=alt)
                assert np.isnan(p).sum() == 0
                assert np.isnan(r).sum() == 0


if __name__ == '__main__':
    unittest.main()
