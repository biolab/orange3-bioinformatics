import unittest

import numpy as np
import scipy.stats

from orangecontrib.bioinformatics.widgets.OWDifferentialExpression import f_oneway


class TestFOneWay(unittest.TestCase):
    def test_f_oneway(self):
        g1 = np.array([0.1, -0.1, 0.2, -0.2])
        g2 = g1 + 1
        g3 = g1

        f1, p1 = scipy.stats.f_oneway(g1, g2)
        f, p = f_oneway(g1, g2)
        np.testing.assert_almost_equal([f, p], [f1, p1])

        f, p = f_oneway(np.c_[g1], np.c_[g2], axis=0)
        np.testing.assert_almost_equal([f[0], p[0]], [f1, p1])

        f1, p1 = scipy.stats.f_oneway(g1, g2, g3)
        f, p = f_oneway(g1, g2, g3)
        np.testing.assert_almost_equal([f, p], [f1, p1])

        g1 = np.random.normal(size=(10, 30))
        g2 = np.random.normal(loc=1, size=(10, 20))
        g3 = np.random.normal(loc=2, size=(10, 10))

        f, p = f_oneway(g1, g2, g3, axis=1)
        self.assertEqual(f.shape, (10,))
        self.assertEqual(p.shape, (10,))

        fp1 = [scipy.stats.f_oneway(g1, g2, g3) for g1, g2, g3 in zip(g1, g2, g3)]

        f1 = [f for f, _ in fp1]
        p1 = [p for _, p in fp1]
        np.testing.assert_almost_equal(f1, f)
        np.testing.assert_almost_equal(p1, p)

        f, p = f_oneway(g1.T, g2.T, g3.T, axis=0)
        np.testing.assert_almost_equal(f1, f)
        np.testing.assert_almost_equal(p1, p)
