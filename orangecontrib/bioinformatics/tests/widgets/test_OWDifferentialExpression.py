import unittest
import scipy.stats
import numpy as np

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

        G1 = np.random.normal(size=(10, 30))
        G2 = np.random.normal(loc=1, size=(10, 20))
        G3 = np.random.normal(loc=2, size=(10, 10))

        F, P = f_oneway(G1, G2, G3, axis=1)
        self.assertEqual(F.shape, (10,))
        self.assertEqual(P.shape, (10,))

        FP1 = [scipy.stats.f_oneway(g1, g2, g3)
               for g1, g2, g3 in zip(G1, G2, G3)]

        F1 = [f for f, _ in FP1]
        P1 = [p for _, p in FP1]
        np.testing.assert_almost_equal(F1, F)
        np.testing.assert_almost_equal(P1, P)

        F, P = f_oneway(G1.T, G2.T, G3.T, axis=0)
        np.testing.assert_almost_equal(F1, F)
        np.testing.assert_almost_equal(P1, P)