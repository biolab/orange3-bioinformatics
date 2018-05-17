import unittest
import numpy as np


from orangecontrib.bioinformatics.utils import statistics


class TestStatistics(unittest.TestCase):

    def test_hypergeometric(self):
        cluster = np.array([0, 1, 2, 3])

        treshold = 1

        # Expression matrix (cells x genes)
        X = np.array([[1.2, 0.38, 0.14, 0.5, 0.94, 0.55, 3.88, 0.49],
                      [0.26, 1.05, 1.83, 0.35, 0.19, 2.62, 0.69, 0.79],
                      [0.28, 0.91, 1.78, 0.38, 0.22, 3.25, 0.66, 0.61],
                      [0.26, 1.02, 1.61, 0.35, 0.2, 2.67, 0.56, 0.62],
                      [0.81, 1.45, 0.48, 0.85, 3.57, 0.63, 0.16, 0.52],
                      [0.78, 2.04, 0.53, 0.83, 3.05, 0.71, 0.14, 0.63],
                      [0.65, 1.84, 0.54, 0.98, 3.11, 0.58, 0.16, 0.64],
                      [1.34, 0.43, 0.16, 0.57, 0.73, 0.55, 4.11, 0.4],
                      [1.49, 0.41, 0.15, 0.6, 0.73, 0.44, 3.49, 0.48],
                      [0.29, 0.87, 1.57, 0.36, 0.18, 2.67, 0.65, 0.65]])

        scores = statistics.hypergeometric_test(X, cluster, treshold)
        self.assertIsNotNone(scores)


if __name__ == '__main__':
    unittest.main()
