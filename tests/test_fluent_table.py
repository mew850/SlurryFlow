import unittest

import numpy as np

from src.export.fluent_table import interpolate_table


class FluentTableTests(unittest.TestCase):
    def test_exact_nodes_and_log_bilinear_midpoint(self):
        phi = np.array([0.0, 1.0])
        log_gamma = np.array([0.0, 1.0])
        log_eta = np.array([[0.0, 1.0], [1.0, 2.0]])
        self.assertAlmostEqual(interpolate_table(0, 1, phi, log_gamma, log_eta), 1.0)
        self.assertAlmostEqual(interpolate_table(1, 10, phi, log_gamma, log_eta), 100.0)
        self.assertAlmostEqual(interpolate_table(0.5, 10**0.5, phi, log_gamma, log_eta), 10.0)

    def test_out_of_range_is_clamped(self):
        phi = np.array([0.1, 0.2]); log_gamma = np.array([0.0, 1.0])
        log_eta = np.array([[0.0, 0.0], [1.0, 1.0]])
        self.assertAlmostEqual(interpolate_table(-1, 1e-9, phi, log_gamma, log_eta), 1.0)
        self.assertAlmostEqual(interpolate_table(2, 1e9, phi, log_gamma, log_eta), 10.0)


if __name__ == "__main__":
    unittest.main()
