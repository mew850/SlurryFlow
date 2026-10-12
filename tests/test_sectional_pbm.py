import unittest

import numpy as np

from src.models.sectional_pbm import (
    Config,
    relative_viscosity,
    relative_viscosity_limits,
    run_case,
)


class ReproductionTests(unittest.TestCase):
    def test_viscosity_limits_and_interpolation(self):
        phi = np.array([1e-8, 0.05, 0.2])
        m0, minf = relative_viscosity_limits(phi)
        self.assertTrue(np.all(m0 >= 1.0) and np.all(minf >= 1.0))
        self.assertTrue(np.allclose(relative_viscosity(phi, 0), m0))
        self.assertTrue(np.allclose(relative_viscosity(phi, 1e20), minf))

    def test_population_mass_is_conserved(self):
        cfg = Config(n_bins=30, agg_scale=1e13, break_scale=0.5, t_end=0.1)
        model, _, (_, _, _, _, c) = run_case(cfg)
        mass = np.sum(c * model.mass[:, None], axis=0)
        self.assertLess(np.max(np.abs(mass - 1.0)), 1e-10)

    def test_expected_direction_of_structure_change(self):
        _, _, (_, df_a, _, _, _) = run_case(Config(n_bins=30, agg_scale=1e13, t_end=0.1))
        _, _, (_, df_s, _, _, _) = run_case(Config(n_bins=30, shear_compact=2.0, t_end=0.1))
        self.assertLess(df_a[-1], df_a[0])
        self.assertGreater(df_s[-1], df_s[0])


if __name__ == "__main__":
    unittest.main()
