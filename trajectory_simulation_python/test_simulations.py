# Copyright © 2018, Alexander Serov


from constants import L, q, D_ratio, D_0
from D_func import D_func
import numpy as np
import unittest


class D_func_test(unittest.TestCase):
    """
    Test class for diffusivity simulations
    """

    def test_case1(self):
        # Constants
        tol = 1e-8
        # Load data
        period = L / q
        x = np.asarray([0, 0.03, 0.5, 0.6, 1.95]) * period
        # Ground truth
        D_true = D_0 + np.asarray([0.0, 0.03, 0.5, 0.4, 0.05]) * 2.0 * (D_0 * D_ratio - D_0)
        D_grad_abs_true = 2 * q * (D_ratio - 1.0) * D_0
        D_grad_true = D_grad_abs_true * np.asarray([1.0, 1.0, 1.0, -1.0, -1.0])

        # >>> D_case = 1: saw-tooth profile <<<
        D_case_number = 1
        y = 0
        D, D_grad = D_func(D_case_number, x, y, L)

        self.assertTrue(np.all(np.isclose(D_true, D, atol=tol, rtol=tol)),
                        "Simulated diffusivity failed to match the expected values. Expected:\n%s\nObtained:\n%s" % (D_true, D))
        self.assertTrue(np.all(np.isclose(D_grad_true, D_grad, atol=tol, rtol=tol)),
                        "Simulated diffusivity gradient failed to match the expected values. Expected:\n%s\nObtained:\n%s" % (D_grad_true, D_grad))


# Launch
unittest.main()
