


from calculate_bayes_factors import calculate_bayes_factors
from calculate_marginalized_integral import calculate_marginalized_integral
import numpy as np
import unittest

class bayes_test(unittest.TestCase):
	"""Tests for Bayes factor calculations"""

	def setUp(self):
		"""Initialize data for different tests"""
		self.dim = 2
		self.n_pi = 4
		self.tol = 1e-8
		self.B_threshold = 10.0
		self.t1_zeta_sps = np.asarray([[1.0, 1.0], [2.0, 3.0]])
		self.t1_zeta_ts = np.asarray([[2.0, 1.54], [-4.0, 0.0]])
		self.t1_ns = np.asarray([10, 5])
		self.t1_Vs = np.asarray([1, 0.2])
		self.t1_Vs_pi = np.asarray([0.3, 0.5])
		self.t1_us = np.divide(self.t1_Vs_pi, self.t1_Vs)


	def test_marginalized_integral(self):
		# Load data
		tol = self.tol
		dim = self.dim
		n_pi = self.n_pi

		# >> Test case E = 0 <<
		zeta_t = [0.7, 0.4]
		zeta_sp = [-0.3, -0.2] 
		n = 20
		u = 0.95
		eta = np.sqrt(n_pi / (n + n_pi))
		pow = dim * (n + n_pi + 1.0) / 2.0 - 2.0
		v0 = 1 + n_pi / n * u
		res = calculate_marginalized_integral(zeta_t = zeta_t, zeta_sp = zeta_sp,
			p = pow, v = v0, E = 0.0)
		self.assertTrue(res - v0 ** (-pow) < tol, "Marginalized integral calculation test failed for E = 0")


		# >> Check result with 0 break points <<
		zeta_t = [0.7, 0.4]
		zeta_sp = [-0.3, -0.2] 
		n = 20
		u = 0.95
		eta = np.sqrt(n_pi / (n + n_pi))
		pow = dim * (n + n_pi + 1.0) / 2.0 - 2.0
		v0 = 1 + n_pi / n * u
		res = calculate_marginalized_integral(zeta_t = zeta_t, zeta_sp = zeta_sp,
			p = pow, v = v0, E = eta ** 2.0)
		true_res = 0.00111261
		self.assertTrue(np.abs(res - true_res) < tol, 
			"Marginalized integral calculation test failed for 0 break points. The obtained value %.2g does not match the expected %.2g" % (res, true_res))


		# >> Check result with 1 break points <<
		zeta_t = [0.7, 0.4]
		zeta_sp = [-0.3, 0.6] 
		n = 20
		u = 0.95
		eta = np.sqrt(n_pi / (n + n_pi))
		pow = dim * (n + n_pi + 1.0) / 2.0 - 2.0
		v0 = 1 + n_pi / n * u
		res = calculate_marginalized_integral(zeta_t = zeta_t, zeta_sp = zeta_sp,
			p = pow, v = v0, E = eta ** 2.0)
		true_res = 0.00183427
		self.assertTrue(np.abs(res - true_res) < tol, 
			"Marginalized integral calculation test failed for 1 break points. The obtained value %.2g does not match the expected %.2g" % (res, true_res))


		# >> Check result with 2 break points <<
		zeta_t = [0.7, 0.4]
		zeta_sp = [0.8, 0.6] 
		n = 20
		u = 0.95
		eta = np.sqrt(n_pi / (n + n_pi))
		pow = dim * (n + n_pi + 1.0) / 2.0 - 2.0
		v0 = 1 + n_pi / n * u
		res = calculate_marginalized_integral(zeta_t = zeta_t, zeta_sp = zeta_sp,
			p = pow, v = v0, E = eta ** 2.0)
		true_res = 0.011847121
		self.assertTrue(np.abs(res - true_res) < tol, 
			"Marginalized integral calculation test failed for 2 break points. The obtained value %.8g does not match the expected %.8g" % (res, true_res))


		# >> Check the result with zeta_sp = 0 <<
		zeta_t = [0.7, 0.4]
		zeta_sp = [0.0, 0.0] 
		res = calculate_marginalized_integral(zeta_t = zeta_t, zeta_sp = zeta_sp,
			p = pow, v = v0, E = eta ** 2.0)
		true_res = 0.00246671
		self.assertTrue(np.abs(res - true_res) < tol, 
			"Marginalized integral calculation test failed for zeta_sp = 0. The obtained value %.2g does not match the expected %.2g" % (res, true_res))



	def test_bayes_factors(self):
		# Load data
		tol = self.tol
		B_threshold = self.B_threshold

		# >> Test one input bin <<
		zeta_ts = np.asarray([[0.7, 0.4]])
		zeta_sps = np.asarray([[0.8, 0.6]])
		ns = np.asarray([20])
		Vs = np.asarray([0.8 ** 2.0])
		us = np.asarray([0.95])
		Vs_pi = us * Vs
		Bs, bl_forces = calculate_bayes_factors(zeta_ts = zeta_ts, zeta_sps = zeta_sps, ns = ns, Vs = Vs, Vs_pi = Vs_pi)
		true_B = 0.3584912872575023

		# Check value
		self.assertTrue(np.isclose(Bs[0],true_B, rtol = tol), 
			"Bayes factor calculation failed for one bin. The obtained B = %.8g does not match the expected B = %.8g" % (Bs[0, 0], true_B))
		# Check force presence
		self.assertTrue((true_B >= B_threshold) == bl_forces[0], "Boolean conservative force return incorrect for the case of one bin")


		# >> Test three input bins <<
		zeta_ts = np.asarray([[0.7, 0.5], [-1.0, 2.3], [-1.1, -0.33]])
		zeta_sps = np.asarray([[0.8, 0.6], [-0.45, 0.67], [6.32, 0.115]])
		ns = np.asarray([20, 100, 3])
		Vs = np.asarray([0.8 ** 2.0, 1.25 ** 2.0, 0.3 ** 2.0])
		us = np.asarray([0.95, 8.7, 2.45])
		ns = ns.T
		Vs = Vs.T
		us = us.T
		Vs_pi = us * Vs
		N = len(ns)
		Bs, bl_forces = calculate_bayes_factors(zeta_ts = zeta_ts, zeta_sps = zeta_sps, ns = ns, Vs = Vs, Vs_pi = Vs_pi)
		true_Bs = [0.3163940070088798, 7.049757711102560e47, 1.559486019347759]

		for i in range(N):
			# Check value
			self.assertTrue(np.isclose(Bs[i], true_Bs[i], rtol = tol), 
				"Bayes factor calculation failed for one bin. For bin no. %i, the obtained B = %.8g does not match the expected B = %.8g" % (i + 1, Bs[i], true_Bs[i]))
			# Check force presence
			self.assertTrue((true_Bs[i] >= B_threshold) == bl_forces[i], 
				"For bin %i, the boolean force prediction (%r) did not correspond to expected value (%r)" % (i, true_Bs[i] >= B_threshold, bl_forces[i]))

unittest.main()