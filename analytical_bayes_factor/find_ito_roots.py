


from ito_roots_equation import ito_roots_equation
import numpy as np

def find_ito_roots(n_0, n, n_pi):
	"""
	Find zeta_t roots of the equation for lambda = 0.
	There may not always exist roots...
	zeta_t = 0 is a special value for this function symmetric against zeta_t.
	It is then enough to find only the positive root.
	"""

	# First find the regions with sign changes. There must be 2 of those
	x_mesh = np.linspace(0, 10, 100)
	y = ito_roots_equation(x_mesh, n_0, n, n_pi)

	print(y)

	





