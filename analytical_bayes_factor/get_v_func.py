


def get_v_func(ns, n_pi, us, dim, zeta_t_perp):

	def v(s):
		return(1.0 + n_pi / ns * us + (dim - 1.0) * s * zeta_t_perp ** 2)
	return v