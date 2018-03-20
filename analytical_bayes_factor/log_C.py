

import numpy as np
from scipy.special import gammaln	# for gammaln


def log_C(n_0, n_pi):
	res = (gammaln((n_pi - 2.0) / 2.0) 
		+ gammaln((n_0 - 3.0)/2.0) 
		- gammaln((n_pi - 3.0) / 2.0) 
		- gammaln((n_0 - 2.0) / 2.0))

	return(res)