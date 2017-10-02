## The function calculates local diffusivity D in um^2/s and b'b in um/s


import math
import numpy as np


def D_func(D_case_number, x, L):

	## Select f function
	def D_func_local(x):
		if D_case_number == 1:
			D_min = 1.0
			D_max = 2.0
			return D_min + (D_max - D_min) * (2.0*x/L)**2
		elif D_case_number == 2:
			D_min = 1.0
			D_max = 2.0
			return D_min + (x/L + 0.5) * (D_max - D_min)
		elif D_case_number == 3:
			D_min = 1.0
			D_max = 2.0
			return D_min * int(x<0.0) + D_max * int(x>=0.0)
		elif D_case_number == 4: # Linearly growing 
			D_min = 1.0
			D_max = 2.0
			alpha = 7.0 * math.log(2.0)
			A = (D_max - D_min) / (math.exp(alpha) - 1.0)
			B = D_min - A
			return A * math.exp(alpha * (x + 0.5)) + B
		elif D_case_number == 5: # Linearly growing
			D_shift = 1.5
			D_slope = 1.0
			return D_shift + D_slope * x/float(L)
		elif D_case_number == 6: # Oscillating
			D_0 = 1.0e-2	# in um^2/s
			w = 10.0	# in um^-1
			return D_0/2.0 * (2 + np.sin(np.pi * w * x))
			


	## Select D' function
	def D_prime_func(x):
		if D_case_number == 1:
			D_min = 1.0
			D_max = 2.0
			return (D_max - D_min) * 8.0 * x / L**2
		elif D_case_number == 2:
			D_min = 1.0
			D_max = 2.0
			return (D_max - D_min)/L
		elif D_case_number == 3:
			D_min = 1.0
			D_max = 2.0
			return 0.0 * x
		elif D_case_number == 4: # Linearly growing force 
			D_min = 1.0
			D_max = 2.0
			alpha = 7.0 * math.log(2.0)
			A = (D_max - D_min) / (math.exp(alpha) - 1.0)
			B = D_min - A
			return alpha * A * math.exp(alpha * (x + 0.5))
		elif D_case_number == 5: # Linearly growing
			D_shift = 1.5
			D_slope = 1.0
			return D_slope * 1.0 / float(L)
		elif D_case_number == 6: # Oscillating
			D_0 = 1.0e-2	# in um^2/s
			w = 10.0	# in um^-1
			return D_0/2.0 * np.pi * w * np.cos(np.pi * w * x)


	## Calculate the result
	D_func_value = D_func_local(x);
	D_prime_func_value = D_prime_func(x);


	return [D_func_value, D_prime_func_value]







