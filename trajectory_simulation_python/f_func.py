


def f_func(f_case_number, x, L):
	## Local constans
	f_shift_7 = 10.0


	## Select f function
	def f_func_local(x):
		if f_case_number == 1:
			f0 = 0.0
			return 0.0 * x + f0
		elif f_case_number == 2:
			f0 = 10.0
			return 0.0 * x + f0
		elif f_case_number == 3:
			f0 = -10.0
			return 0.0 * x + f0
		elif f_case_number == 4: # Linearly growing force 
			f_slope = 10.0
			return f_slope * (2.0*x/L)
		elif f_case_number == 5: # Linearly decreasing force
			f_slope = -5.0
			f_shift = 5.0
			return f_shift +  f_slope * (2*x/L)
		elif f_case_number == 6: # Linear
			f_slope = 1.0
			f_shift = 1.5
			return f_shift +  f_slope * (x/float(L))
		elif f_case_number == 7: # Constant
			f_shift = f_shift_7
			return f_shift +  0.0 * x	

	## Select U function: f(x) = - grad U(x)
	def U_func(x):
		if f_case_number == 1:
			f0 = 0.0
			return -f0 * x
		elif f_case_number == 2:
			f0 = 10.0
			return -f0 * x
		elif f_case_number == 3:
			f0 = -10.0
			return -f0 * x
		elif f_case_number == 4: # Linearly growing force 
			f_slope = 10.0
			return -f_slope/L * x**2
		elif f_case_number == 5: # Linearly decreasing force
			f_slope = -5.0
			f_shift = 5.0
			return -f_slope/L * x**2
		elif f_case_number == 6: # Linear
			f_shift = 1.5
			f_slope = 1.0
			return - f_shift * x - f_slope/L * x**2.0 / 2.0
		elif f_case_number == 7: # Constant
			f_shift = f_shift_7
			return - f_shift * x		


	## Calculate the result
	f_func_value = f_func_local(x);
	U_func_value = U_func(x);


	return [f_func_value, U_func_value]

