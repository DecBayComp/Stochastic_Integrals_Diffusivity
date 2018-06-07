# The function calculates local diffusivity D in um^2/s and b'b in um/s


import math
import numpy as np
from constants import *


def D_func(D_case_number, x, y, L):

    # Load constants
    D_min = D_0
    D_max = D_0 * D_ratio

    # Initialize
    period = L / q
    D_grad_abs = 2.0 * q / L * D_0 * (D_ratio - 1.0)

    def rel_x_in_period(x):
        return (x % period) / period

    # Select f function
    def D_func_local(x):

        # Saw-tooth
        if D_case_number == 1:
            return D_max - (D_max - D_min) * np.abs(rel_x_in_period(x) - 0.5) * 2.0

        # Round well
        elif D_case_number == 2:
            return D_min + (np.sqrt((x - L / 2.0)**2 + (y - L / 2.0)**2) >= well_radius).astype(np.float16) * (D_max - D_min)

        # # Linear
        # elif D_case_number == 2:
        # 	return D_0 * (2.0 + k * x)

        # # Diffusivity jump
        # elif D_case_number == 3:
        # 	D_min = 1.0
        # 	D_max = 2.0
        # 	return D_min * int(x<0.0) + D_max * int(x>=0.0)
        # elif D_case_number == 4: # Linearly growing
        # 	D_min = 1.0
        # 	D_max = 2.0
        # 	alpha = 7.0 * math.log(2.0)
        # 	A = (D_max - D_min) / (math.exp(alpha) - 1.0)
        # 	B = D_min - A
        # 	return A * math.exp(alpha * (x + 0.5)) + B
        # elif D_case_number == 5: # Linearly growing
        # 	D_shift = 1.5
        # 	D_slope = 1.0
        # 	return D_shift + D_slope * x/float(L)

        # elif D_case_number == 6: # Oscillating
        # 	return D_0/2.0 * (2.0 + np.sin(np.pi * omega * x))

        # # Saw-tooth, two parts
        # elif D_case_number == 7:
        # 	return D_0 * (1.0 + k * (L/2.0 - np.abs(x)))

    # Select D' function

    def D_prime_func(x):
        # Saw-tooth
        if D_case_number == 1:
            sign = (rel_x_in_period(x) <= 0.5) * 2.0 - 1.0
            return D_grad_abs * sign

        # Round well
        elif D_case_number == 2:
            return x * 0.0

        # # Linear
        # elif D_case_number == 2:
        # 	return D_0 * k

        # elif D_case_number == 3:
        # 	D_min = 1.0
        # 	D_max = 2.0
        # 	return 0.0 * x
        # elif D_case_number == 4: # Linearly growing force
        # 	D_min = 1.0
        # 	D_max = 2.0
        # 	alpha = 7.0 * math.log(2.0)
        # 	A = (D_max - D_min) / (math.exp(alpha) - 1.0)
        # 	B = D_min - A
        # 	return alpha * A * math.exp(alpha * (x + 0.5))
        # elif D_case_number == 5: # Linearly growing
        # 	D_shift = 1.5
        # 	D_slope = 1.0
        # 	return D_slope * 1.0 / float(L)

        # elif D_case_number == 6: # Oscillating
        # 	return D_0/2.0 * np.pi * omega * np.cos(np.pi * omega * x)

        # # Saw-tooth, two parts
        # elif D_case_number == 7:
        # 	return - D_0 * k * np.sign(x)

    # Calculate the result
    D_func_value = D_func_local(x)
    D_prime_func_value = D_prime_func(x)

    return [D_func_value, D_prime_func_value]
