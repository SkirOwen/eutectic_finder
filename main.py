import numpy as np
from scipy.optimize import fsolve
from time import time, monotonic_ns


def equations(variables):
	# xi, T = variables
	xi = variables[:-1]
	T = variables[-1]

	# Constants
	R = 8.314510  # Gas constant

	# Equations
	equations = [
		np.log(xi[i]) + Hi[i] / (R * T) - Hi[i] / (R * Ti[i]) for i in range(len(xi))
	]
	equations.append(np.sum(xi) - 1)

	return equations


# Set the values for Hi and Ti
Hi = np.array([11300, 46000])  # Example values
Ti = np.array([1234, 1688])  # Example values

Hi = np.array([25500, 13400, 28500])  # Example values
Ti = np.array([772, 610, 808]) + 273.15  # Example values

# Hi = np.array([11700, 25500, 15900])  # Example values
# Ti = np.array([337, 254, 310]) + 273.15  # Example values

# Initial guess for xi and T
# initial_guess = np.array([[0.5, 0.5], [1234, 1688]])
start = time()
initial_guess = np.ones(len(Hi) + 1) / (len(Hi) + 1)
initial_guess[-1] = Ti[0]

# Solve the equations
solution = fsolve(equations, x0=initial_guess)
stop = time()
print(stop - start)

xi = solution[:-1]
T = solution[-1]

print("xi:", xi)
print("T:", T - 273.15)
