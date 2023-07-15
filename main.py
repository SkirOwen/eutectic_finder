from __future__ import annotations

import numpy as np


RELAX = 0.1         # Under - relaxation
Rg = 8.314510       # R
MAX_CONST = 20       # Max number of components
PRECISION = 1e-5    # Accuracy
To = 273.15         # Temperature reference

# Jc jacobian matrix inversion then multiply by F vector


def SolvJc(Jc, delat_x, F, n_const):
	for i in range(n_const):

		X = Jc[i, i]
		Jc[i, i] = 1.0

		for j in range(n_const):
			Jc[i, j] = Jc[i, j] / X

		for k in range(n_const):
			if k != i:				
				X = Jc[k, i]
				Jc[k, i] = 0.0
				for l in range(n_const):
					Jc[k, l] = Jc[k, l] - X * Jc[i, l] 
				
	for i in range(n_const):
		delat_x[i] = 0.0
		for j in range(n_const):
			delat_x[i] = delat_x[i] + Jc[j, i] * F[j]

	return Jc, delat_x


def eutectic(n_const: int, caract_h: np.ndarray, caract_t: np.ndarray) -> None:
	print(".....................................................................")
	print(" N Components Eutectic")
	print(" L.Brunet 2002")
	print(".....................................................................")
	print("\n")

	Jc = np.zeros((MAX_CONST, MAX_CONST), dtype=float)
	delta_x = np.zeros(MAX_CONST, dtype=float)
	X = np.zeros(MAX_CONST, dtype=float)
	F = np.zeros(MAX_CONST, dtype=float)

	print("Running...")
	# Create jacobian matrix
	for i in range(n_const):
		for j in range(n_const):
			Jc[i, j] = 0.0
		
		delta_x[i] = 0.0
		F[i] = 0.0
		X[i] = 1.0 / n_const

	X[n_const - 1] = caract_t[0]

	while True:
		for i in range(n_const):
			if i < n_const - 1:
				Jc[i, i] = 1.0 / X[i]
			else:
				coef = 1.0
				for j in range(n_const - 1):
					coef = coef - X[j]
				
				for j in range(n_const - 1):
					Jc[j, i] = -1.0 / coef

			Jc[n_const - 1, i] = -caract_h[i] / Rg / (X[n_const - 1] ** 2)

		# Create vector F
		for i in range(n_const):
			if i < n_const - 1:
				g = np.log(X[i]) + caract_h[i] / Rg * (1.0 / X[n_const - 1] - 1.0 / caract_t[i])
			else:
				g = np.log(coef) + caract_h[i] / Rg * (1.0 / X[n_const - 1] - 1.0 / caract_t[i])
			F[i] = g

		# Comprinte DeltaX
		Jc, delta_x = SolvJc(Jc, delta_x, F, n_const)
		# Newton-Raphson method
		for i in range(n_const):
			X[i] = X[i] - RELAX * delta_x[i]

		residu = 0.0
		for i in range(n_const - 1):
			residu = residu + abs(delta_x[i])
		print(residu)
		if residu == np.nan:
			break
		
		if residu < PRECISION:
			break
	# Display results
	print("-----------RESULT--------")

	coef = 1.0
	for j in range(n_const - 1):
		coef = coef - X[j]
		print(f"X {j} = {X[j]}")
	
	print(f"X {n_const} = {coef}")
	print(f"T(K)= {X[n_const - 1]} = {X[n_const - 1] - To:5.2f} øC")


def main() -> None:
	n_const = int(input("N components = "))

	# caract_h : ARRAY (1..MaxConst ) OF float;
	caract_h = np.zeros(MAX_CONST, dtype=float)
	# caract_t : ARRAY (1..MaxConst ) OF float;
	caract_t = np.zeros(MAX_CONST, dtype=float)

	for i in range(n_const):
		print(f"Enter Hfus and Tfus for {i}:")
		caract_h[i] = float(input(" H(kJ/mol) = "))
		caract_t[i] = float(input(" T(K) = "))

	# 1234
	# 1688

	eutectic(n_const, caract_h, caract_t)


if __name__ == '__main__':
	main()
