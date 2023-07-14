from __future__ import annotations

import numpy as np


RELAX = 0.1; # Under - relaxation
Rg = 8.314510; # R
MaxConst = 20; # Max number of components
Precision = 1e-5; # Accuracy
To = 273.15; # Temperature reference

# Jc jacobian matrix inversion then multiply by F vector

def SolvJc(Jc, Dx, F, N_const):
	for i in range(N_const):

		X = Jc[i, i]
		Jc[i, i] = 1.0

		for j in range(N_const):
			Jc[i, j] = Jc[i, j] / X

		for k in range(N_const):
			if k != i:				
				X = Jc[k, i]
				Jc[k, i] = 0.0
				for l in range(N_const):
					Jc[k, l] = Jc[k, l] - X * Jc[i, l] 
				
	for i in range(N_const):
		
		Dx[i] = 0.0;
		for j in range(N_const):
			Dx[i] = Dx[i] + Jc[j, i] * F[j]   

	return Jc


def main():
	print(".....................................................................\n")
	print(" \n")
	print(" N Components Eutectic \n")
	print(" \n")
	print(" L.Brunet 2002 \n")
	print(".....................................................................\n")
	N_const = int(input("N components = "))

	Jc = np.zeros((MaxConst, MaxConst), dtype=float)
	# CaractH : ARRAY (1..MaxConst ) OF float;
	CaractH = np.zeros(MaxConst, dtype=float)
	# CaractT : ARRAY (1..MaxConst ) OF float;
	CaractT = np.zeros(MaxConst, dtype=float)

	Dx = np.zeros(MaxConst, dtype=float)
	X = np.zeros(MaxConst, dtype=float)
	F = np.zeros(MaxConst, dtype=float)

	for i in range(N_const):
		print(f"Enter Hfus and Tfus for {i}:")
		CaractH[i] = float(input(" H(kJ/mol) = "))
		CaractT[i] = float(input(" T(K) = "))

	print("Running...")
	# Create jacobian matrix
	for i in range(N_const):
		for j in range(N_const):
			Jc[i, j] = 0.0
		
		Dx[i] = 0.0
		F[i] = 0.0
		X[i] = 1.0 / N_const

	X[N_const] = CaractT[1]
	while True:
		for i in range(N_const):
			if i < N_const - 1:
				Jc[i, i] = 1.0 / X[i]
			else:
				coef = 1.0
				for j in range(N_const - 1):
					coef = coef - X[j]
				
				for j in range(N_const - 1):
					Jc[j, i] = -1.0 / coef

			Jc[N_const, i] = -CaractH[i] / Rg / (X[N_const]**2)

		# Create vector F
		for i in range(N_const):
			if i < N_const - 1:
				g = np.log(X[i]) +	CaractH[i] / Rg * (1.0 / X[N_const] - 1.0 / CaractT[i])
			else:
				g = np.log(coef) + CaractH[i] / Rg * (1.0 / X[N_const] - 1.0 / CaractT[i])

			F[i] = g
		# Comprinte DeltaX
		Jc = SolvJc(Jc, Dx, F, N_const)
		# Newton-Raphson method
		for i in range(N_const):
			X[i] = X[i] - RELAX * Dx[i]

		# Calcul du residu
		residu = 0.0
		for i in range(N_const-1):
			residu = residu + abs(Dx[i])

		print(residu, Precision)
		
		if residu < Precision:
			break
	# Display results
	print("-----------RESULT--------")

	coef = 1.0
	for j in range(N_const - 1):
		coef = coef - X[j]
		print(f"X {j} = {X[j]}")
	
	print(f"X {N_const} = {coef}")
	print(f"T(K)= {X[N_const]} = {X[N_const] - To:5.2f} øC");


if __name__ == '__main__':
	main()