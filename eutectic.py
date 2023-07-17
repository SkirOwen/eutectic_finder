from __future__ import annotations

import numpy as np
from scipy.optimize import fsolve
from rich import print
from rich.panel import Panel
from rich.console import Console
from rich.table import Table

from utils import get_lines_for_run


def equations(variables, Hi, Ti) -> list:
	"""Setting up the equations to solve."""
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


def equations_fix_value(variables, Hi, Ti, xf: None | np.ndarray = None) -> list:
	# xi, T = variables
	xi = variables[:-1]
	T = variables[-1]
	# TODO: Hi, Ti, (and xi) needs to be len(xf) smaller.
	# TODO: Need to remove the element.s fixed by xf.

	sum_xf = 0 if xf is None else np.sum(xf)

	# Constants
	R = 8.314510  # Gas constant

	# Equations
	equations = [
		np.log(xi[i]) + Hi[i] / (R * T) - Hi[i] / (R * Ti[i]) for i in range(len(xi))
	]
	equations.append(np.sum(xi) + sum_xf - 1)

	return equations


def find_eutectic(Hi: np.ndarray, Ti: np.ndarray) -> tuple:
	"""Function to set up the initial guesses and run the solver to find the eutectic.

	Parameters
	----------
	Hi : array-like
		Array like of the enthalpy of each element, in J/mol.
	Ti : array-like
		Array like of the melting temperature of each element.

	Returns
	-------
	tuple
		the solution given by scipy.fsolve

	See Also
	--------
	scipy.fsolve
	"""
	initial_guess = np.ones(len(Hi) + 1) / (len(Hi) + 1)
	initial_guess[-1] = Ti[0]

	# Solve the equations
	solution = fsolve(equations, x0=initial_guess, args=(Hi, Ti))
	return solution


def run_eutectic(Hi: np.ndarray, Ti: np.ndarray, use_celsius: bool, xf: None | list = None) -> None:
	"""Wrapper for running the solver and printing information.

	Parameters
	----------
	Hi : array-like
		Array like of the enthalpy of each element, in J/mol.
	Ti : array-like
		Array like of the melting temperature of each element.
	use_celsius : bool
		If True the temperature is assumed to be in Celsius.
		Otherwise, using Kelvin.
	"""

	temp_unit = "C" if use_celsius else "K"
	print(Panel.fit(
		f"Enthalpies (J/mol):  {', '.join(map(str, Hi))}\n"
		f"Temperatures ({temp_unit}):    {', '.join(map(str, Ti))}",
		title="[green]Eutectic Finder",
	))
	print(f"Using Kelvin: {not use_celsius}\n")

	if xf is None:
		solution = find_eutectic(Hi, Ti)
	else:
		solution = find_eutectic(Hi, Ti)
	xi = solution[:-1]
	T = solution[-1]

	table = Table(title="Eutectic")
	table.add_column("Variable", justify="left", no_wrap=True, style="cyan")
	table.add_column("Value", justify="right", no_wrap=True, style="red")
	for i, x in enumerate(xi):
		table.add_row(
			f"x_{i}",
			f"{x:.5f}"
		)
	table.add_row(
		f"Temp Eutectic (C)",
		f"{T - 273.15:.5f}"
	)
	console = Console()
	console.print(table)


def run_from_csv(filename: str, run: int, use_celsius: bool) -> None:
	"""Function to run the solver from a csv.

	Parameters
	----------
	filename : str
		The filename of the csv to use, support path.
	run : int
		Run number to pick from the csv.
	use_celsius : bool
		If True the temperature is assumed to be in Celsius.
		Otherwise, using Kelvin.

	Raises
	------
	ValueError
		If the temperature chosen with `use_celsius` does not exist in the csv,
		i.e. has a value of -300
	"""
	run_info, mixture = get_lines_for_run(filename, run_number=run)
	print(f"Using Run #{run}")

	print(f"Elements: [green]{', '.join(map(str, run_info[:, 1]))}[/green]")

	Hi = mixture[:, 0]
	if use_celsius:
		Ti = mixture[:, 1]
		if (Ti == -300).any():
			raise ValueError("Celsius temperature not defined")
		Ti += 273.15
	else:
		Ti = mixture[:, 2]
		if (Ti == -300).any():
			raise ValueError("Kelvin temperature not defined")
	run_eutectic(Hi, Ti, use_celsius)


def main():
	run_from_csv("data.csv", 1, False)


if __name__ == "__main__":
	main()
