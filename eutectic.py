import numpy as np
from scipy.optimize import fsolve
from rich import print
from rich.panel import Panel
from rich.console import Console
from rich.table import Table

from utils import get_lines_for_run


def equations(variables, Hi, Ti) -> list:
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


def find_eutectic(Hi, Ti) -> tuple:
	initial_guess = np.ones(len(Hi) + 1) / (len(Hi) + 1)
	initial_guess[-1] = Ti[0]

	# Solve the equations
	solution = fsolve(equations, x0=initial_guess, args=(Hi, Ti))
	return solution


def run_eutectic(Hi, Ti, use_celsius) -> None:
	"""Wrapper for printing information"""

	temp_unit = "C" if use_celsius else "K"
	print(Panel.fit(
		f"Enthalpies (J/mol):  {', '.join(map(str, Hi))}\n"
		f"Temperatures ({temp_unit}):    {', '.join(map(str, Ti))}",
		title="[green]Eutectic Finder",
	))
	print(f"Using Kelvin: {not use_celsius}\n")

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


def run_from_csv(filename: str, run: int, use_celsius: bool):
	run_info, mixture = get_lines_for_run(filename, run_number=run)
	print(f"Using Run #{run}")

	print(f"Elements: [green]{', '.join(map(str, run_info[:, 1]))}[/green]")

	Hi = mixture[:, 0]
	if use_celsius:
		Ti = mixture[:, 1]
		if (Ti == -1).any():
			raise ValueError("Celsius temperature not defined")
		Ti += 273.15
	else:
		Ti = mixture[:, 2]
		if (Ti == -1).any():
			raise ValueError("Kelvin temperature not defined")
	run_eutectic(Hi, Ti, use_celsius)


def main():
	run_from_csv("data.csv", 1, False)


if __name__ == "__main__":
	main()
