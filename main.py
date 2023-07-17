import argparse
from argparse import Namespace
from rich.prompt import Prompt
import numpy as np

from eutectic import run_from_csv, run_eutectic


def parse_cli() -> Namespace:
	parser = argparse.ArgumentParser()

	parser.add_argument(
		"-m", "--manual",
		action="store_true",
		help="run the solver in manual.",
	)

	parser.add_argument(
		"-f", "--file",
		help="Location of the csv file to use."
	)

	parser.add_argument(
		"-r", "--run",
		type=int,
		help="Run number, if using a csv file with the -f flag."
	)

	parser.add_argument(
		"-e", "--enthalpy",
		nargs="*",
		type=float,
		help="Value of the enthalpy in J/mol."
	)

	parser.add_argument(
		"-t", "--temperature",
		nargs="*",
		type=float,
		help="Values of the Temperature."
	)

	parser.add_argument(
		"-c", "--celsius",
		action="store_true",
		help="Flag to set the temperature are given in Celsius (not in Kelvin)."
	)

	parser.add_argument(
		"-s", "--symbols",
		nargs="*",
		help="Specify the symbols of the element to be using"
	)

	return parser.parse_args()


def main() -> None:
	args = parse_cli()

	if args.file:
		run_from_csv(args.file, args.run, args.celsius)
	elif args.symbols:
		Hi = []
		Ti = []
		run_eutectic(Hi, Ti, use_celsius=False)
	else:
		if args.manual:
			temp_unit = "C" if args.celsius else "K"
			print("\nSupports multiple values, separate by spaces\n")
			Hi = list(map(float, Prompt.ask("Enthalpy (J/mol)").split()))
			Ti = list(map(float, Prompt.ask(f"Temperature ({temp_unit})").split()))
		else:
			Hi = args.enthalpy
			Ti = args.temperature
		if args.celsius:
			Ti = np.array(Ti) + 273.15
		run_eutectic(Hi, Ti, args.celsius)


if __name__ == "__main__":
	main()
