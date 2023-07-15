import csv
from typing import Tuple, List, Any


def get_lines_for_run(csv_file: str, run_number: int) -> tuple:
	lines = []
	run_info = []
	with open(csv_file, 'r') as file:
		reader = csv.DictReader(file)
		row: dict
		for row in reader:
			if int(row["run"]) == run_number:
				run_info.append(list(row.values())[:2])
				lines.append(list(map(float, list(row.values())[2:])))
			if int(row["run"]) > run_number:
				break
	return run_info, lines


def main() -> None:
	csv_file = 'data.csv'
	run_number = 2
	run_info, lines = get_lines_for_run(csv_file, run_number)
	# Print the lines for the specified run number
	for line in lines:
		print(line)


if __name__ == "__main__":
	main()

