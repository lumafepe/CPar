import argparse
import re

from typing import AnyStr

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import LogFormatter


def process_log_perf(file_path: str, mode: str) -> list[list[float | int]]:
    def normalize_time(minute_str: str, seconds_str: str) -> float:
        return float(minute_str) * 60 + float(seconds_str)

    with open(file_path, 'r') as file:
        lines: list[AnyStr] = file.readlines()

    start_processing: bool = False
    thread_data: list[list[float | int]] = []

    line: AnyStr
    for line in lines:

        # Check for the line containing "Threads" to start processing
        if "Threads" in line:
            start_processing = True

        if start_processing:
            # Extract numbers after "Threads" and the number before and after +- in the specified format
            match = re.search(r'Threads\s(\d+)', line)
            if match:
                thread_data.append([int(match.group(1))])

            if mode == "perf":
                match = re.search(r'(\d+\.\d+)\s\+\-\s(\d+\.\d+)', line)
                if match:
                    thread_data[-1].append(float(match.group(1)))
                    thread_data[-1].append(float(match.group(2)))

            elif mode == "time":
                match = re.search(r'real\t(\d+)m(\d+\.\d+)s', line)  # real    0m0.775s

                if match:
                    normalized_time: float = normalize_time(match.group(1), match.group(2))
                    thread_data[-1].append(normalized_time)

    return thread_data


def plot_time(result: list[list[float | int]], plot_name: str = "plot.png"):
    data_array = np.array(result)

    x = data_array[:, 0]
    y = data_array[:, 1]

    outlier_threshold: int = 90

    # Replace outliers with the mean of the adjacent values
    for i in range(len(y)):
        if 0 < i < len(y) - 1:
            if y[i] > outlier_threshold:
                y[i] = (y[i - 1] + y[i + 1]) / 2.0

    # Create the plot
    plt.figure(figsize=(8, 6))
    plt.scatter(x, y, color='blue', label='Execution Time')
    plt.plot(x, y, linestyle='-', color='red')  # Connect points with lines

    for i, txt in enumerate(y):
        plt.annotate(f'{txt}s', (x[i], y[i]), textcoords="offset points", xytext=(0, 12), ha='center', rotation=90)

    plt.title('Number of Threads vs. Execution Time')
    plt.xlabel('Number of Threads')
    plt.ylabel('Execution Time')

    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.yscale('log')

    class LogFormatterIntegers(LogFormatter):
        def __call__(self, xx, pos=None):
            return "{:.0f}".format(xx)  # Display integers for the tick positions

    plt.gca().yaxis.set_major_formatter(LogFormatterIntegers())

    # Show the plot
    plt.savefig(plot_name)


def main(file_path: str, mode: str) -> None:
    result: list[list[float | int]] = process_log_perf(file_path, mode)

    if mode == "perf" and result:

        x = np.array([i[0] for i in result])
        y = np.array([i[1] for i in result])

        err = np.array(i[2] for i in result)

        plt.plot(x, y, 'k-')
        plt.fill_between(x, y - err, y + err)
        plt.show()

    elif mode == "time" and result:
        plot_time(result)

    else:
        print("No relevant data found in the log file. Exiting...")


if __name__ == "__main__":

    log_file_path = "out.txt"

    parser = argparse.ArgumentParser(description='mode in which the log_file_path is set to')
    parser.add_argument('argument', choices=['time', 'perf'], help='specify either "time" or "perf"')

    args = parser.parse_args()

    # Check if the argument provided is 'time' or 'perf'
    if args.argument == 'time':
        SystemExit(main(log_file_path, "time"))
    elif args.argument == 'perf':
        SystemExit(main(log_file_path, "perf"))
