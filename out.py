import re
from matplotlib import pyplot as plt
import numpy as np



def process_log(log_file_path):
    with open(log_file_path, 'r') as file:
        lines = file.readlines()

    start_processing = False
    thread_data = []

    for line in lines:
        # Check for the line containing "Threads" to start processing
        if "Threads" in line:
            start_processing = True

        if start_processing:
            # Extract numbers after "Threads" and the number before and after +- in the specified format
            match = re.search(r'Threads\s(\d+)', line)
            if match:
                thread_data.append([int(match.group(1))])
            match = re.search(r'(\d+\.\d+)\s\+\-\s(\d+\.\d+)', line)
            if match:
                thread_data[-1].append(float(match.group(1)))
                thread_data[-1].append(float(match.group(2)))

    return thread_data

if __name__ == "__main__":
    log_file_path = "out.txt"
    result = process_log(log_file_path)
    x = np.array([i[0] for i in result])
    y = np.array([i[1] for i in result])
    err = np.array([i[2] for i in result])

    if result:
        plt.plot(x, y, 'k-')
        plt.fill_between(x, y-err,y+err)
        plt.show()
    else:
        print("No relevant data found in the log file.")


