from cProfile import label
from curses.ascii import isdigit
from matplotlib import pyplot as plt
import matplotlib
import numpy as np
from sys import argv

def main():
    if len(argv) < 2:
        print("Usage: plot.py <log_file>")
        exit(1)

    with open(argv[1], "r") as f:
        lines = f.readlines()
        i = 0
        fig_i = 0

        matplotlib.use("QT5Agg")
    
        plt.ion()
        plt.figure(argv[1])
        while i < len(lines):
            if lines[i].startswith("/*"):
                while i < len(lines) and not lines[i].startswith("*/"):
                    i += 1
            if lines[i].startswith("Height"):
                data: np.ndarray = np.empty(0)
                title = lines[i-1]
                metrics = lines[i].split()
                i += 1
                while i < len(lines) and isdigit(lines[i][0]):
                    data = np.append(data, np.asarray(lines[i].split()).astype(float))
                    i += 1
                data = data.reshape(-1, len(metrics)).T
                for j in range(6, len(metrics)-1):
                    plt.plot(data[0, :], data[j, :],
                             label=f"{title}", marker="o")

            else:
                i += 1
        plt.legend()
        plt.xlabel("Tree Height")
        plt.ylabel("Time (ms)")
        plt.grid()
        plt.ioff()
        plt.show()


if __name__ == "__main__":
    main()
