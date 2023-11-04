from typing import List, Dict, Tuple
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


def read_times_from_file(filename: str) -> Tuple[List[int], Dict[int, List[float]]]:
    times_per_size: Dict[int, List[float]] = {}

    with open(filename) as fp:
        strides: List[int] = list(map(int, fp.readline().rstrip().split(", ")[1:]))

        while line := fp.readline().rstrip().split(", "):
            if line == [""]:
                break
            vec_size, times = int(line[0]), list(map(float, line[1:]))

            times_per_size[vec_size] = times

    return strides, times_per_size


def plot_times(strides: List[int], times_per_size: Dict[int, List[float]]) -> None:
    for vec_size in times_per_size.keys():
        times: List[float] = times_per_size[vec_size]
        plt.loglog(strides[: len(times)], times, label=f"{vec_size} bytes", marker="x")

    ax = plt.gca()
    fig = plt.gcf()
    fig.set_size_inches(16, 6)
    ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
    plt.legend()
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    plt.xlabel("Stride in bytes")
    plt.xscale("log", base=2)
    plt.ylabel("Time per memory access (ns)")
    plt.tight_layout()
    plt.savefig("ex1_plot.png")


if __name__ == "__main__":
    strides, times_per_size = read_times_from_file("out_30**2_O1.txt")
    plot_times(strides, times_per_size)
