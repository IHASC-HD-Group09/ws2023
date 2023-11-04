from typing import List, Dict, Tuple
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update(matplotlib.rcParamsDefault)


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
        plt.plot(strides[: len(times)], times, label=vec_size)

    plt.legend()
    plt.yscale("log")
    plt.xlabel("Stride")
    plt.ylabel("Time (s)")
    plt.savefig("ex1_plot.png")


if __name__ == "__main__":
    strides, times_per_size = read_times_from_file("out.txt")
    plot_times(strides, times_per_size)
