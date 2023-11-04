import matplotlib.pyplot as plt

def removeNewLine(line):
    line[-1] = line[-1].replace('\n', '')

def readTimesFromFile(filename):
    timesPerSize = {}
    with open(filename) as fp:
        strides = fp.readline().split(", ")[1:]
        removeNewLine(strides)
        while line := fp.readline().split(", "):
            if line == ['']:
                break
            removeNewLine(line)
            timesPerSize[line[0]] = line[1:]

    return strides, timesPerSize

def plotTimes(strides, timesPerSize):
    for vecSize in timesPerSize.keys():
        times = timesPerSize[vecSize]
        plt.loglog(strides[:len(times)], times, label=vecSize)

    plt.legend()
    plt.savefig('ex1_plot.png')



if __name__ == "__main__":
    strides, timesPerSize = readTimesFromFile("out.txt")
    plotTimes(strides, timesPerSize)

