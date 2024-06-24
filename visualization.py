from main import *
import matplotlib.pyplot as plt

def plot_peaks(peaklist):
    """
    Plot the peaks in a peaklist.
    """
    fig, ax = plt.subplots()
    ax.scatter(peaklist['Position F1'], peaklist['Position F2'])
    ax.set_xlabel('Position F1')
    ax.set_ylabel('Position F2')
    ax.set_title(f'Peaks in {peaklist.name}')
    plt.show()


if __name__ == "__main__":
    NCACB = PeakList('NCACB', ('13C', '15N'))
    plot_peaks(NCACB)
