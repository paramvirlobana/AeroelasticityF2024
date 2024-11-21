import matplotlib.pyplot as plt

def genFlightEnvelope(showPlot=False):
    x_points = [0, 252, 302, 302, 112.5, 0, 0]
    y_points = [0, 0, 0.8, 3, 3, 1.2, 0]

    # Create the figure and axis
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot the envelope outline
    ax.plot(x_points, y_points, 'k-', linewidth=2)

    plt.tight_layout()
    if showPlot:
        plt.show()

if __name__ == "__main__":
    genFlightEnvelope(showPlot=True)
