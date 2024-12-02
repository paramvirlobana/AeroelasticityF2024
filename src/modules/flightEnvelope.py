import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def conditions(showPlot:bool=False, saveResults:bool=False, showConditions:bool=False):

    x_points = [0, 252/3.6, 302/3.6, 302/3.6, 112.5/3.6, 0, 0]
    y_points = [0, 0, 800, 3000, 3000, 1200, 0]
    xy = (x_points, y_points)
    points, mask = flightEnvelopePoints(*xy)
    if showPlot:
        plotEnvelope(points, *xy, showConditions)
    if saveResults:
        np.savetxt('results/flight_envelope_points.csv', points, delimiter=',', header='velocity,altitude', comments='')

    return points

def flightEnvelopePoints(x_points, y_points, v_resolution=2, h_resolution=500):
    """
    Generate points within the flight envelope with exact boundary points.
    
    Parameters:
    v_resolution: Resolution for velocity points (km/h)
    h_resolution: Resolution for altitude points (m)
    
    Returns:
    tuple: (valid_points, mask)
        - valid_points: List of [velocity, altitude] points within the envelope
        - mask: 2D boolean array indicating valid points
    """
    
    bottom_interp = interp1d([0, 252/3.6], [0, 0])
    
    right_diag_interp = interp1d([252/3.6, 302/3.6], [0, 800])
    
    top_interp = interp1d([112.5/3.6, 302/3.6], [3000, 3000])
    
    left_diag_interp = interp1d([0, 112.5/3.6], [1200, 3000])
    
    v_points = np.arange(0, 302/3.6 + v_resolution, v_resolution)
    h_points = np.arange(0, 3000 + h_resolution, h_resolution)
    V, H = np.meshgrid(v_points, h_points)
    
    mask = np.zeros_like(V, dtype=bool)
    
    for i, h in enumerate(h_points):
        for j, v in enumerate(v_points):
            if v <= 252/3.6:
                if v <= 112.5/3.6:
                    upper_bound = left_diag_interp(v)
                else:
                    upper_bound = 3000
                mask[i,j] = (h >= 0 and h <= upper_bound)
            
            elif v <= 302/3.6:
                lower_bound = right_diag_interp(v)
                mask[i,j] = (h >= lower_bound and h <= 3000)
    
    valid_points = np.column_stack((V[mask], H[mask]))
    
    return valid_points, mask


def plotEnvelope(points, x_points, y_points, showConditions:bool=False):
    fig, ax = plt.subplots(figsize=(8, 5))
    
    # Plot the boundary lines first
    ax.plot(x_points, y_points, color='black', linewidth=2, alpha=0.7, label='Flight Envelope')

    if showConditions:
        ax.scatter(points[:,0], points[:,1], s=2, alpha=0.5, label='Analysis Conditions')
    
    ax.set_xlabel('True speed (m/s)')
    ax.set_ylabel('Altitude (m)')
    # ax.set_title('Flight ')
    ax.fill(x_points, y_points, color='lightgray', alpha=0.5)

    # ax.grid(True)
    ax.set_xlim(-5, max(x_points)*1.15)
    ax.set_ylim(-100, max(y_points)*1.15)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


    ax.text(87, 30, "252 km/h", ha='right')
    ax.text(85, 800, "302 km/h", ha='left')
    ax.text(82, 800, "800 m", ha='right')
    ax.text(14, 1200, "1200 m", ha='right')
    ax.text(85, 3000, "3000 m", ha='left')


    plt.savefig("figs/FlightEnvelope.eps", format="eps")
    plt.show()

def plotSettings():
    plt.rcParams.update({
    'font.family': 'serif', 
    'font.serif': ['Times New Roman'],  
    'font.size':       16,  
    'axes.titlesize':  16,  
    'axes.labelsize':  16,  
    'legend.fontsize': 16, 
    'xtick.labelsize': 16, 
    'ytick.labelsize': 16  
    })

if __name__ == "__main__":
    plotSettings()
    conditions(showPlot=True)
