import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def conditions(showPlot:bool=False, saveResults:bool=False):

    x_points = [0, 252/3.6, 302/3.6, 302/3.6, 112.5/3.6, 0, 0]
    y_points = [0, 0, 800, 3000, 3000, 1200, 0]
    xy = (x_points, y_points)
    points, mask = flightEnvelopePoints(*xy)
    if showPlot:
        plotEnvelope(points, *xy)
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


def plotEnvelope(points, x_points, y_points):
    plt.figure(figsize=(8, 5))
    
    # Plot the boundary lines first
    plt.plot(x_points, y_points, 'r-', linewidth=2, alpha=0.7, label='Flight Envelope')
    plt.scatter(points[:,0], points[:,1], s=2, alpha=0.5, label='Analysis Conditions')
    
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Altitude (m)')
    plt.title('Flight Envelope Points')
    plt.grid(True)
    plt.legend()
    
    plt.xlim(-10/3.6, 320/3.6)
    plt.ylim(-100, 3200)
    
    plt.show()


if __name__ == "__main__":
    conditions()
