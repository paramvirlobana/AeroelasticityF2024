import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from ambiance import Atmosphere

def conditions(showPlot:bool=False, saveResults:bool=False, showConditions:bool=False):

    x_points = [0, 252/3.6, 302/3.6, 302/3.6, 112.5/3.6, 0, 0]
    y_points = [0, 0, 800, 3000, 3000, 1200, 0]
    EAS      = [0, 70.0000005178368, 80.70061865507049, 72.27350059255646, 26.923075551862922, 0, 0]
    EAS15    = [i * 1.15 for i in EAS]

    xy = (EAS15, y_points)
    XY = (x_points, y_points)
    points, mask = flightEnvelopePoints(*xy)
    if showPlot:
        plotEnvelope(points, *XY, EAS, EAS15, showConditions)
    if saveResults:
        np.savetxt('results/flight_envelope_points.csv', points, delimiter=',', header='velocity,altitude', comments='')

    return points

def flightEnvelopePoints(x_points, y_points, v_resolution=0.05, h_resolution=10):
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


def plotEnvelope(points, x_points, y_points, EAS, EAS15, showConditions:bool=False):
    fig, ax = plt.subplots(figsize=(10, 5))

    df_full = pd.read_csv("modules/full.csv")
    df_empty = pd.read_csv("modules/empty.csv")

    h_results_full = df_full['h']
    V_results_full = df_full['VV']

    h_results_empty = df_empty['h']
    V_results_empty = df_empty['VV']
    
    # Plot the boundary lines first
    ax.plot(x_points, y_points, color='black', linewidth=2, alpha=0.7)


    # EAS AND EAS*1.15 plotting
    x_pts = [x_points[1], EAS[1], EAS[2], EAS[3], x_points[3], x_points[2], x_points[1]]
    y_pts = [y_points[1], y_points[1],y_points[2], y_points[3], y_points[3], y_points[2], y_points[1]]
    x_pts15 = [EAS[1], EAS15[1], EAS15[2], EAS15[3], EAS[3], EAS[2], EAS[1]]

    #ax.plot(x_pts, y_pts, color='black', linewidth=2, alpha=0.7)
    #ax.plot(x_pts15, y_pts, color='black', linewidth=2, alpha=0.7)

    if showConditions:  
        ax.scatter(points[:,0], points[:,1], s=2, alpha=0.5, label='Analysis Conditions')

    ax.plot(V_results_full, h_results_full, color='black', label=r'$U_F$ with Max Fuel')
    ax.plot(V_results_empty, h_results_empty, color='red', linestyle='--', label=r'$U_F$ with Min Fuel')
    
    ax.set_xlabel('True speed (m/s)')
    ax.set_ylabel('Altitude (m)')
    # ax.set_title('Flight ')
    ax.fill(x_points, y_points, color='white', alpha=0.5, label='Proposed Flight Envelope')
    #ax.fill(x_pts, y_pts, color='green', alpha=0.5, label='EAS')
    ax.fill(x_pts15, y_pts, color='lightgray', alpha=0.5, label='EAS with Safety Margin')

    # ax.grid(True)
    ax.set_xlim(-5, max(x_points)*1.75)
    ax.set_ylim(-100, max(y_points)*1.15)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.legend(loc='upper left', bbox_to_anchor=(0, 1), ncol=1)

    # ax.text(87, 30, "252 km/h", ha='right')
    # ax.text(85, 800, "302 km/h", ha='left')
    # ax.text(82, 800, "800 m", ha='right')
    # ax.text(14, 1200, "1200 m", ha='right')
    # ax.text(85, 3000, "3000 m", ha='left')
    #

    plt.savefig("figs/FlightEnvelopeDetailed.eps", format="eps")
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
