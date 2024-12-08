import numpy as np
from ambiance import Atmosphere

# Initialize x_points (True Air Speed in m/s) and y_points (altitude in m)
x_points = [0, 252/3.6, 302/3.6, 302/3.6, 112.5/3.6, 0, 0]
y_points = [0, 0, 800, 3000, 3000, 1200, 0]

# Calculate density at each altitude
density = np.array([Atmosphere(y).density[0] for y in y_points])

# Sea level standard air density
rho_0 = 1.225  # kg/m³

# Calculate Equivalent Air Speed (EAS)
# EAS = TAS * sqrt(rho_0 / rho)
eas = np.array(x_points) * np.sqrt(rho_0 / density)

# Calculate EAS15 by multiplying EAS by 1.15
eas15 = 1.15 * eas

# Print results
print("True Air Speed (m/s):", x_points)
print("Altitude (m):", y_points)
print("Air Density (kg/m³):", density)
print("Equivalent Air Speed (m/s):", eas)
print("EAS15 (m/s):", eas15)
