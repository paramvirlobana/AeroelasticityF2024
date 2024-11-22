import numpy as np
from modules.methods import *

# Define all static variables:
# PROGRAM      VALUE                UNITS   DESCRIPTION
# VARIABLE
# -----------------------------------------------------------------
#                           AIRCRAFT
MTOW    =     997.9032          #   [kg]    Maximum Take-off Weight
SHP     =     145               #   [hp]    Horsepower
ENDR    =     4*60*60           #   [ s]    Aircraft Endurance
ALTMAX  =     3000              #   [ m]    Service Ceiling
UMAXSL  =     70                #   [m/s]   Maximum Speed at Sea Level (SL)
UMAXSC  =     83.8889           #   [m/s]   Maximum Speed at Service Ceiling (SC)
    
#                             WING
AR      =     1.17              #   [NA]    Wing Aspect Ratio
S       =     15.0              #   [m^2]   Wing Area
WFUEL   =     80                #   [kg]    Weight of Fuel in Each Wing
MMIF    =     7                 #   [kg.m^2]Mass Moment of Inertia (MMI) at MTOW
MMIE    =     4                 #   [kg.m^2]Mass Moment of Inertia (MMI) No Fuel
EI      =     2e5               #   [N.m^2] Bending Rigidity
GJ      =     1e5               #   [N.m^2] Torsional Rigidity

# CALCULATED VALUES
MWING   =     26.91*S           #   [kg]    Wing Mass
C       =     np.sqrt(S*AR)/2   #   [ m]    Wing chord length
CMF     =     0.35 * C
CME     =     0.45 * C

def part1():
    computeDimensionlessParameters(UMAXSL, MWING, 1.25, EI, GJ, MMIF, C)

def validation():
    """
    NOTE:
        -- This function calls the validation subroutine.
        -- This section defines all the test case variables. These do not
            correspond to the aircraft analysis problem. (Local function variables)
    """
    a:      float   =   -0.25
    e:      float   =   -0.1
    mu:     int     =   20
    rs:     float   =   6/25
    sigma:  float   =   0.4
    xTheta: float   =   e - a
    initK:  float   =   1.0
    case:   str     =   'validation'

    input_vars = {
        'a': a,
        'e': e,
        'mu': mu,
        'rs': rs,
        'sigma': sigma,
        'xTheta': xTheta
    }
    
    V_vec = np.arange(0, 4 + 0.05, 0.005)
    V_vec = V_vec[1:]

    # PRINT INITIALIZATION STATEMENT
    print("[INFO]   Validation subroutine has been called.")

    flutter_results, roots = pkmethod(sigma,   V_vec,  mu,
                                               a,       xTheta, rs,
                                               initK)

    flutter_points = findFlutter(V_vec, roots)

    # Print flutter results
    print("-"*65)
    for mode_key, flutter_info in flutter_points.items():
        if flutter_info is not None:
            V_flutter = flutter_info['V_flutter']
            freq_flutter = flutter_info['freq_flutter']
            print(f"\nMode {mode_key}:")
            print(f"  Flutter Velocity (V_flutter): {V_flutter}")
            print(f"  Flutter Frequency (freq_flutter): {freq_flutter}")
        else:
            print(f"\nMode {mode_key}: No flutter detected in the given velocity range.")
    print("-"*65)

    writeResults(input_vars, flutter_points)

    plotDimensionlessFrequency(V_vec, roots["R1"], roots["R2"], roots["R3"], roots["R4"], case)
    plotDimensionlessDamping(V_vec, roots["R1"], roots["R2"], roots["R3"], roots["R4"], case)

