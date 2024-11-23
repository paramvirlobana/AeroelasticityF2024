import os
import numpy as np
from datetime import datetime
from modules.methods import *
import modules.flightEnvelope as flEnv
from modules.utils import *

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
AR      =     7.15              #   [NA]    Wing Aspect Ratio
S       =     15.0              #   [m^2]   Wing Area
WFUEL   =     80                #   [kg]    Weight of Fuel in Each Wing
MMIF    =     7                 #   [kg.m^2]Mass Moment of Inertia (MMI) at MTOW
MMIE    =     4                 #   [kg.m^2]Mass Moment of Inertia (MMI) No Fuel
EI      =     2e5               #   [N.m^2] Bending Rigidity
GJ      =     1e5               #   [N.m^2] Torsional Rigidity

# CALCULATED VALUES
MWINGE  =     (26.91*S)/2       #   [kg]    Wing Mass
MWINGF  =     MWINGE + 80       #   [kg]    Wing Mass
B       =     np.sqrt(S*AR)     #   [ m]    Wing span total length
C       =     S / B             #   [ m]    Wing chord
CMF     =     0.35
CME     =     0.45


"""
TODO:
    -- Define the flight conditions where we need to evaluate for flutter. Two approaches:
        --- Introduce case in a files and check for each case.
        --- Use linspace to generate points over the entire flight envelope and validate the flutter criteria.


    NOW HOW TO DO IT?
    -- For each flight condition, we obtain a velocity vector.
    -- For example, at lets say, h=1000m, we have the same ambient conditions, so we can check for a single velocity vector.
"""


def part1(dv:float=0.05, dt:int = 100, velocityRange:str='FlightEnvelope',case:str='part1Flutter',
          showPlots:bool=False):

    # Step 1: get velocity and altitude conditions:
    #   -- for each height, get a velocity vector.
    #   -- this velocity vector become the input.
    # will debug the problem initially assuming that the aircraft is full.

    # Conditions returned from the flight envelope for altitude and velocity.
    conditions = flEnv.conditions()
    arrALT =  np.unique(conditions[:,1])

    # Return flight conditions for the aircraft mission.
    IP_sweep = np.linspace(MMIE, MMIF, dt)
    chordCoeffCOM_sweep = np.linspace(CME, CMF, dt)
    mass_sweep = np.linspace(MWINGE, MWINGF, dt)
    time_sweep = np.linspace(0, ENDR, dt)
    
    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    outfilename = f"results/{case}_{current_time}.csv"
    velocity_dict = {}

    for i in range(dt):
        MMI = IP_sweep[i]
        MWING = mass_sweep[i]
        chordCoeffCOM = chordCoeffCOM_sweep[i]
        time = time_sweep[i]

        for h in arrALT:
            U_vec = conditions[conditions[:,1] == h][:,0]
            velocity_dict[h] = U_vec

            # Create tuples for the aircraft parameters
            amb         = (U_vec, MWING, h)
            structural  = (EI, GJ, MMI)
            geometric   = (C/2, B/2)

            parameters = computeDimensionlessParameters(*amb, *structural, *geometric, dv=dv, velocityRange=velocityRange)
            sectionModel = computeSectionModel(chordCoeffCOM, 0.4)

            flutter_results, roots = pkmethod(*parameters, *sectionModel, initK=1.0, printAlt=h)
            flutter_points = findFlutter(parameters[3], roots)

            # Print flutter results
            print("-"*65)
            for mode_key, flutter_info in flutter_points.items():
                if flutter_info is not None:
                    V_flutter = flutter_info['V_flutter']
                    freq_flutter = flutter_info['freq_flutter']
                    print(f"\nMode {mode_key}:")
                    print(f"  Flutter Velocity (V_flutter): {V_flutter}")
                    print(f"  Flutter Frequency (freq_flutter): {freq_flutter}")
                    if showPlots:
                        plotDimensionlessFrequency(parameters[3], roots["R1"], roots["R2"], roots["R3"], roots["R4"], case)
                        plotDimensionlessDamping(parameters[3], roots["R1"], roots["R2"], roots["R3"], roots["R4"], case)
                else:
                    print(f"\nMode {mode_key}: No flutter detected in the given velocity range.")
            print("-"*65)

            input_vars_flutter = {
                'time[s]': time,
                'a': parameters[0],
                'e': 0.45,
                'mu': parameters[2],
                'rs': parameters[0],
                'sigma': parameters[1],
                'xTheta': parameters[1],
                'MMI': MMI,
                'MWING': MWING,
                'chordCoeffCOM': chordCoeffCOM,
                'h': h
            }
            
            writeResults(input_vars_flutter, flutter_points, filename=outfilename)
    print("")
    print_green(f"Results written to {outfilename}")
    print("")

def validation(showPlot:bool):
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
    parameters = (rs, sigma, mu, V_vec)
    flutter_results, roots = pkmethod(*parameters, a, xTheta, initK, printAlt=0.0)

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

    writeResultsValidation(input_vars, flutter_points)
    
    if showPlot:
        plotDimensionlessFrequency(V_vec, roots["R1"], roots["R2"], roots["R3"], roots["R4"], case)
        plotDimensionlessDamping(V_vec, roots["R1"], roots["R2"], roots["R3"], roots["R4"], case)
