import os
import pandas as pd
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
EI      =     2*1e5             #   [N.m^2] Bending Rigidity
GJ      =     1*1e5             #   [N.m^2] Torsional Rigiditycasescases

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

def part1(dv:float=0.005, dt:int = 2, velocityRange:str='FlightEnvelope',case:str='part1Flutter',
          showPlots:bool=False):

    # Step 1: get velocity and altitude conditions:
    #   -- for each height, get a velocity vector.
    #   -- this velocity vector become the input.
    # will debug the problem initially assuming that the aircraft is full.

    # Conditions returned from the flight envelope for altitude and velocity.
    conditions = flEnv.conditions()
    arrALT =  np.unique(conditions[:,1])

    # Return flight conditions for the aircraft mission.
    IP_sweep = np.linspace(MMIF, MMIE, dt)
    chordCoeffCOM_sweep = np.linspace(CMF, CME, dt)
    mass_sweep = np.linspace(MWINGF, MWINGE, dt)
    time_sweep = np.linspace(0, ENDR, dt)

    print(chordCoeffCOM_sweep)
    
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
                    #print(f"\nMode {mode_key}:")
                    #print(f"  Flutter Velocity (V_flutter): {V_flutter}")
                    #print(f"  Flutter Frequency (freq_flutter): {freq_flutter}")
                    if showPlots:
                        plotDimensionlessFrequency(parameters[3], roots["R1"], roots["R2"], roots["R3"], roots["R4"], case)
                        plotDimensionlessDamping(parameters[3], roots["R1"], roots["R2"], roots["R3"], roots["R4"], case)
                else:
                    print(f"\nMode {mode_key}: No flutter detected in the given velocity range.")
            print("-"*65)

            input_vars_flutter = {
                'time[s]': time,
                'a': sectionModel[0],
                'e': sectionModel[1] + sectionModel[0],
                'mu': parameters[2],
                'rs': parameters[0],
                'sigma': parameters[1],
                'xTheta': sectionModel[1],
                'MMI': MMI,
                'MWING': MWING,
                'chordCoeffCOM': chordCoeffCOM,
                'h': h
            }
            
            writeResults(input_vars_flutter, flutter_points, filename=outfilename)

    columns_to_check = ['R1_V_flutter', 'R2_V_flutter', 'R3_V_flutter', 'R4_V_flutter']
    df = pd.read_csv(outfilename)
    
    if (df[columns_to_check] == 99999).all().all():
        print("")
        print_green("The design is flutter free for the current mission.")
    else:
        print("")
        print_red("The design has a flutter condition for this configuration")
        df_temp = pd.read_csv(outfilename)

        # Drop rows where all columns in columns_to_check have the value 99999
        df_temp = df_temp[~(df_temp[columns_to_check] == 99999).all(axis=1).to_numpy()]
        
        min_time = df_temp['time[s]'].min()
        max_time = df_temp['time[s]'].max()

        df_min_time = df_temp[df_temp['time[s]'] == min_time]
        df_max_time = df_temp[df_temp['time[s]'] == max_time]

        print(f"Data for the smallest unique value of time[s] ({min_time}):")
        print(df_min_time)

        print(f"Data for the largest unique value of time[s] ({max_time}):")
        print(df_max_time)

    print("")
    print_green(f"Results written to {outfilename}")
    print("")

def part1MDOF(dv:float=0.005, dt:int = 2, velocityRange:str='FlightEnvelope',case:str='part1Flutter',
          showPlots:bool=False):

    conditions = flEnv.conditions()
    arrALT =  np.unique(conditions[:,1])

    IP_sweep = np.linspace(MMIF, MMIE, dt)
    chordCoeffCOM_sweep = np.linspace(CMF, CME, dt)
    mass_sweep = np.linspace(MWINGF, MWINGE, dt)
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

            # Compute dimensionless parameters
            parameters = computeDimensionlessParameters(U_vec, MWING, h, EI, GJ, MMI, C/2, B/2, dv=dv, velocityRange=velocityRange)
            sectionModel = computeSectionModel(chordCoeffCOM, 0.4)

            # Now run the 2DOF method:
            # twoDOFMethod returns results, roots similar to pkmethod but with only R1 and R2.
            flutter_results, roots = multiDOFMethod(*parameters, *sectionModel, initK=1.0, printAlt=h)

            # Since we only have R1 and R2:
            flutter_points = findFlutter(parameters[3], roots) 
            # NOTE: Ensure findFlutter can handle modes=["R1","R2"] only.
            # If not, adapt it: for each mode in that list, check flutter.

            print("-"*65)
            for mode_key, flutter_info in flutter_points.items():
                if flutter_info is not None:
                    V_flutter = flutter_info['V_flutter']
                    freq_flutter = flutter_info['freq_flutter']
                    if showPlots:
                        plotDimensionlessFrequency(parameters[3], roots["R1"], roots["R2"], None, None, case) 
                        plotDimensionlessDamping(parameters[3], roots["R1"], roots["R2"], None, None, case)
                else:
                    print(f"\nMode {mode_key}: No flutter detected in the given velocity range.")
            print("-"*65)

            input_vars_flutter = {
                'time[s]': time,
                'a': sectionModel[0],
                'e': sectionModel[1] + sectionModel[0],
                'mu': parameters[2],
                'rs': parameters[0],
                'sigma': parameters[1],
                'xTheta': sectionModel[1],
                'MMI': MMI,
                'MWING': MWING,
                'chordCoeffCOM': chordCoeffCOM,
                'h': h
            }

            writeResults(input_vars_flutter, flutter_points, filename=outfilename)

    columns_to_check = ['R1_V_flutter', 'R2_V_flutter']
    df = pd.read_csv(outfilename)
    
    if (df[columns_to_check] == 99999).all().all():
        print("")
        print_green("The design is flutter free for the current mission.")
    else:
        print("")
        print_red("The design has a flutter condition for this configuration")
        df_temp = pd.read_csv(outfilename)

        # Filter out rows where all flutter velocities are 99999
        df_temp = df_temp[~(df_temp[columns_to_check] == 99999).all(axis=1).to_numpy()]
        
        min_time = df_temp['time[s]'].min()
        max_time = df_temp['time[s]'].max()

        df_min_time = df_temp[df_temp['time[s]'] == min_time]
        df_max_time = df_temp[df_temp['time[s]'] == max_time]

        print(f"Data for the smallest unique value of time[s] ({min_time}):")
        print(df_min_time)
        print(f"Data for the largest unique value of time[s] ({max_time}):")
        print(df_max_time)

    print("")
    print_green(f"Results written to {outfilename}")
    print("")


def part2(dv:float=0.005, dt:int = 1, velocityRange:str='FlightEnvelope',case:str='part1Flutter',
          showPlots:bool=False):


    # This part of the program is to perfrm flutter analysis on the all electric variant of the aircraft.
    # Conditions returned from the flight envelope for altitude and velocity.
    
    MWING = MWINGF
    MMI = 5.2
    chordCoeffCOM = 0.41

    conditions = flEnv.conditions()
    arrALT =  np.unique(conditions[:,1])

    
    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    outfilename = f"results/{case}_{current_time}.csv"
    velocity_dict = {}

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
                #print(f"\nMode {mode_key}:")
                #print(f"  Flutter Velocity (V_flutter): {V_flutter}")
                #print(f"  Flutter Frequency (freq_flutter): {freq_flutter}")
                if showPlots:
                    plotDimensionlessFrequency(parameters[3], roots["R1"], roots["R2"], roots["R3"], roots["R4"], case)
                    plotDimensionlessDamping(parameters[3], roots["R1"], roots["R2"], roots["R3"], roots["R4"], case)
            else:
                print(f"\nMode {mode_key}: No flutter detected in the given velocity range.")
        print("-"*65)

        input_vars_flutter = {
            'a': sectionModel[0],
            'e': sectionModel[1] + sectionModel[0],
            'mu': parameters[2],
            'rs': parameters[0],
            'sigma': parameters[1],
            'xTheta': sectionModel[1],
            'MMI': MMI,
            'MWING': MWING,
            'chordCoeffCOM': chordCoeffCOM,
            'h': h
        }
        
        writeResults(input_vars_flutter, flutter_points, filename=outfilename)

    columns_to_check = ['R1_V_flutter', 'R2_V_flutter', 'R3_V_flutter', 'R4_V_flutter']
    df = pd.read_csv(outfilename)
    
    if (df[columns_to_check] == 99999).all().all():
        print("")
        print_green("The design is flutter free for the current mission.")
    else:
        print("")
        print_red("The design has a flutter condition for this configuration")
        df_temp = pd.read_csv(outfilename)

        # Drop rows where all columns in columns_to_check have the value 99999
        df_temp = df_temp[~(df_temp[columns_to_check] == 99999).all(axis=1).to_numpy()]
        
        min_time = df_temp['time[s]'].min()
        max_time = df_temp['time[s]'].max()

        df_min_time = df_temp[df_temp['time[s]'] == min_time]
        df_max_time = df_temp[df_temp['time[s]'] == max_time]

        print(f"Data for the smallest unique value of time[s] ({min_time}):")
        print(df_min_time)

        print(f"Data for the largest unique value of time[s] ({max_time}):")
        print(df_max_time)

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
    flutter_results, roots = multiDOFMethod(*parameters, a, xTheta, initK, printAlt=0.0)

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
