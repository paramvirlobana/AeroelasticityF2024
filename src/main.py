#===========================================================
#Program written for  Aeroelasticity MECH 6481 - Fall 2024
#                        PROJECT
# Authors:
#   -- Paramvir Lobana -- Jeremy Burg -- Charles Gauthier --
#===========================================================

import numpy as np
import time
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import modules.cases as case
from modules.utils import *

__version__ = 1.0
"""
 Define all static variables:
 PROGRAM      VALUE                UNITS   DESCRIPTION
 VARIABLE
 -----------------------------------------------------------------
                           AIRCRAFT
MTOW    =     997.9032          #   [kg]    Maximum Take-off Weight
SHP     =     145               #   [hp]    Horsepower
ENDR    =     4*60*60           #   [ s]    Aircraft Endurance
ALTMAX  =     3000              #   [ m]    Service Ceiling
UMAXSL  =     70                #   [m/s]   Maximum Speed at Sea Level (SL)
UMAXSC  =     83.8889           #   [m/s]   Maximum Speed at Service Ceiling (SC)
   
                             WING
AR      =     1.17              #   [NA]    Wing Aspect Ratio
S       =     15.0              #   [m^2]   Wing Area
WFUEL   =     80                #   [kg]    Weight of Fuel in Each Wing
MMIF    =     7                 #   [kg.m^2]Mass Moment of Inertia (MMI) at MTOW
MMIE    =     4                 #   [kg.m^2]Mass Moment of Inertia (MMI) No Fuel
EI      =     2e5               #   [N.m^2] Bending Rigidity
GJ      =     1e5               #   [N.m^2] Torsional Rigidity

 CALCULATED VALUES
MWING   =     26.91*S           #   [kg]    Wing Mass
C       =     np.sqrt(S*AR)/2   #   [ m]    Wing chord length
CMF     =     0.35 * C
CME     =     0.45 * C
"""
"""
TODO:
    -- 
    --
    --
"""

def main():
    startTime = time.time()

    # SET THE ENVIRONMENT FOR ANALYSIS
    parser = argparse.ArgumentParser(description="Flutter Analysis tool | P-K METHOD | P METHOD")
    parser.add_argument('-v', '--validation', action='store_true', help="If passed, the validation logic is called along with part 1.")
    parser.add_argument('-z', '--validationo', action='store_true', help="If passed, the validation standalone is called")
    parser.add_argument('-p', '--plot', action='store_true', help="Shows plot for the validation case.")
    arguments = parser.parse_args()

    v = arguments.validation
    z = arguments.validationo
    p = arguments.plot

    """
    -- CALL THE REQUESTED METHOD.
    -- NOTE: MULTIPLE METHODS CAN BE CALLED AT THE SAME TIME.
    -- PK method is ALWAYS called since it is quite reliable.
    """
    if z == False:
        case.FLUTTER()

    if v == True:
        case.validation(p)

    if z == True:
        case.validation(p)

    endTime = time.time()
    print("STATS:")
    print("-"*6)
    print(f"Program took {(endTime - startTime) *10**3:10.03f}ms to execute.")


if __name__ == "__main__":
    printHead()
    main()
