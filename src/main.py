#===========================================================
#Program written for  Aeroelasticity MECH 6481 - Fall 2024
#                        PROJECT
# Authors:
#   -- Paramvir Lobana -- Jeremy Burg -- Charles Gauthier --
#===========================================================

import time
import argparse
import modules.cases as case
from modules.utils import *

__version__ = 1.0

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
        case.part1(velocityRange='FlightEnvelope')

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
