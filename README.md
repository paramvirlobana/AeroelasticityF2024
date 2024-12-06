# Flutter Analysis Tool 
Authors: ```Paramvir Lobana```, ```Jeremy Burg```, ```Charles Gauthier``` \
Program written for the course: ```Aeroelasticity``` ```Fall 2024``` ```Concordia University, Montreal```

## Description
This program uses PK method to evaluate the flutter condition in a given flight envelope. 

## Usage
The main entry point of the program is the ```main.py``` file. This file contains the arguments to run different cases within the program. This files does not contain any logical part and all the computations are carried out within the methods present in the ```modules``` directory. See **Directory Structure** section for more information. 

The ```case.py``` file contains the cases that a user can define and run to solve for a given flutter problem. The cases are defined within functions and then invoked from the ```main.py``` file. In the program, there are 3 predefined cases:
1. Validation case
2. Part 1 - Generic Aircraft
3. Part 2 - All Electric Aircraft - TO BE COMPLETED

### Example Usage
The program can be run directy using the main file. User can either run a single case using the predefind commands, or run mutiple cases by adding cases in the ```main()``` function. The following arguments can be passed for the predefined cases.

1. Running the validation case: \
The validation case is define to validate the PK method from the refernce it has been adopted from (Hodges and Pierce, 2011). The command below will run only the validation case. It can be run as follows:
```bash
python3 main.py -z
```


2. Running user defined cases: \
To run the user defined cases, given the user has invoked the case function in the ```main()``` function, no argument is required and can be invoked as follows:
```bash
python3 main.py
```

## Sample Case Setup
The following section shows how a user can define a sample case defined in ```cases.py```:
```python
def customcase1(dv:float=0.05):

    a:      float   =   -0.25
    e:      float   =   -0.1
    mu:     int     =   20
    rs:     float   =   6/25
    sigma:  float   =   0.4
    xTheta: float   =   e - a
    initK:  float   =   1.0
    case:   str     =   'validation'

    # Define reduced velocity range
    V_vec = np.arange(0, 4 + dv, dv)
    V_vec = V_vec[1:]

    parameters = (rs, sigma, mu, V_vec)
    flutter_results, roots = pkmethod(*parameters, a, xTheta, initK, printAlt=0.0)

    flutter_points = findFlutter(V_vec, roots)
```
Another more complex example is below. In this case, all the dimensionless parameters are calculated for each flight condition instead of defining for a single flight condition by the user.
```python
def customcase2(dv:float=0.05, dt:int = 100, velocityRange:str='FlightEnvelope',case:str='part1Flutter', showPlots:bool=False):

    conditions = flEnv.conditions()
    arrALT =  np.unique(conditions[:,1])

    IP_sweep = np.linspace(MMIE, MMIF, dt)
    chordCoeffCOM_sweep = np.linspace(CME, CMF, dt)
    mass_sweep = np.linspace(MWINGE, MWINGF, dt)
    time_sweep = np.linspace(0, ENDR, dt)

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

            parameters = computeDimensionlessParameters(*amb, *structural, *geometric, dv, velocityRange)
            sectionModel = computeSectionModel(chordCoeffCOM, 0.4)

            flutter_results, roots = pkmethod(*parameters, *sectionModel, initK=1.0, printAlt=h)
            flutter_points = findFlutter(parameters[3], roots)
```

The above case then must be called to the ```main.py``` file as follows:

```python
from modules.cases import *

def main():
    customcase1()
    customcase2()

if __name__ == "__main__":
    main()
```

## Directory Structure
```bash
AeroelasticityF2024/src
.
├── figs
├── main.py
├── modules
│   ├── __init__.py
│   ├── cases.py
│   ├── flightEnvelope.py
│   ├── methods.py
│   ├── utils.py
│   └── validation.csv
└── results
```


