import numpy as np
import cmath
import matplotlib.pyplot as plt
from modules.utils import *

def pkmethod(sigma:float,   V_vec:float,    mu:float,
             a:float,       xTheta:float,   rs:float,
             initK:float,    tol:float=0.001):
    """
     This function holds the main logic for the pk method.   
    """
    print("STARTING PK METHOD")
    print("-"*65)

    results = []
    roots = {"R1": [], "R2": [], "R3": [], "R4": []}

    for V in V_vec:

        k = initK
        eigenvaluesList = []

        for mode in range(4):
            k_mode = k
            while True:
                eigenvalues = pkMethodLogic(sigma, V, mu,
                                            k_mode, a, xTheta,
                                            rs)

                eigenvaluesSorted = sorted(eigenvalues, key=lambda p: p.imag)
                imaginaryP = eigenvaluesSorted[mode]
                kUpdated = imaginaryP.imag

                if np.isclose(kUpdated, k_mode, atol=tol):
                    break

                k_mode = kUpdated

            eigenvaluesList.append(imaginaryP)

        # The roots can be stored now:
        roots["R1"].append(eigenvaluesList[0])
        roots["R2"].append(eigenvaluesList[1])
        roots["R3"].append(eigenvaluesList[2])
        roots["R4"].append(eigenvaluesList[3])
        
        results.append((V, eigenvaluesList))

    # Convert roots to numpy arrays for plotting
    for key in roots:
        roots[key] = np.array(roots[key])

    return results, roots

def computeDimensionlessParameters(U, m, rho_inf, EI, GJ, I_P, b):
    l = 5
    omega_h     = (1.8751**2) * np.sqrt(EI / (m * l**3))
    omega_theta = (np.pi / 2) * np.sqrt(GJ / I_P * l)

    sigma = omega_h / omega_theta
    
    rs = I_P / (m * b**2)

    mu = m / (rho_inf * np.pi * b**2)

    V = U / (b * omega_theta)
    print(rs, sigma, mu, V)

    return rs, sigma, mu, V



def TheodorsenFunction(k:float):
    numerator = complex(0.01365, 0.2808 * k) - (k**2 / 2)
    denominator = complex(0.01365, 0.3455 * k) - k**2
    
    ck = numerator / denominator

    return ck

def pkMethodLogic(sigma:float,   V:float,    mu:float,
                 k:float,       a:float,    xTheta:float,
                 rs:float):
    """
    This method defines the logic to be used in pk method from the p method perspective.
    """
    # Define components of the characteristic equation.

    i  = complex(0, 1)
    Ck = TheodorsenFunction(k)

    F11 = (sigma**2 / V**2) - (k**2 / mu) + ((2 * i * k * Ck) / (mu))
    F12 = (k * (i + a*k) + (2 + i*k*(1 - 2*a)) * Ck) / (mu)
    F21 = (a * k**2 - ((i * k * (1 + 2*a)) * Ck)) / (mu)
    F22 = ((8 * mu * rs / V**2) + (4 * i * (1 + 2*a) * (2*i - k*(1 - 2*a)) * Ck) - k * (k - 4*i + 8*a*(i + a*k))) / (8*mu)

    p = [rs - xTheta**2,
        0,
        (F22 + F11 * rs - F21 * xTheta - F12 * xTheta),
        0,
        ((F11 * F22) - (F12 * F21))]

    roots = np.roots(p)
    return roots

def findFlutter(V_vec, roots):
    """
    Analyzes the damping data to find the flutter velocity and frequency.
    """
    flutter_results = {}
    for mode_key in roots.keys():

        damping = np.real(roots[mode_key])              # Modal damping     -> Real Part
        frequency = V_vec * np.imag(roots[mode_key])    # Modal frequency   -> Imaginary Part

        # Find the index where the damping terms cross the x axis.
        # np.sign -> This did not work with p-method but works well with pk-method.
        sign_changes = np.where(np.diff(np.sign(damping)) > 0)[0]

        if sign_changes.size > 0:
            idx = sign_changes[0]                           # The first crossing point is considered as the flutter point.
            V1, V2 = V_vec[idx], V_vec[idx + 1]             # Linear interpolation to get the exact point.
            d1, d2 = damping[idx], damping[idx + 1]         # |
            V_flutter = V1 - d1 * (V2 - V1) / (d2 - d1)     # |

            freq1, freq2 = frequency[idx], frequency[idx + 1]                       # Get frequency using the flutter point index
            freq_flutter = freq1 + (freq2 - freq1) * (V_flutter - V1) / (V2 - V1)   # Linear interpolation 

            """
            NOTE:
                -- When using the previous logic, the program passed the mode
                    with negative frequency as an instable mode.
                -- The following logic checks if negative frequency is obtained.
                -- If yes, the mode is ignored.
            """
            if freq_flutter < 0:
                print_red(f"Negative frequency at flutter point for {mode_key}, ignoring.")
                flutter_results[mode_key] = None        # Set mode key to NONE.
            else:
                flutter_results[mode_key] = {
                    'V_flutter': V_flutter,
                    'freq_flutter': freq_flutter
                }
                print(f"Flutter detected for {mode_key}:")
                print(f"  Flutter Velocity (V_flutter): {V_flutter}")
                print(f"  Flutter Frequency (freq_flutter): {freq_flutter}")
        else:
            print(f"No flutter detected for {mode_key} in the given velocity range.")
            flutter_results[mode_key] = None

    return flutter_results

def kMethodLogic():
    print_green("STARTING K METHOD")

def plotDimensionlessFrequency(V_vec, r1, r2, r3, r4, case):
    plt.figure()
    plt.plot(V_vec, V_vec * np.imag(r1), 'ro', markersize=5, label='Root 1')
    plt.plot(V_vec, V_vec * np.imag(r2), 'bd', markersize=5, label='Root 2')
    plt.plot(V_vec, V_vec * np.imag(r3), 'k*', markersize=5, label='Root 3')
    plt.plot(V_vec, V_vec * np.imag(r4), 'gs', markersize=5, label='Root 4')
    plt.xlabel('V', fontsize=12)
    plt.ylabel(r'$\Omega/\omega_\theta$', fontsize=12)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    #plt.savefig(f"figs/{case}_Frequency.eps", format="eps")

def plotDimensionlessDamping(V_vec, r1, r2, r3, r4, case):
    plt.figure()
    plt.plot(V_vec, V_vec * np.real(r1), 'ro', markersize=5, label='Root 1')
    plt.plot(V_vec, V_vec * np.real(r2), 'bd', markersize=5, label='Root 2')
    plt.plot(V_vec, V_vec * np.real(r3), 'k*', markersize=5, label='Root 3')
    plt.plot(V_vec, V_vec * np.real(r4), 'gs', markersize=5, label='Root 4')
    plt.xlabel('V', fontsize=12)
    plt.ylabel(r'$\Gamma/\omega_\theta$', fontsize=12)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    #plt.savefig(f"figs/{case}_Damping.eps", format="eps")
    plt.show()
