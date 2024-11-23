
#===========================================================
#Program written for  Aeroelasticity MECH 6481 - Fall 2024
#                        UTILITIES
# Author:
#   -- Paramvir Lobana --
#===========================================================

import platform
import csv
import matplotlib.pyplot as plt

def writeResults(input_vars, flutter_points, filename="results/validation.csv"):
    """
    Writes the flutter analysis results and input variables to a CSV file.
    """
    fieldnames = ['a', 'e', 'mu', 'rs', 'sigma', 'xTheta']
    modes = ['R1', 'R2', 'R3', 'R4']
    for mode in modes:
        fieldnames.extend([f'{mode}_V_flutter', f'{mode}_freq_flutter'])

    data = input_vars.copy()
    for mode in modes:
        flutter_info = flutter_points.get(mode)
        if flutter_info is not None:
            V_flutter = flutter_info['V_flutter']
            freq_flutter = flutter_info['freq_flutter']
        else:
            V_flutter = 99999
            freq_flutter = 99999
        data[f'{mode}_V_flutter'] = V_flutter
        data[f'{mode}_freq_flutter'] = freq_flutter

    # Write to CSV
    with open(filename, mode='w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        writer.writerow(data)
    print("")
    print_green(f"Results written to {filename}")
    print("")


def print_green(text):
    """
    Function to print a string in Green color
    Args:
        text (str): String to colorize green
    """
    attr = ['42', '1']
    if platform.system() == "Windows":
        print(text)
    else:
        print('\x1b[%sm%s\x1b[0m' % (';'.join(attr), text))


def text_green(text):
    """
    Function to print a string in Green color
    Args:
        text (str): String to colorize green
    """
    attr = ['42', '1']
    if platform.system() == "Windows":
        return text
    else:
        return '\x1b[%sm%s\x1b[0m' % (';'.join(attr), text)


def print_red(text):
    """
    Function to print a text in Red for error notifications
    Args:
        text (str): String to colorize in Red
    """
    attr = ['41', '1']

    if platform.system() == "Windows":
        print(text)
    else:
        print('\x1b[%sm%s\x1b[0m' % (';'.join(attr), text))


def text_red(text):
    """
    Function to print a text in Red for error notifications
    Args:
        text (str): String to colorize in Red
    """
    attr = ['41', '1']

    if platform.system() == "Windows":
        return text
    else:
        return '\x1b[%sm%s\x1b[0m' % (';'.join(attr), text)

def plotSettings():
    plt.rcParams.update({
    'font.family': 'serif',  
    'font.serif': ['Times New Roman'],  
    'font.size':       11,  
    'axes.titlesize':  11,  
    'axes.labelsize':  11,  
    'legend.fontsize': 11, 
    'xtick.labelsize': 11, 
    'ytick.labelsize': 11  
    })

def printHead():
    head="""
    .         .                   .                                                      
       .                                                                                     
               .   .                   .                                                 
    Flutter Analysis Tool | PK METHOD | P METHOD |
    ------                      .       
    PROGRAMMED FOR AEROELASTICITY FALL 2024                           .
           .                 .                               .                            
      
                         .                .                     . .                      
                    .    .         .                                                     
                                                                                        .
                  .                      .      .                                       .
                   .                                  .                                  
    """

    print(head)

