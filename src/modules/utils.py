
#===========================================================
#Program written for  Aeroelasticity MECH 6481 - Fall 2024
#                        UTILITIES
# Author:
#   -- Paramvir Lobana --
#===========================================================

import platform
import matplotlib.pyplot as plt

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
    Flutter Analysis Tool | PK METHOD | P METHOD | K METHOD |
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

