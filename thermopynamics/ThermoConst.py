""" 
Constants file

This file contains thermodynamic or units
constants used throughout the package

    
    1-16-2013
    sjpfu
"""

__version__ = "0.1"

import math

R_BAR = 8314.0
FARADAY = 96418.0

def Cubiccv(therm):
    """
    Function to calculate a cubic polynomial fit
    for cv data

    cv = c3 T**3 + c2 T**2 + c1 T + c0

    cvint = INT(cv dT)|{T, Tcrit}

    cvtint = INT(cv/T dT)|{T,Tcrit}

    The parameters ci are specified in the dictionary therm

    """
    try:
        c0 = therm['cx0']
        c1 = therm['cx1']
        c2 = therm['cx2']
        c3 = therm['cx3']
        T = therm['Temperature']
        Tcrit = therm['TCrit']
        cv = cvint = cvtint = 0
        cv = c3 * T**3 + c2 * T**2 + c1 * T + c0
        cvint = (c3/4)*(T**4 - Tcrit**4) + (c2/3)*(T**3 - Tcrit**3) + (c1/2)*(T*T - Tcrit*Tcrit) + (c0*(T - Tcrit))
        cvtint = (c3*(T**3 - Tcrit**3)/3.0)  + (c2*(T**2 - Tcrit**2)/2.0) + (c1*(T - Tcrit)) + (c0*math.log(T/Tcrit))
    except:
        cv = 0
        cvint = 0
        cvtint = 0
        print ('Exception in Cubiccv')
    return cv, cvint, cvtint

