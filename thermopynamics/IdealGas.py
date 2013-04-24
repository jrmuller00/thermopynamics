""" Ideal gas functions library

The functions in this library will be used to calcuate thermodynamic properties 
of an ideal gas

    P v = R T

"""

import ThermoConst as thermc
import math
import scipy.optimize


__version__ = "0.1"


#
# Calculate the pressure from Temperature and specific volume
# P = R T / v

def pvt(thermconsts):
    Rgas = thermc.R_BAR/thermconsts['MW']
    return Rgas * thermconsts['Temperature'] / thermconsts['SpVol']

#
# Calculate the Temperature from pressure and spcific volume
# T = P v / R
def tpv(thermconsts):
    Rgas = thermc.R_BAR/thermconsts['MW']
    print (Rgas, thermconsts['Pressure'] ,thermconsts['SpVol'])
    return thermconsts['Pressure'] * thermconsts['SpVol'] / Rgas


#
# Calculate the specific volume from pressure and spcific volume
# v = R T / P
def vtp(thermconsts):
    Rgas = thermc.R_BAR/thermconsts['MW']
    T = thermconsts['Temperature']
    P = thermconsts['Pressure']
    return Rgas *  T / P


#
# Calculate dp/dv|t = -R T / v^2

def dptdv(thermconsts):
    Rgas = thermc.R_BAR/thermconsts['MW']
    v = thermconsts['SpVol']
    T = thermconsts['Temperature']
    return -(Rgas * T / (v * v))


#
# Calculate dp/dt|v = R / v

def dpvdt(thermconsts):
    Rgas = thermc.R_BAR/thermconsts['MW']
    v = thermconsts['SpVol']
    return (Rgas / v)

# 
# calculate the dv integrals, udv = 0, sdv = R ln (v / vcrit)
# note returns 2 values 

def dvint(thermconsts):
    Rgas = thermc.R_BAR/thermconsts['MW']
    v = thermconsts['SpVol']
    T = thermconsts['Temperature']
    vcrit = thermconsts['vCrit']
    print ('sp vol = ',v)
    return 0, Rgas * math.log(v / vcrit)

def vts(thermo):
    thermo['sActual'] = thermo['SpEntropy']
    v0 = 2 * thermo['vCrit']
    thermtuple = 'Ideal vts',thermo
    #print (thermtuple, type(thermtuple))
    v = scipy.optimize.newton(vfts, v0, None, thermtuple)
    return v

def vfts(v,str, thermo):
    thermo['SpVol'] = v
    #print ('Ideal in vfts', v)
    thermo['Pressure'] = pvt(thermo)
    (cx, cxint, cxtint) = thermc.Cubiccv(thermo)
    (udvint, sdvint) = thermo['dvint'](thermo)
    scalc = cxtint + sdvint + thermo['sRef'];
    return scalc - thermo['sActual']


def tph(thermo):
    thermo['hActual'] = thermo['SpEnthalpy']
    T0 = thermo['TCrit']
    thermtuple = 'Ideal tph', thermo
    T = scipy.optimize.newton(tfph, T0, None, thermtuple)
    return T

def tfph(T, str, thermo):
    thermo['Temperature'] = T
    v = vtp(thermo)
    (cx, cxint, cxtint) = thermc.Cubiccv(thermo)
    (udvint, sdvint) = thermo['dvint'](thermo)
    u = cxint +udvint + thermo['xRef']
    h = u + thermo['Pressure']*v
    return h - thermo['hActual']
    
def tps(thermo):
    thermo['sActual'] = thermo['SpEntropy']
    T0 = thermo['TCrit']
    thermtuple = 'Ideal tps', thermo
    T = scipy.optimize.newton(tfps, T0, None, thermtuple)
    return T

def tfps(T, str, thermo):
    thermo['Temperature'] = T
    v = vtp(thermo)
    (cx, cxint, cxtint) = thermc.Cubiccv(thermo)
    (udvint, sdvint) = thermo['dvint'](thermo)
    scalc = cxtint + sdvint + thermo['sRef'];
    return scalc - thermo['sActual']

