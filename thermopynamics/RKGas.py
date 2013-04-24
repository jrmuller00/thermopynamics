""" Redlich-Kwong gas functions library

The functions in this library will be used to calcuate thermodynamic properties 
of an ideal gas

     P(T,v) = RT/v-b  - a/(v(v+b)T^(1/2))

"""

import ThermoConst as thermc
import math
import scipy.optimize
import IdealGas as ideal


__version__ = "0.1"


#
# Calculate the pressure from Temperature and specific volume
# P(T,v) = RT/v-b  - a/(v(v+b)T^(1/2))

def pvt(thermconsts):
    Rgas = thermc.R_BAR/thermconsts['MW']
    T = thermconsts['Temperature']
    v = thermconsts['SpVol']
    a = thermconsts['RKa']
    b = thermconsts['RKb']
    return (Rgas*T/(v-b)) - (a/(v*(v+b)*math.sqrt(T))) # Return Pressure = [R T/(v-b)] - a/v(v+b)t^1/2


def vtp(thermo):
    thermo['pActual'] = thermo['Pressure']
    v0 = ideal.vtp(thermo)
    thermtuple = 'RK vtp',thermo
    #print (thermtuple, type(thermtuple))
    v = scipy.optimize.newton(vftp, v0, None, thermtuple)
    return v


def vftp(v,str, thermo):
    thermo['SpVol'] = v
    #print ('RK in vftp', v)
    pcalc = pvt(thermo)
    return pcalc - thermo['pActual']


def vts(thermo):
    thermo['sActual'] = thermo['SpEntropy']
    v0 = 2 * thermo['vCrit']
    thermtuple = 'RK vts',thermo
    #print (thermtuple, type(thermtuple))
    v = scipy.optimize.newton(vfts, v0, None, thermtuple)
    return v

def vfts(v,str, thermo):
    thermo['SpVol'] = v
    #print ('RK in vftp', v)
    thermo['Pressure'] = pvt(thermo)
    (cx, cxint, cxtint) = thermc.Cubiccv(thermo)
    (udvint, sdvint) = thermo['dvint'](thermo)
    scalc = cxtint + sdvint + thermo['sRef'];
    return scalc - thermo['sActual']

def tpv(thermo):
    thermo['pActual'] = thermo['Pressure']
    T0 = ideal.tpv(thermo)
    thermtuple = 'RK tpv',thermo
    #print (thermtuple, type(thermtuple))
    T = scipy.optimize.newton(tfpv, T0, None, thermtuple)
    return T


def tfpv(T,str, thermo):
    thermo['Temperature'] = T
    #print ('RK in vftp', v)
    pcalc = pvt(thermo)
    return pcalc - thermo['pActual']




# 
# calculate the dv integrals, udv = 0, sdv = R ln (v / vcrit)
# note returns 2 values 

def dvint(thermconsts):
    Rgas = thermc.R_BAR/thermconsts['MW']
    T = thermconsts['Temperature']
    P = thermconsts['Pressure']
    v = thermconsts['SpVol']
    vcrit = thermconsts['vCrit']
    a = thermconsts['RKa']
    b = thermconsts['RKb']
    #print ('sp vol = ',v)

    vpb = v + b
    vcpb = vcrit + b
    vmb = v - b
    vcmb = vcrit - b
    lnv = math.log(v*vcpb/(vcrit*vpb))
    sqt = math.sqrt(T)
    udvint = 3*a*lnv/(2*b*sqt);

    # Calculate Integral sdvint

    sdvint = (a*lnv/(2*b*sqt*T)) + (Rgas*math.log(vmb/vcmb));
    return udvint, sdvint


def tph(thermo):
    thermo['hActual'] = thermo['SpEnthalpy']
    T0 = 0.95*thermo['TCrit']
    thermtuple = 'RK tph', thermo
    T = scipy.optimize.newton(tfph, T0, None, thermtuple)
    return T

def tfph(T, str, thermo):
    thermo['Temperature'] = T
    # print ('T = ', T)
    v = vtp(thermo)
    (cx, cxint, cxtint) = thermc.Cubiccv(thermo)
    (udvint, sdvint) = thermo['dvint'](thermo)
    u = cxint +udvint + thermo['xRef']
    h = u + thermo['Pressure']*v
    return h - thermo['hActual']

def tps(thermo):
    thermo['sActual'] = thermo['SpEntropy']
    T0 = thermo['TCrit']
    thermtuple = 'RK tps', thermo
    T = scipy.optimize.newton(tfps, T0, None, thermtuple)
    return T

def tfps(T, str, thermo):
    thermo['Temperature'] = T
    v = vtp(thermo)
    (cx, cxint, cxtint) = thermc.Cubiccv(thermo)
    (udvint, sdvint) = thermo['dvint'](thermo)
    scalc = cxtint + sdvint + thermo['sRef'];
    return scalc - thermo['sActual']


