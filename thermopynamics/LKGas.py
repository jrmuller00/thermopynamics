""" Lee-Kessler gas functions library

The functions in this library will be used to calcuate thermodynamic properties 
of an ideal gas

     Z = Pr vr'/ Tr = 1 + B/vr' + C/(vr')^2 + D/(vr')^2 + [c4/(Tr^3 vr'^2)](beta + gamma / vr'^2) exp(- gamma/ vr'^2)

     Pressure = Pc [Tr/vrp] Z = [T*R/v] Z

"""

import ThermoConst as thermc
import math
import scipy.optimize
import IdealGas as ideal


__version__ = "0.1"

#
# Function calcLKexp
def pvt(thermconsts):
    """ Function CalcLKexp will calculate the value for
  exponential term in the Lee Kessler EOS.  The function
  is general in that it will calculate either the value
  for the exponential term for the simple fluid or the
  reference fluid depending on the portion
  of the constants array that is passed.  IE, the function
  simply takes the first four doubles from the list and uses
  them for the calculation

      exp term =  [c4/(Tr^3 vr'^2)] [beta + (gamma/vr'^2)] exp(-gamma/vr'^2)

  INPUT:
      DOUBLE      *expvals          array containing the parameters
                                  b1 ... b4
      DOUBLE      Tr              reduced temperature
      DOUBLE      vrp             reduced specific volume as defined by
                                  Lee-Kessler EOS, ie vr' = v/[R Tc / Pc]

  OUTPUT:
      DOUBLE      exp term        value of the exponential term [returned]

    """
    Rgas = thermc.R_BAR/thermconsts['MW']
    T = thermconsts['Temperature']
    Tc = thermconsts['TCrit']
    Pc = thermconsts['PCrit']
    omega = thermconsts['Accentric']
    omegaf = v = thermconsts['LKomegaref']
    v = thermconsts['SpVol']
    thermconsts['LKvrp'] = v*Pc/(Rgas*Tc)
    thermconsts['LKTr'] = T/Tc
    thermconsts['LKFluid'] = 'simple'
    Z = calcLKZ(thermconsts)
    thermconsts['LKFluid'] = 'ref'
    Zf = calcLKZ(thermconsts)
    Z1 = (Zf - Z)/omegaf
    Z = Z + omega*Z1

    return T*Rgas*Z/v # Return Pressure = Pc [Tr/vrp] Z = [T*R/v] Z

def vtp(thermo):
    thermo['pActual'] = thermo['Pressure']
    v0 = ideal.vtp(thermo)
    thermtuple = 'LK vtp',thermo
    #print (thermtuple, type(thermtuple))
    v = scipy.optimize.newton(vftp, v0, None, thermtuple)
    return v

def vftp(v,str, thermo):
    thermo['SpVol'] = v
    #print ('RK in vftp', v)
    pcalc = pvt(thermo)
    return pcalc - thermo['pActual']


#
# Function calcLKZ
def calcLKZ(thermconsts):
    """ Function CalcLKZ will calculate the value for the 
      compressibility using the Lee Kessler EOS.  The function
      is general in that it will calculate either the value
      for the compressibility for the simple fluid or the
      reference fluid depending on the portion
      of the constants array that is passed.  IE, the function
      simply takes the first value of the array to be the
      beginning of the constants for the Lee Kessler EOS
      and uses them for the calculation


            Z = 1 + B/vr' + C/(vr')^2 + D/(vr')^2 + [c4/(Tr^3 vr'^2)](beta + gamma / vr'^2) exp(- gamma/ vr'^2)

      INPUT:
          DOUBLE      *zvals          array containing the parameters
                                      b1 ... b4
          DOUBLE      Tr              reduced temperature
          DOUBLE      vrp             reduced specific volume as defined by
                                      Lee-Kessler EOS, ie vr' = v/[R Tc / Pc]

      OUTPUT:
          DOUBLE      Z               compressibility [returned]
    """
    vrp = thermconsts['LKvrp']
    v2 = vrp**2
    v5 = vrp**5
    Tr = thermconsts['LKTr']
    Z = 1 + calcB(thermconsts)/vrp
    Z = Z + calcC(thermconsts)/v2
    Z = Z + calcD(thermconsts)/v5
    Z = Z + calcexp(thermconsts)
    return Z

# 
# Function calcLKB
def calcB(thermconsts):
    """  Function CalcLKB will calculate the value for B in the
  Lee Kessler EOS.  The function is general in that it
  will calculate either the value for B for the simple
  fluid or the reference fluid depending on the portion
  of the constants array that is passed.  IE, the function
  simply takes the first four doubles from the list and uses
  them for the calculation

      B(Tr) = b1 - b2/Tr - b3/Tr^2 - b4/Tr^3

  INPUT:
      DOUBLE      *bvals          array containing the parameters
                                  b1 ... b4
      DOUBLE      Tr              reduced temperature

   OUTPUT:
      DOUBLE      B               value of the function B [returned]


    """
    if thermconsts['LKFluid'] == 'ref':
        b1 = thermconsts['LKb1ref']
        b2 = thermconsts['LKb2ref']
        b3 = thermconsts['LKb3ref']
        b4 = thermconsts['LKb4ref']
    elif thermconsts['LKFluid'] == 'simple':
        b1 = thermconsts['LKb1simple']
        b2 = thermconsts['LKb2simple']
        b3 = thermconsts['LKb3simple']
        b4 = thermconsts['LKb4simple']
    else:
        print ('LK error in calcLKB')

    Tr = thermconsts['LKTr']

    return b1 - (b2/Tr) - (b3/Tr**2) - (b4/Tr**3)

#
# function calcLKC
def calcC(thermconsts):
    """  Function CalcLKC will calculate the value for C in the
  Lee Kessler EOS.  The function is general in that it
  will calculate either the value for C for the simple
  fluid or the reference fluid depending on the portion
  of the constants array that is passed.  IE, the function
  simply takes the first three doubles from the list and uses
  them for the calculation

      C(Tr) = c1 - c2/Tr + c3/Tr^2

  INPUT:
      DOUBLE      *cvals          array containing the parameters
                                  b1 ... b4
      DOUBLE      Tr              reduced temperature

  OUTPUT:
      DOUBLE      C               value of the function C[returned]
    """
    if thermconsts['LKFluid'] == 'ref':
        c1 = thermconsts['LKc1ref']
        c2 = thermconsts['LKc2ref']
        c3 = thermconsts['LKc3ref']    
    elif thermconsts['LKFluid'] == 'simple':
        c1 = thermconsts['LKc1simple']
        c2 = thermconsts['LKc2simple']
        c3 = thermconsts['LKc3simple']
    else:
        print ('LK error in calcLKC')

    Tr = thermconsts['LKTr']
    return c1 - (c2/Tr) + (c3/Tr**3)


#
# function calcLKD
def calcD(thermconsts):
    """ Function CalcLKD will calculate the value for D in the
  Lee Kessler EOS.  The function is general in that it
  will calculate either the value for D for the simple
  fluid or the reference fluid depending on the portion
  of the constants array that is passed.  IE, the function
  simply takes the first two doubles from the list and uses
  them for the calculation

      D(Tr) = d1 + d2/Tr

  INPUT:
      DOUBLE      *dvals          array containing the parameters
                                  b1 ... b4
      DOUBLE      Tr              reduced temperature

  OUTPUT:
      DOUBLE      D               value of the function C[returned]
    """
    if thermconsts['LKFluid'] == 'ref':
        d1 = thermconsts['LKd1ref']
        d2 = thermconsts['LKd2ref']
    elif thermconsts['LKFluid'] == 'simple':
        d1 = thermconsts['LKd1simple']
        d2 = thermconsts['LKd2simple']
    else:
        print ('LK error in calcLKD')

    Tr = thermconsts['LKTr']
    return d1 + (d2/Tr)

#
# function calcLKexp
def calcexp(thermconsts):
    """ Function CalcLKexp will calculate the value for
  exponential term in the Lee Kessler EOS.  The function
  is general in that it will calculate either the value
  for the exponential term for the simple fluid or the
  reference fluid depending on the portion
  of the constants array that is passed.  IE, the function
  simply takes the first four doubles from the list and uses
  them for the calculation

      exp term =  [c4/(Tr^3 vr'^2)] [beta + (gamma/vr'^2)] exp(-gamma/vr'^2)

  INPUT:
      DOUBLE      *expvals          array containing the parameters
                                  b1 ... b4
      DOUBLE      Tr              reduced temperature
      DOUBLE      vrp             reduced specific volume as defined by
                                  Lee-Kessler EOS, ie vr' = v/[R Tc / Pc]

  OUTPUT:
      DOUBLE      exp term        value of the exponential term [returned]
    """
    if thermconsts['LKFluid'] == 'ref':
        c4 = thermconsts['LKc4ref']
        beta = thermconsts['LKbetaref']
        gamma = thermconsts['LKgammaref']
    elif thermconsts['LKFluid'] == 'simple':
        c4 = thermconsts['LKc4simple']
        beta = thermconsts['LKbetasimple']
        gamma = thermconsts['LKgammasimple']
    else:
        print ('LK error in calcLKexp')

    Tr = thermconsts['LKTr']
    vrp = thermconsts['LKvrp']
    return c4*(beta + (gamma/vrp**2))*math.exp(-gamma/vrp**2)/(Tr**3*vrp**2)


# 
# calculate the dv integrals, udv = 0, sdv = R ln (v / vcrit)
# note returns 2 values 

def dvint(thermconsts):
    """
        Function dvint calculates the following integrals

                udvint  =       INT {T(dP/dT)|v - P} dv [vcrit,v]
                        =       RT^2 INT {(1/v)(dZ/dT|v)} dv {vcrit, v}

                sdvint  =       INT {(dP/dT)|v dv [vcrit,v]

        based on the Lee Kessler EOS.

        These integrals are used to calculate the thermodynamic
        properties u, h, and s.

        INPUT
                DOUBLE          *consts         critical constants associated with R-503
                DOUBLE          t               temperature                                             K
                DOUBLE          v               specific volume                                         m^3/kg
                DOUBLE          *udvint         see above
                DOUBLE          *sdvint         see above

        OUTPUT
                DOUBLE          *udvint         see above
                DOUBLE          *sdvint         see above

        

    """
    Rgas = thermc.R_BAR/thermconsts['MW']
    T = thermconsts['Temperature']
    P = thermconsts['Pressure']
    v = thermconsts['SpVol']
    Tc = thermconsts['TCrit']
    Pc = thermconsts['PCrit']
    vc = thermconsts['vCrit']

    c4 = thermconsts['LKc4simple']
    beta = thermconsts['LKbetasimple']
    gamma = thermconsts['LKgammasimple']
    
    Tr = T/Tc;
    Tr3 = Tr*Tr*Tr;
    RTP = Rgas*Tc/Pc;
    RTP2 = RTP*RTP;
    RTP4 = RTP2*RTP2;
    A = gamma*RTP2;
    Ainv = 1/A;
    Av = A/(v*v);
    Avc = A/(vc*vc);
    Kappa1 = c4*RTP2*beta/(Tr3);
    Kappa2 = c4*RTP4*gamma/(Tr3);
    Kappa3 = 3*c4*RTP2*beta/(Tr3*T);
    Kappa4 = 3*c4*RTP4*gamma/(Tr3*T);
    vinv = 1/v;
    vinv2 = vinv*vinv;
    vinv5 = vinv2*vinv2*vinv;
    vcinv = 1/vc;
    vcinv2 = vcinv*vcinv;
    vcinv5 = vcinv2*vcinv2*vcinv;
    
    B = calcB(thermconsts)
    C = calcC(thermconsts)
    D = calcD(thermconsts)

    dBdT = calcdBdT(thermconsts)
    dCdT = calcdCdT(thermconsts)
    dDdT = calcdDdT(thermconsts)

    #
    # // calculate portion of du integral
    # //
    # //  du = INT {(1/v)(dZ/dT)|v)} {vc, v}
    # //
    
    du = RTP*(dBdT*(vcinv - vinv) + 0.5*RTP*dCdT*(vcinv2 - vinv2) + 0.2*RTP4*dDdT*(vcinv5 - vinv5)) - \
        Kappa3*(0.5*Ainv*(math.exp(-Av) - math.exp(-Avc))) - Kappa4*(Ainv*(0.5*((vinv2*math.exp(-Av)) - (vcinv2*math.exp(-Avc))) + \
        (0.5*Ainv*(math.exp(-Av) - math.exp(-Avc)))))

    ds = math.log(v/vc) + RTP*((B*(vcinv - vinv)) + (0.5*C*RTP*(vcinv2 - vinv2)) + \
                        (0.2*D*RTP4*(vcinv5 - vinv5))) + \
                        0.5*Kappa1*Ainv*(math.exp(-Av) - math.exp(-Avc)) + \
                        (0.5*Kappa2/A)*((vinv2*math.exp(-Av) - vcinv2*math.exp(-Avc)) + \
                        (0.5*Ainv*(math.exp(-Av) - math.exp(-Avc))));

    udvint = Rgas*T*T*du
    sdvint = Rgas*ds + T*Rgas*du

    return udvint, sdvint


#
# Function calcdBdT
def calcdBdT(thermconsts):
    """
    Function CalcLKdBdT will calculate the derivative of the
  function B(T) at constant v in the Lee Kessler EOS

      dB/dT|v = (b2 Tc/ T^2) + (b3 Tc^2/T^3) + (b4 Tc^3/T^4)

  INPUT:
      DOUBLE      *bvals      parameters in the
                              B function
      DOUBLE      T           temperature             [K]
      DOUBLE      Tc          critical temperature    [K]

  OUTPUT:
      DOUBLE          derivative of function
    """
    if thermconsts['LKFluid'] == 'ref':
        b2 = thermconsts['LKb2ref']
        b3 = thermconsts['LKb3ref']
        b4 = thermconsts['LKb4ref']
    elif thermconsts['LKFluid'] == 'simple':
        b2 = thermconsts['LKb2simple']
        b3 = thermconsts['LKb3simple']
        b4 = thermconsts['LKb4simple']
    else:
        print ('LK error in calcLKB')

    T = thermconsts['Temperature']
    Tc = thermconsts['TCrit']

    t2 = T*T;
    t3 = t2*T;
    t4 = t3*T;

    tc2 = Tc*Tc;
    tc3 = tc2*Tc;

    return (b2*Tc/t2) + (b3*tc2/t3) + (b4*tc3/t4)


#
# function calcLKdCdT
def calcdCdT(thermconsts):
    """  
    Function CalcLKdCdT will calculate the derivative of the
//  function C(T) at constant v in the Lee Kessler EOS
//
//      dC/dT|v = (c2 Tc/ T^2) + (c3 Tc^3/T^4)
//
//  INPUT:
//      DOUBLE      *cvals      parameters in the
//                              C function
//      DOUBLE      T           temperature             [K]
//      DOUBLE      Tc          critical temperature    [K]
//
//  OUTPUT:
//      DOUBLE      dCdt        derivative of function
    """
    if thermconsts['LKFluid'] == 'ref':
        c1 = thermconsts['LKc1ref']
        c2 = thermconsts['LKc2ref']
        c3 = thermconsts['LKc3ref']    
    elif thermconsts['LKFluid'] == 'simple':
        c1 = thermconsts['LKc1simple']
        c2 = thermconsts['LKc2simple']
        c3 = thermconsts['LKc3simple']
    else:
        print ('LK error in calcLKC')

    T = thermconsts['Temperature']
    Tc = thermconsts['TCrit']
    t2 = T*T;
    t4 = t2*t2;
    tc3 = Tc*Tc*Tc;
    
    return (c2*Tc/t2) - (c3*tc3/t4)


#
# function calcdDdT
def calcdDdT(thermconsts):
    """ Function CalcLKD will calculate the value for D in the
  Lee Kessler EOS.  The function is general in that it
  will calculate either the value for D for the simple
  fluid or the reference fluid depending on the portion
  of the constants array that is passed.  IE, the function
  simply takes the first two doubles from the list and uses
  them for the calculation

      D(Tr) = d1 + d2/Tr

  INPUT:
      DOUBLE      *dvals          array containing the parameters
                                  b1 ... b4
      DOUBLE      Tr              reduced temperature

  OUTPUT:
      DOUBLE      D               value of the function C[returned]
    """
    if thermconsts['LKFluid'] == 'ref':
        d2 = thermconsts['LKd2ref']
    elif thermconsts['LKFluid'] == 'simple':
        d2 = thermconsts['LKd2simple']
    else:
        print ('LK error in calcLKD')

    T = thermconsts['Temperature']
    Tc = thermconsts['TCrit']
    t2 = T*T
    return -(d2*Tc/t2)


def vts(thermo):
    thermo['sActual'] = thermo['SpEntropy']
    v0 = 2 * thermo['vCrit']
    thermtuple = 'LK vts',thermo
    #print (thermtuple, type(thermtuple))
    v = scipy.optimize.newton(vfts, v0, None, thermtuple)
    return v

def vfts(v,str, thermo):
    thermo['SpVol'] = v
    #print ('LK in vftp', v)
    thermo['Pressure'] = pvt(thermo)
    (cx, cxint, cxtint) = thermc.Cubiccv(thermo)
    (udvint, sdvint) = thermo['dvint'](thermo)
    scalc = cxtint + sdvint + thermo['sRef'];
    return scalc - thermo['sActual']

def tpv(thermo):
    thermo['pActual'] = thermo['Pressure']
    T0 = ideal.tpv(thermo)
    thermtuple = 'LK tpv',thermo
    #print (thermtuple, type(thermtuple))
    T = scipy.optimize.newton(tfpv, T0, None, thermtuple)
    return T


def tfpv(T,str, thermo):
    thermo['Temperature'] = T
    #print ('LK in vftp', v)
    pcalc = pvt(thermo)
    return pcalc - thermo['pActual']

def tph(thermo):
    thermo['hActual'] = thermo['SpEnthalpy']
    T0 = 0.95*thermo['TCrit']
    thermtuple = 'LK tph', thermo
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
    thermtuple = 'LK tps', thermo
    T = scipy.optimize.newton(tfps, T0, None, thermtuple)
    return T

def tfps(T, str, thermo):
    thermo['Temperature'] = T
    v = vtp(thermo)
    (cx, cxint, cxtint) = thermc.Cubiccv(thermo)
    (udvint, sdvint) = thermo['dvint'](thermo)
    scalc = cxtint + sdvint + thermo['sRef'];
    return scalc - thermo['sActual']
