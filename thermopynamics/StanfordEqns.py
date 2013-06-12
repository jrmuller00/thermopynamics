"""Stanford EOS functions library

// This file contains the implementation of all 
// the equations in the Stanford thermo 
// properties book
//

"""

import ThermoConst as thermc
import math
import scipy.optimize
import IdealGas as ideal
import sys


__version__ = "0.1"


#
# Stanford props eqn driver function
# 

def Stanford(thermconst):
    """
Function Stanford will calculate the
  thermodynamic properties of a fluid based on the
  job number.  The job number will indicate the two
  properties passed to the function.  The equation
  of state and any pertinent functions will be passed
  via the function pointer FunctionHead.  This structure
  is defined in the jrmprops.h file.  Any constants
  associated with the fluid will be in the array consts
  
    
      UNITS:
      Pressure                p                       Pa
      Temperature             t                       K
      Specific volume         v                       m^3/kg
      Specific Energy         u                       J/kg
      Specific Enthalpy       h                       J/kg
      Specific Entropy        s                       J/kg K
      
        
          INPUT:
          INT             job     indicates properties passed to routine
          
            =       1   T,P are inputs
            =       2   T,v are inputs
            =       3   T,s are inputs
            =       4   T,x are inputs
            =       5   P,v are inputs
            =       6   P,h are inputs
            =       7   P,s are inputs
            =       8   P,x are inputs
            =       9   T,h are inputs
            =       10  h,s are inputs
            
              
                DOUBLE          *consts         constants used in EOS
                DOUBLE          *p              pressure
                DOUBLE          *t              temperature
                DOUBLE          *v              specific voulme
                DOUBLE          *u              specific energy
                DOUBLE          *h              specific enthalpy
                DOUBLE          *s              specific entropy
                DOUBLE          *cvr            returned specific heat at constant volume
                DOUBLE          *dpdvr          returned derivative (dP/dv)|T
                DOUBLE          *dpdtr          returned derivative (dP/dT)|v
                
                  
                    FUNCPOINT       *FunctionHead   Function list pointer
    """

    
    job = thermconst['job']
    Tc = thermconst['TCrit']
    Pc = thermconst['PCrit']
    vc = thermconst['vCrit']
    Rgas = thermc.R_BAR/thermconst['MW']
    setStanfordEqns(thermconst)

    #
    # case 1, Inputs are P, T

    if job == 1:
        T = thermconst['Temperature']
        P = thermconst['Pressure']
        vtp = thermconst['vtp']

        #
        # set the region
        if T > Tc:
            region = 'SC'
        else:
            psat = thermconst['psat'](thermconst)
            if P < psat:
                region = 'SHV'
            elif P == psat:
                region = '2phase'
                print ('2 phase region, T,P not independent')
            else:
                region = 'CL'

        thermconst['SpVol'] = ideal.vtp(thermconst)
        thermconst['SpVol']  = vtp(thermconst)
        thermconst['Region'] = region
        calcuhsx(thermconst)



    #
    # case 2, Inputs are T, v

    elif job == 2:
        #
        # Case 2, Inputs are T, v
        T = thermconst['Temperature']
        v = thermconst['SpVol']
        P = thermconst['pvt'](thermconst)
        vtp = thermconst['vtp']

        #
        # set the region
        if T > Tc:
            region = 'SC'
            thermconst['Pressure'] = P
        else:
            psat = thermconst['psat'](thermconst)
            rhosat = thermconst['rhosat'](thermconst)
            vf = 1/rhosat
            #
            # calculate vg
            thermconst['Pressure'] = psat
            thermconst['SpVol'] = ideal.vtp(thermconst)
            vg = vtp(thermconst)
            #
            # reset v
            thermconst['SpVol'] = v
            #
            # reset P
            thermconst['Pressure'] = P
            if v > vg:
                region = 'SHV'
            elif (v > vf) & (v < vg):
                region = '2phase'
                thermconst['Pressure'] = psat
            else:
                region = 'CL'

        thermconst['Region'] = region
        #
        # if 2 phase set variable
        thermconst['2phaseprop'] = 'v'
        calcuhsx(thermconst)

    #
    # case 3, inputs are T, s

    elif job == 3:
        T = thermconst['Temperature']
        s = thermconst['SpEntropy']
        vtp = thermconst['vtp']
                #
        # set the region
        if T > Tc:
            region = 'SC'
            vg = thermconst['vCrit']
        else:
            psat = thermconst['psat'](thermconst)
            rhosat = thermconst['rhosat'](thermconst)
            vf = 1/rhosat
            #
            # calculate vg
            thermconst['Pressure'] = psat
            thermconst['SpVol'] = ideal.vtp(thermconst)
            vg = vtp(thermconst)
            (cx, cxint, cxtint) = thermconst['cvfunc'](thermconst)
            (udvint, sdvint) = thermconst['dvint'](thermconst)
            ug = cxint + udvint + thermconst['xRef']
            sg = cxtint + sdvint + thermconst['sRef']
            hg = ug + psat*vg
            (hfg, sfg) = calchfgsfg(thermconst)
            vfg = vg - vf

            if s > sg:
                region = 'SHV'
                thermconst['SpVol'] = vg
                v = vts(thermconst)
                thermconst['SpVol'] = v
            elif (s > (sg - sfg)):
                region = '2phase'
            else:
                region = 'CL'
        thermconst['Region'] = region
        #
        # if 2 phase set variable
        thermconst['2phaseprop'] = 's'
        calcuhsx(thermconst)

    #
    # case 4, Inputs are T, x

    elif job == 4:
        T = thermconst['Temperature']
        x = thermconst['Quality']
        vtp = thermconst['vtp']
        region = '2phase'
        psat = thermconst['psat'](thermconst)
        thermconst['Pressure'] = psat
        thermconst['SpVol'] = ideal.vtp(thermconst)
        vg = vtp(thermconst)
        thermconst['SpVol'] = vg
        thermconst['Region'] = region
        #
        # if 2 phase set variable
        thermconst['2phaseprop'] = 'x'
        calcuhsx(thermconst)

    #
    # case 5, Inputs are P, v

    elif job == 5:
        #
        # case 5, inputs are P, v
        v = thermconst['SpVol']
        P = thermconst['Pressure']
        Pc = thermconst['PCrit']
        Tc = thermconst['TCrit']
        tpv = thermconst['tpv']
        vtp = thermconst['vtp']
        #tsat = themconsts['tsat']

        #
        # set the region
        if P > Pc:
            region = 'SC'
            thermconst['SpVol'] = v
            thermconst['Temperature'] = ideal.tpv(thermconst)
            T = tpv(thermconst)
            thermconst['Temperature'] = T
        else:
            thermconst['Temperature'] = 0.75*Tc
            tsat = thermconst['tsat'](thermconst)
            thermconst['Temperature'] = tsat
            rhosat = thermconst['rhosat'](thermconst)
            vf = 1/rhosat
            #
            # calculate vg
            vg = vtp(thermconst)
            if v > vg:
                region = 'SHV'
                thermconst['SpVol'] = v
                thermconst['Temperature'] = ideal.tpv(thermconst)
                T = tpv(thermconst)
                thermconst['Temperature'] = T
            elif (v > vf) & (v < vg):
                region = '2phase'
                thermconst['Temperature'] = tsat
            else:
                region = 'CL'
                thermconst['Temperature'] = tsat

        thermconst['Region'] = region
        #
        # if 2 phase set variable
        thermconst['2phaseprop'] = 'v'
        calcuhsx(thermconst)

    #
    # case 6, inputs are P, h
    elif job == 6:
        h = thermconst['SpEnthalpy']
        P = thermconst['Pressure']
        Pc = thermconst['PCrit']
        Tc = thermconst['TCrit']
        vtp = thermconst['vtp']
        thermconst['2phaseprop'] = 'h'

        if P > Pc:
            region = 'SC'
            #
            # use a guess for temperature
            thermconst['Temperature'] = 0.75*Tc
            T = tph(thermconst)
            thermconst['temperature'] = T
            thermconst['SpVol'] = ideal.vtp(thermconst)
            v = vtp(thermconst)
        else:
            thermconst['Temperature'] = 0.75*Tc
            tsat = thermconst['tsat'](thermconst)
            thermconst['Temperature'] = tsat
            rhosat = thermconst['rhosat'](thermconst)
            vf = 1/rhosat
            #
            # calculate vg at Tsat
            thermconst['SpVol'] = ideal.vtp(thermconst)
            vg = vtp(thermconst)
            thermconst['SpVol'] = vg
            thermconst['Region'] = 'SHV'
            calcuhsx(thermconst)
            hg = thermconst['SpEnthalpy']
            #
            # reset h to original value
            thermconst['SpEnthalpy'] = h
            (hfg, sfg) = calchfgsfg(thermconst)
            if h > hg:
                region = 'SHV'
                #
                # use a guess for temperature
                thermconst['Temperature'] = 0.75*Tc
                T = tph(thermconst)
                thermconst['Temperature'] = T
                thermconst['SpVol'] = ideal.vtp(thermconst)
                v = vtp(thermconst)
            elif (h > (hg - hfg)):
                region = '2phase'
                thermconst['Temperature'] = tsat
                thermconst['SpVol'] = vg
            else:
                region = 'CL'
                thermconst['Temperature'] = tsat
                thermconst['SpVol'] = vf

        thermconst['Region'] = region        
        calcuhsx(thermconst)


    #
    # case 7, inputs are P, s
    elif job == 7:
        s = thermconst['SpEntropy']
        P = thermconst['Pressure']
        Pc = thermconst['PCrit']
        Tc = thermconst['TCrit']
        vtp = thermconst['vtp']
        thermconst['2phaseprop'] = 's'

        #
        # set the region
        if P > Pc:
            region = 'SC'
            #
            # use a guess for temperature
            thermconst['Temperature'] = 0.75*Tc
            T = tps(thermconst)
            thermconst['temperature'] = T
            thermconst['SpVol'] = ideal.vtp(thermconst)
            v = vtp(thermconst)
        else:
            thermconst['Temperature'] = 0.75*Tc
            tsat = thermconst['tsat'](thermconst)
            thermconst['Temperature'] = tsat
            rhosat = thermconst['rhosat'](thermconst)
            vf = 1/rhosat
            #
            # calculate vg
            thermconst['SpVol'] = ideal.vtp(thermconst)
            vg = vtp(thermconst)
            thermconst['Region'] = 'SHV'
            calcuhsx(thermconst)
            sg = thermconst['SpEntropy']
            #
            # reset s to original value
            thermconst['SpEntropy'] = s
            (hfg, sfg) = calchfgsfg(thermconst)
            if s > sg:
                region = 'SHV'
                #
                # use a guess for temperature
                thermconst['Temperature'] = 0.75*Tc
                T = tps(thermconst)
                thermconst['Temperature'] = T
                thermconst['SpVol'] = ideal.vtp(thermconst)
                v = vtp(thermconst)
            elif (s > (sg - sfg)):
                region = '2phase'
                thermconst['Temperature'] = tsat
                thermconst['SpVol'] = vg
            else:
                region = 'CL'
                thermconst['Temperature'] = tsat
                thermconst['SpVol'] = vg
        
        thermconst['Region'] = region        
        calcuhsx(thermconst)

    #
    # case 8, Inputs are P, x

    elif job == 8:
        P = thermconst['Pressure']
        x = thermconst['Quality']
        vtp = thermconst['vtp']
        region = '2phase'
        thermconst['Temperature'] = 0.75*Tc
        tsat = thermconst['tsat'](thermconst)
        thermconst['Temperature'] = tsat
        #
        # calculate vg
        thermconst['SpVol'] = ideal.vtp(thermconst)
        vg = vtp(thermconst)
        thermconst['SpVol'] = vg
        thermconst['Region'] = region
        #
        # if 2 phase set variable
        thermconst['2phaseprop'] = 'x'
        calcuhsx(thermconst)

    #
    # case 9, inputs are T, h

    elif job == 9:
        T = thermconst['Temperature']
        h = thermconst['SpEnthalpy']
        vtp = thermconst['vtp']
                #
        # set the region
        if T > Tc:
            region = 'SC'
            vg = thermconst['vCrit']
        else:
            psat = thermconst['psat'](thermconst)
            rhosat = thermconst['rhosat'](thermconst)
            vf = 1/rhosat
            #
            # calculate vg
            thermconst['Pressure'] = psat
            thermconst['SpVol'] = ideal.vtp(thermconst)
            vg = vtp(thermconst)
            (cx, cxint, cxtint) = thermconst['cvfunc'](thermconst)
            (udvint, sdvint) = thermconst['dvint'](thermconst)
            ug = cxint + udvint + thermconst['xRef']
            sg = cxtint + sdvint + thermconst['sRef']
            hg = ug + psat*vg
            (hfg, sfg) = calchfgsfg(thermconst)
            print ('Job 9: T = ',T,' h = ',h,' hg = ',hg)
            
            if h > hg:
                region = 'SHV'
                thermconst['SpVol']  = vg
                v = vth(thermconst)
                thermconst['SpVol'] = v
            elif (h > (hg - hfg)):
                region = '2phase'
            else:
                region = 'CL'
        thermconst['Region'] = region
        #
        # if 2 phase set variable
        thermconst['2phaseprop'] = 'h'
        calcuhsx(thermconst)

    return

#
#  P-rho-T EQUATIONS
#

#
# P-1

def pvt1(thermconst):
    """
//
//
//
// P =  R T/(v-b) + SUM{i=2,5} [(1 /(v-b)^i) (A_i + B_i T + C_i exp(-k T/ Tc))]
//              + ((A_6 + B_6 T + C_6 exp(-k T/ Tc)/(exp(a v) (1 + c exp(a v)))
        //
        // DECLARATIONS
        //
        //      INTEGERS
        //      
        //              i                       loop counter, index for arrays
        //              offsetA         offset to A parameters
        //              offsetB         offset to B parameters
        //              offsetC         offset to C parameters
        //
        //      DOUBLES
        //
        //              Tc                      critical temperature
        //              R                       gas constant
        //              A6                      equation parameter
        //              B6                      equation parameter
        //              C6                      equation parameter
        //              b                       equation parameter
        //              kappa           equation parameter
        //              alpha           equation parameter
        //              P                       pressure
        //              term            intermediate calculation term
        //              term2           intermediate calculation term
        //              expterm1        exponential term exp(-kT/Tc)
        //              expterm2        exponential term exp(a v)
        //

    """
    #
    # create the A,B, C arrays

    A = []
    B = []
    C = []
    #for i in range(2, 7):
    #    key = 'StanfordA' + str(i)
    #    A.append(thermconst[key])
    #    key = 'StanfordB' + str(i)
    #    B.append(thermconst[key])
    #    key = 'StanfordC' + str(i)
    #    C.append(thermconst[key])

    loadArray(thermconst,A,'StanfordA',1,7)
    loadArray(thermconst,B,'StanfordB',1,7)
    loadArray(thermconst,C,'StanfordC',1,7)

    try:
        A6 = thermconst['StanfordA6']
    except:
        A6 = 0
    try:
        B6 = thermconst['StanfordB6']
    except:
        B6 = 0
    try:
        C6 = thermconst['StanfordC6']
    except:
        C6 = 0
    try:
        b = thermconst['StanfordAb']
    except:
        b = 0
    try:
        c = thermconst['StanfordAc']
    except:
        c = 0

    try:
        alpha = thermconst['StanfordAalpha']
    except:
        alpha = 0
    try:
        kappa = thermconst['StanfordAkappa']
    except:
        kappa = 0
        
    # print (A)

    v = thermconst['SpVol']
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    R = thermc.R_BAR/thermconst['MW']

    P = R * T /(v - b)

    expterm = math.exp(-kappa * T / Tc)
    expterm2 = math.exp(-alpha*v)
    for i in range(1,5):
        P = P + (A[i] + B[i] * T + C[i] * expterm)/((v - b)**(i+1))


    term = 0

    if (c == 0):
        term = expterm2*(A6 + B6*T + C6*expterm)
    else:
        term = expterm2*expterm2*(A6 + B6*T + C6*expterm)/(expterm2+c)

        P = P + term
    return P


def dvint1(thermconst):
    """

    """

    #
    # create the A,B, C arrays

    A = []
    B = []
    C = []
    #for i in range(2, 7):
    #    key = 'StanfordA' + str(i)
    #    A.append(thermconst[key])
    #    key = 'StanfordB' + str(i)
    #    B.append(thermconst[key])
    #    key = 'StanfordC' + str(i)
    #    C.append(thermconst[key])

    loadArray(thermconst,A,'StanfordA',2,7)
    loadArray(thermconst,B,'StanfordB',2,7)
    loadArray(thermconst,C,'StanfordC',2,7)

    try:
        A6 = thermconst['StanfordA6']
    except:
        A6 = 0
    try:
        B6 = thermconst['StanfordB6']
    except:
        B6 = 0
    try:
        C6 = thermconst['StanfordC6']
    except:
        C6 = 0
    try:
        b = thermconst['StanfordAb']
    except:
        b = 0
    try:
        c = thermconst['StanfordAc']
    except:
        c = 0

    try:
        alpha = thermconst['StanfordAalpha']
    except:
        alpha = 0
    try:
        kappa = thermconst['StanfordAkappa']
    except:
        kappa = 0

        
    v = thermconst['SpVol']
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    R = thermc.R_BAR/thermconst['MW']

    #
    # calculate exponential terms

    expterm1 = math.exp(-kappa*T/Tc)
    expterm2 = math.exp(-alpha*v)
    #expterm3 = 1/expterm2

    #//
    #// calculate the du integral

    #//
    #// series term
        
    term = 0
    term2 = (T*kappa/Tc) + 1
        
    for i in range(0,4):
        term = term + (A[i] + C[i] * term2 * expterm1)/ ((i+1)*((b-v)**(i+1)))
        

    udvint = -term

    #//
    #// Now, add the final term based on the value of alpha

    if (alpha != 0):
        if ( c != 0):
            term = c*math.log(1 + (expterm2/c)) - expterm2
            term = term/alpha
        else:
             term = -expterm2/alpha
        
        udvint = udvint + (-A6 - (C6*expterm1*((T*kappa/Tc) + 1)))*term



    #//
    #// Now do the entropy integral

    #//
    #// series term
        
    term = 0
    term2 = (T*kappa/Tc) + 1
        
    for i in range(0,4):
        term = term + (B[i] - C[i] * (kappa/Tc) * expterm1)/ ((i+1)*((b-v)**(i+1)))

    sdvint = term;

    term2 = B6 + (C6*(-kappa/Tc)*expterm1);

    #//
    #// Now, add the final term based on the value of alpha

    if (alpha != 0):
        if ( c != 0):
            term = c*math.log(1 + (expterm2/c)) - expterm2
            term = term/alpha
        else:
             term = -expterm2/alpha

    sdvint = sdvint + term2 * term
        
    sdvint = sdvint + (R*math.log(v-b));
               
    return udvint, sdvint 

#
# P-2
def pvt2(thermconst):
    """

 STANFORD PRESSURE-DENSITY-TEMPERATURE EQUATION (P-2)



 P =  rho R T + [B_0 R T - A_0 - (C_0/T^2) + (D_0/T^3) - (E_0/T^4)] rho^2

                      + [b R T - a - (d/T)] rho^3 + alpha  [a + (d/T)] rho^6

                      + c (rho^3 / T^2) [1 + (gamma rho^2)] exp (-gamma rho^2)

        //
        // DECLARATIONS
        //
        //      INTEGERS
        //      
        //              i               loop counter, index for arrays
        //
        //      DOUBLES
        //
        //              Tc              critical temperature
        //              R               gas constant
        //              A0              equation constant
        //              B0              equation constant
        //              C0              equation constant
        //              D0              equation constant
        //              E0              equation constant
        //              a               equation constant
        //              b               equation constant
        //              c               equation constant
        //              d               equation constant
        //              alpha   equation constant
        //              gamma   equation constant
        //              P               pressure
        //              term    intermediate calculation term
        //              expterm intermediate calculation term
        //              rho3    rho^3
        //              rho4    rho^4
        //              rho5    rho^5
        //              rho6    rho^6

    """
    rho = 1.0 / thermconst['SpVol']
    T = thermconst['Temperature']
    Rgas = thermc.R_BAR/thermconst['MW']

    A0 = thermconst['StanfordAA0']
    B0 = thermconst['StanfordAB0']
    C0 = thermconst['StanfordAC0']
    D0 = thermconst['StanfordAD0']
    E0 = thermconst['StanfordAE0']
    a = thermconst['StanfordAa']
    b = thermconst['StanfordAb']
    c = thermconst['StanfordAc']
    d = thermconst['StanfordAd']
    alpha = thermconst['StanfordAalpha']
    gamma = thermconst['StanfordAgamma']

    rho2 = rho*rho
    rho3 = rho2*rho
    rho4 = rho3*rho
    rho5 = rho4*rho
    rho6 = rho5*rho
    T2 = T*T
    T3 = T2*T
    T4 = T3*T

    #
    # 1st term in equation

    P = rho*Rgas*T

    #
    # 2nd term in equation

    term = (B0*Rgas*T - A0 - (C0/T2) + (D0/T3) - (E0/T4))

    P = P + term*rho2

    #
    # 3rd term in equation

    term = (b*Rgas*T - a - (d/T))

    P = P + term*rho3

    #
    # 4th term in equation
        
    P = P + alpha*(a + (d/T))*rho6

    #
    # 5th term in equation

    expterm = math.exp(-gamma*rho2)

    P = P + (c*rho3*(1 + gamma*rho2)*expterm/T2)

    return P

def dvint2(thermconst):
    """
    // Function dvint2 calculates the
// integrals for the "dv" or "d rho" portion 
// of the internal energy and the entropy
// calculations.
    """

    rho = 1.0 / thermconst['SpVol']
    T = thermconst['Temperature']
    Rgas = thermc.R_BAR/thermconst['MW']

    A0 = thermconst['StanfordAA0']
    B0 = thermconst['StanfordAB0']
    C0 = thermconst['StanfordAC0']
    D0 = thermconst['StanfordAD0']
    E0 = thermconst['StanfordAE0']
    a = thermconst['StanfordAa']
    b = thermconst['StanfordAb']
    c = thermconst['StanfordAc']
    d = thermconst['StanfordAd']
    alpha = thermconst['StanfordAalpha']
    gamma = thermconst['StanfordAgamma']

    rho2 = rho*rho
    rho3 = rho2*rho
    rho4 = rho3*rho
    rho5 = rho4*rho
    rho6 = rho5*rho
    T2 = T*T
    T3 = T2*T
    T4 = T3*T
    T5 = T4*T

    #
    # Calculate the u dv integral
        
    term = (-A0 - (3*C0/T2) + (4*D0/T3) - (5*E0/T4))

    udvint = term*rho

    term = (-a - (2*d/T))

    udvint = udvint + (term*rho2/2)

    udvint = udvint + (alpha*(a + (2*d/T))*rho5/5)

    expterm = math.exp(-gamma*rho2)

    udvint = udvint + (3*c*(1 - expterm)/(2*T2*gamma))

    udvint = udvint + (3*c*gamma*((-rho2*expterm/(2*gamma)) + ((expterm - 1)/(2*gamma*gamma)))/T2)

    #
    # now calcualte the entropy integral

    term = -(B0*Rgas + (2*C0/T3) - (3*D0/T4) + (4*E0/T5))
        
    sdvint = term*rho

    term = b*Rgas + (d/T2)
        
    sdvint = sdvint - 0.5*term*rho2

    sdvint = sdvint + (alpha*d*rho5/(5*T2))

    expterm = math.exp(-gamma*rho2)

    sdvint = sdvint + (2*c*(1 - expterm)/(2*T3*gamma))

    sdvint = sdvint + (2*c*gamma*((-rho2*expterm/(2*gamma)) + ((expterm - 1)/(2*gamma*gamma)))/T3)

    sdvint = sdvint - Rgas*math.log(rho)

    return udvint, sdvint


#
# P-3

def pvt3(thermconst):
    """
// STANFORD PRESSURE-DENSITY-TEMPERATURE EQUATION (P-3)
//
//
//
// P =  rho R T + rho^2 SUM[i=1,5]{A_i T^(2-i)} + rho^3 SUM[i=6,8]{A_i T^(7-i)}
//              
//              + rho^4 (A_9 T + A_10) + rho^5 (A_11 T + A_12) + rho^6 A_13
//
//              + { rho^3 SUM[i=14,16]{A_i T^(12-i)} 
//
//              + rho^5 SUM[i=17,19]{A_i T^(15-i)} } exp(-gamma rho^2)
// 
//
// INPUT:
//      
//              double                  *consts         constants for the equations
//              double                  rho                     density of the fluid            [kg/m^3]
//              double                  T                       temperature of the fluid        [K]
//              FunctionPointer *List           pointer to full equation list
//
//      OUTPUT:
//              
//              VOID
//
//      RETURN:
//
//              double                  P                       pressure of the fluid           [Pa]
    """

    #
    # create the A array

    A = []
    #for i in range(1, 20):
    #    key = 'StanfordA' + str(i)
    #    A.append(thermconst[key])
    loadArray(thermconst, A, 'StanfordA', 1, 20)
    # print (Atest)

    A9 = thermconst['StanfordA9']
    A10 = thermconst['StanfordA10']
    A11 = thermconst['StanfordA11']
    A12 = thermconst['StanfordA12']
    A13 = thermconst['StanfordA13']
    gamma = thermconst['StanfordAgamma']
        
    # print (A)

    rho = 1.0 / thermconst['SpVol']
    T = thermconst['Temperature']
    Rgas = thermc.R_BAR/thermconst['MW']

    rho2 = rho*rho;
    rho3 = rho2*rho;
    rho4 = rho3*rho;
    rho5 = rho4*rho;
    rho6 = rho5*rho;

    #
    # 1st term in equation

    P = rho*Rgas*T

    #
    # 2nd term in equation

    term = 0;
    
    for i in range(0,5):
        term = term + A[i] * (T**(1-i))
        #print ('2nd term = ',i)
    term = term * rho2
    P = P + term

    #
    # 3rd term in equation

    term = 0;

    for i in range(5, 8):
        term = term + A[i] * (T**(6-i))
        #print ('3rd term = ',i)
    term = term * rho3
    P = P + term
    
    #
    # 4th term in equation

    P = P + rho4*(A9*T + A10)

    #
    # 5th term in equation

    P = P + rho5*(A11*T + A12)

    #
    # 6th term in equation

    P = P + rho6*A13

    #
    # 7th term in equation

    term = 0;
    for i in range(13,16):
        term = term + A[i] * (T**(11-i))
        
    term = rho3*term

    term2 = 0;

    for i in range(16,19):
        term2 = term2 + A[i] * (T**(14-i))

    term2 = rho5*term2

    P = P + math.exp(-gamma*rho2)*(term + term2)

    return P



def dvint3(thermconst):
    """
// Function StanfordPvTdv3 calculates the
// integrals for the "dv" or "d rho" portion 
// of the internal energy and the entropy
// calculations.
//
// INPUT:
//              DOUBLE          *consts         critical constants
//                                                                              associated with fluid
//              DOUBLE          t               temperature                             [K]
//              DOUBLE          rho             density                                 [kg/m^3]
//              DOUBLE          *udvint         pointer to internal 
//                                                                              energy integral
//              DOUBLE          *sdvint         pointer to entropy 
//                                                                              integral
// OUTPUT
//              DOUBLE          *udvint         internal energy 
//                                                                              integral                                [J/kg]
//              DOUBLE          *sdvint         entropy integral                [J/kg K]
//
//
//              
//
//      NOTE:   The "d rho" integral vanishes at rho = 0.  This is 
//                      so the equation will approach ideal gas behavior
//                      at low density.  Therefore, the derived integrals
//                      are evaluated at rho only since, at rho = 0, ideal 
//                      gas behavior predicts P = 0.
    """

       #
    # create the A array

    A = []
    #for i in range(1, 20):
    #    key = 'StanfordA' + str(i)
    #    A.append(thermconst[key])

    loadArray(thermconst,A,'StanfordA',1,20)

    A9 = thermconst['StanfordA9']
    A10 = thermconst['StanfordA10']
    A11 = thermconst['StanfordA11']
    A12 = thermconst['StanfordA12']
    A13 = thermconst['StanfordA13']
    gamma = thermconst['StanfordAgamma']
        
    # print (A)

    rho = 1.0 / thermconst['SpVol']
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    Rgas = thermc.R_BAR/thermconst['MW']

    rho2 = rho*rho;
    rho3 = rho2*rho;
    rho4 = rho3*rho;
    rho5 = rho4*rho;
    rho6 = rho5*rho;

    #
    # Calculate the first integral, Pint

    Pint = 0.0

    #
    # 2nd term in equation

    term = 0.0

    for i in range(0,5):
        term = term + i * A[i] * T**(1-i)

    term = term * rho

    Pint = Pint + term

    #
    # 3rd term in equation

    term = 0.0
    for i in range(5,9):
        term = term + (i - 5) * A[i] * (T**(6-i))

    term = term * ( rho2 / 2.0)

    Pint = Pint + term 

    #
    # 4th term in equation

    Pint = Pint + ((rho3*A10)/3.0)

    #
    # 5th term in equation

    Pint = Pint + ((rho4*A12)/4.0)

    #
    # 6th term in equation

    Pint = Pint + (A13*rho5/5.0)

    #
    # Calculate expterm

    expterm = math.exp(-gamma*rho2)
    #
    # 7th term in equation

    term = 0
    for i in range(13,16):
        term = term + A[i] * (i-10)*(T**(11-i))

    Pint = Pint + term * (1/(2*gamma)) * ( 1 - expterm)

    term2 = 0
    for i in range(16,19):
        term2 = term2 + A[i] * (i-13) * (T**(14-i))

    Pint = Pint + term2*((-rho2/(2*gamma))*expterm + (1/gamma) * ((1/(2*gamma)) * ( 1 - expterm)))
    udv = Pint

    sdv = 0

    #
    # 2nd term in equation

    term = 0.0

    for i in range(0,5):
        term = term + (1 - i) * A[i] * T**(-i)

    sdv = rho * term

    #
    # 3rd term in equation

    term = 0.0
    for i in range(5,9):
        term = term + (6 - i) * A[i] * (T**(5-i))

    term = term * ( rho2 / 2.0)

    sdv = sdv + term

    #
    # 4th term in equation

    sdv = sdv + A9*(rho3/3);

    #
    # 5th term in equation

    sdv = sdv + A11*(rho4/4);

        #
    # 7th term in equation

    term = 0
    for i in range(13,16):
        term = term + A[i] * (11 - i)*(T**(10-i))

    sdv = sdv + term * (1/(2*gamma)) * ( 1 - expterm)

    term2 = 0
    for i in range(16,19):
        term2 = term2 + A[i] * (14 - i) * (T**(13-i))

    sdv = sdv + term2*((-rho2/(2*gamma))*expterm + (1/gamma) * ((1/(2*gamma)) * ( 1 - expterm)))
    
    sdv = -Rgas * math.log(rho) - sdv

    return udv, sdv


#
# P-4
#
def pvt4(thermconst):
    """
    //
// STANFORD PRESSURE-DENSITY-TEMPERATURE EQUATION (P-4)
//
//
//
// P =  rho R T + rho^2 [A_1 T + A_2 T^(1/2) + SUM[i=3,5]{A_i T^(3-i)} + rho^3 SUM[i=6,9]{A_i T^(7-i)}
//              
//              + rho^4  SUM[i=10,12]{A_i T^(11-i)} + rho^5 (A_13) + rho^6 (A_14/T + A_15/T^2)
//
//              + rho^7 A_16/T + rho^8 (A_17/T + A_18/T^2) + rho^9 A19/T^2
//
//              + { rho^3 (A_20/T^2 + A_21/T^3) + rho^5 (A_22/T^2 + A_23/T^4)
//      
//              +   rho^7 (A_24/T^2 + A_25/T^3) + rho^9 (A_26/T^2 + A_27/T^4)
//
//              +       rho^11 (A_28/T^2 + A_29/T^3)
//
//              +       rho^13 (A_30/T^2 + A_31/T^3 + A_32/T^4)} exp(-gamma rho^2)
//

/////////////////////////////////////////////////
// Function StanfordPvT3 is the implemntation 
// of the P-rho-t equation (P-4) in the Stanford
// book.  This version of the equation lists
// the density first and the temperature second.



    """
    #
    # create the A array

    A = []
    #for i in range(1, 20):
    #    key = 'StanfordA' + str(i)
    #    A.append(thermconst[key])
    loadArray(thermconst, A, 'StanfordA', 1, 33)
    # print (Atest)

    A1 = thermconst['StanfordA1']
    A2 = thermconst['StanfordA2']

    A13 = thermconst['StanfordA13']
    A14 = thermconst['StanfordA14']
    A15 = thermconst['StanfordA15']
    A16 = thermconst['StanfordA16']
    A17 = thermconst['StanfordA17']
    A18 = thermconst['StanfordA18']
    A19 = thermconst['StanfordA19']
    A20 = thermconst['StanfordA20']
    A21 = thermconst['StanfordA21']
    A22 = thermconst['StanfordA22']
    A23 = thermconst['StanfordA23']
    A24 = thermconst['StanfordA24']
    A25 = thermconst['StanfordA25']
    A26 = thermconst['StanfordA26']
    A27 = thermconst['StanfordA27']
    A28 = thermconst['StanfordA28']
    A29 = thermconst['StanfordA29']
    A30 = thermconst['StanfordA30']
    A31 = thermconst['StanfordA31']
    A32 = thermconst['StanfordA32']
    gamma = thermconst['StanfordAgamma']
        
    # print (A)

    rho = 1.0 / thermconst['SpVol']
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    Rgas = thermc.R_BAR/thermconst['MW']

    #//
    #// calcualte powers of rho

    rho2 = rho*rho
    rho3 = rho2*rho
    rho4 = rho3*rho
    rho5 = rho4*rho
    rho6 = rho5*rho
    rho7 = rho6*rho
    rho8 = rho7*rho
    rho9 = rho8*rho
    rho11 = rho9*rho2
    rho13 = rho11*rho2

    #//
    #// powers of T

    T2 = T*T
    T3 = T2*T
    T4 = T3*T

    #//
    #// 1st term in equation

    P = rho*Rgas*T

    #//
    #// 2nd term in equation

    term = A1*T + A2*math.sqrt(T)

    for i in range(2,5):
        term = term + A[i]*(T**(2-i))

    P = P + term * rho2

    #//
    #// 3rd term in equation

    term = 0

    for i in range(5,9):
        term = term + A[i] * (T**(6-i))

    term = term*rho3

    P = P + term

    #//
    #// 4th term in equation

    term = 0

    for i in range(9,12):
        term = term + A[i] * (T**(10-i))

    term = term*rho4

    P = P + term

    #//
    #// 5th term in equation

    P = P + rho5*A13

    #//
    #// 6th term in equation

    P = P + rho6*((A14/T) + (A15/T2))

    #//
    #// 7th term in equation

    P = P + (rho7*A16/T)

    #//
    #// 8th term

    P = P + (rho8*((A17/T) + (A18/T2)))

    #//
    #// 9th term

    P = P + (rho9*A19/T2)

    #//
    #// exponential term

    expterm = math.exp(-gamma*rho2)

    term = rho3*((A20/T2) + (A21/T3))
    term = term + rho5*((A22/T2) + (A23/T4))
    term = term + rho7*((A24/T2) + (A25/T3))
    term = term + rho9*((A26/T2) + (A27/T4))
    term = term + rho11*((A28/T2) + (A29/T3))
    term = term + rho13*((A30/T2) + (A31/T3) + (A32/T4))

    P = P + expterm*term

    return P

def dvint4(thermconst):
    """
    /////////////////////////////////////////////////
// Function StanfordPvTdv4 calculates the
// integrals for the "dv" or "d rho" portion 
// of the internal energy and the entropy
// calculations.

    """
    #
    # create the A array

    A = []
    #for i in range(1, 20):
    #    key = 'StanfordA' + str(i)
    #    A.append(thermconst[key])
    loadArray(thermconst, A, 'StanfordA', 1, 33)
    # print (Atest)

    A1 = thermconst['StanfordA1']
    A2 = thermconst['StanfordA2']

    A13 = thermconst['StanfordA13']
    A14 = thermconst['StanfordA14']
    A15 = thermconst['StanfordA15']
    A16 = thermconst['StanfordA16']
    A17 = thermconst['StanfordA17']
    A18 = thermconst['StanfordA18']
    A19 = thermconst['StanfordA19']
    A20 = thermconst['StanfordA20']
    A21 = thermconst['StanfordA21']
    A22 = thermconst['StanfordA22']
    A23 = thermconst['StanfordA23']
    A24 = thermconst['StanfordA24']
    A25 = thermconst['StanfordA25']
    A26 = thermconst['StanfordA26']
    A27 = thermconst['StanfordA27']
    A28 = thermconst['StanfordA28']
    A29 = thermconst['StanfordA29']
    A30 = thermconst['StanfordA30']
    A31 = thermconst['StanfordA31']
    A32 = thermconst['StanfordA32']
    gamma = thermconst['StanfordAgamma']
    twogamma = 2 * gamma
        
    # print (A)

    rho = 1.0 / thermconst['SpVol']
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    Rgas = thermc.R_BAR/thermconst['MW']

    #//
    #// calcualte powers of rho

    rho2 = rho*rho
    rho3 = rho2*rho
    rho4 = rho3*rho
    rho5 = rho4*rho
    rho6 = rho5*rho
    rho7 = rho6*rho
    rho8 = rho7*rho
    rho9 = rho8*rho
    rho10 = rho9*rho
    rho11 = rho10*rho
    rho12 = rho11*rho
    rho13 = rho11*rho2

    #//
    #// powers of T

    T2 = T*T
    T3 = T2*T
    T4 = T3*T
    T5 = T4*T

    #//
    #// exponential term

    expterm = math.exp(-gamma*rho2)
    expotgam = expterm/twogamma

    C = []
    Cp = []
    I = []

    #//
    #// u dv integral

    C.append(A1*T + A2*math.sqrt(T))

    Cp.append(A1 + (0.5*A2/math.sqrt(T)))

    for i in range(2,5):
        C[0] = C[0] + (A[i]*(T**(2-i)))
        Cp[0] = Cp[0] + (2-i)*(A[i]*(T**(1-i)))


    C.append(0)
    Cp.append(0)

    
    for i in range(5,9):
        C[1] = C[1] + (A[i]*(T**(6-i)))
        Cp[1] = Cp[1] + (6-i)*(A[i]*(T**(5-i)))

    C.append(0) 
    Cp.append(0)

    for i in range(9,12):
        C[2] = C[2] + (A[i]*(T**(10-i)))
        Cp[2] = Cp[2] + (10-i)*(A[i]*(T**(9-i)))

    C.append(A13)
    Cp.append(0)

    C.append((A14/T) + (A15/T2))
    Cp.append((-A14/T2) - (2*A15/T3))

    C.append(A16/T)
    Cp.append(-A16/T2)

    C.append((A17/T) + (A18/T2))
    Cp.append(-(A17/T2) - (2*A18/T3))

    C.append(A19/T2)
    Cp.append(-2*A19/T3)
    
    C.append((A20/T2) + (A21/T3))
    Cp.append(-(2*A20/T3) - (3*A21/T4))

    C.append((A22/T2) + (A23/T4))
    Cp.append(-(2*A22/T3) - (4*A23/T5))

    C.append((A24/T2) + (A25/T3))
    Cp.append(-(2*A24/T3) - (3*A25/T4))

    C.append((A26/T2) + (A27/T4))
    Cp.append((-2*A26/T3) - (4*A27/T5))

    C.append((A28/T2) + (A29/T3))
    Cp.append(-(2*A28/T3) - (3*A29/T4))

    C.append((A30/T2) + (A31/T3) + (A32/T4))
    Cp.append(-(2*A30/T3) - (3*A31/T4) - (4*A32/T5))


    I.append(rho)
    I.append(0.5*rho2)
    I.append(rho3/3)
    I.append(0.25*rho4)
    I.append(0.2*rho5)
    I.append(rho6/6)
    I.append(rho7/7)
    I.append(rho8/8)
    I.append((1/(2*gamma))-(expotgam))

    I.append(-(rho2*expotgam) + (I[8]/gamma))

    I.append(-(rho4*expotgam) + (2*I[9]/gamma))

    I.append(-(rho6*expotgam) + (3*I[10]/gamma))
    I.append(-(rho8*expotgam) + (4*I[11]/gamma))

    I.append(-(rho10*expotgam) + (5*I[12]/gamma))

    udvint = 0
    sdvint = 0

    for i in range(0,14):
        udvint = udvint + (C[i] - T*Cp[i])*I[i]
        sdvint = sdvint + Cp[i]*I[i]

    sdvint = -Rgas*math.log(rho) - (sdvint);
    return udvint, sdvint

#
# PVT - 5
#

def pvt5(thermconst):
    """
//
// STANFORD PRESSURE-DENSITY-TEMPERATURE EQUATION (P-5)
//
//
//
// P =  rho R T + rho^2 SUM[i=1,9]{A_i T^(2 - ((i-1)/2))} + rho^3 SUM[i=10,17]{A_i T^(1-((i-10)/2))}
//              
//              + rho^4  SUM[i=18,21]{A_i T^((1/2)-(i-18))} + rho^5 SUM[i=22,27]{A_i T^((1/2)-((i-22)/4)))}
//
//              + rho^6 (A_28 + (A_29/T))
//
//              + { rho^3 (A_30 + A_31/T + A_32/T^2) + rho^5 (A_33 + A_34/T + A_35/T^2)} exp(-gamma rho^2)
//
    """
    #
    # create the A array

    A = []
    loadArray(thermconst, A, 'StanfordA', 1, 33)

    A28 = thermconst['StanfordA28']
    A29 = thermconst['StanfordA29']
    A30 = thermconst['StanfordA30']
    A31 = thermconst['StanfordA31']
    A32 = thermconst['StanfordA32']
    A33 = thermconst['StanfordA33']
    A34 = thermconst['StanfordA34']
    A35 = thermconst['StanfordA35']
    gamma = thermconst['StanfordAgamma']
        
    # print (A)

    rho = 1.0 / thermconst['SpVol']
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    Rgas = thermc.R_BAR/thermconst['MW']

    #//
    #// calcualte powers of rho

    rho2 = rho*rho
    rho3 = rho2*rho
    rho4 = rho3*rho
    rho5 = rho4*rho
    rho6 = rho5*rho

    #//
    #// powers of T

    T2 = T*T

    #//
    #// 1st term in equation

    P = rho*Rgas*T

    for i in range (0,9):
        term = term + (A[i])/(T**((i-1)/2)-2)


    term = term*rho2

    P = P + term

    #//
    #// 3rd term in equation

    term = 0

    for i in range(9,17):
        term = term + (A[i])*(T**(1-((i-10)/2)))
  

    term = term*rho3

    P = P + term

    #//
    #// 4th term in equation

    term = 0

    for i in range(17,21):
            term = term + (A[i])*(T**(0.5-(i-18)))

    term = term*rho4

    P = P + term


    #//
    #// 5th term in equation

    term = 0

    for i in range(21,27):
        term = term + (A[i])*(T**(0.5-((i-22)/4)))

    term = term*rho5

    P = P + term

    #//
    #// 6th term in equation

    P = P + rho6*(A28 + (A29/T2))

    #//
    #// exponential term

    expterm = math.exp(-gamma*rho2)

    term = rho3*(A30 + (A31/T) + (A32/T2))
    term = term + rho5*(A33 + (A34/T) + (A35/T2))

    P = P + expterm*term

    return P

def dvint5(thermconst):
    """
/////////////////////////////////////////////////
// Function dvint5 calculates the
// integrals for the "dv" or "d rho" portion 
// of the internal energy and the entropy
// calculations.
//        
//
//      NOTE:   The "d rho" integral vanishes at rho = 0.  This is 
//                      so the equation will approach ideal gas behavior
//                      at low density.  Therefore, the derived integrals
//                      are evaluated at rho only since, at rho = 0, ideal 
//                      gas behavior predicts P = 0.
    """



    #
    # create the A array

    A = []
    loadArray(thermconst, A, 'StanfordA', 1, 33)

    A28 = thermconst['StanfordA28']
    A29 = thermconst['StanfordA29']
    A30 = thermconst['StanfordA30']
    A31 = thermconst['StanfordA31']
    A32 = thermconst['StanfordA32']
    A33 = thermconst['StanfordA33']
    A34 = thermconst['StanfordA34']
    A35 = thermconst['StanfordA35']
    gamma = thermconst['StanfordAgamma']
        
    # print (A)

    rho = 1.0 / thermconst['SpVol']
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    Rgas = thermc.R_BAR/thermconst['MW']

    #//
    #// calcualte powers of rho

    rho2 = rho*rho
    rho3 = rho2*rho
    rho4 = rho3*rho
    rho5 = rho4*rho
    rho6 = rho5*rho

    #//
    #// powers of T

    T2 = T*T
    T3 = T*T2
    T4 = T*T3

    #//
    #// exponential term

    expterm = math.exp(-gamma*rho2)

    #//
    #// u dv integral
        
    for i in range(0,9):
        term = term + (((i-1)/2)-1)*(A[i])/(T**((i-1)/2)-2)

    udvint = term*rho

    term = 0
        
    for i in range(9,17):
        term = term + (((i-10)/2))*(A[i])*(T**(1-((i-10)/2.0)))

    udvint = udvint + 0.5*term*rho2

    term = 0

    for i in range(17,21):
        term = term + (0.5+(i-18))*(A[i])*(T**(0.5-(i-18)))

    udvint = udvint + (term*rho3/3)

    term = 0

    for i in range(21,27):
        term = term + (0.5+((i-22)/4))*(A[i])*(T**(0.5-((i-22)/4)))
            
    udvint = udvint + 0.25*rho4

    udvint = udvint + (A28 + (2*A29/T2))*(rho5/5)

    term2 = 2*gamma

    udvint = udvint + (A30 + (2*A31/T) + (3*A32/T2))*((-1/term2)*(expterm - 1));

    udvint = udvint + (A33 + (2*A34/T) + (3*A35/T2))*((-rho2*expterm/term2) - ((expterm-1)/(term2*gamma)));


    #//
    #// s dv integral

    sdvint = -Rgas*math.log(rho)

    for i in range(0,9):
        term = term + (2-((i-1)/2))*(A[i])/(T**((i-1)/2.0)-1)

    sdvint = sdvint - term*rho

    term = 0

    for i in range(9,17):
        term = term + (1-((i-10)/2))*(A[i])*(T**(-((i-10)/2.0)))

    sdvint = sdvint - 0.5*term*rho2

    term = 0

    for i in range(17,21):
        term = term + (0.5-(i-18))*(A[i])*(T**(-0.5-(i-18)))

    sdvint = sdvint - (term*rho3/3)

    term = 0

    for i in range(21,27):
        term = term + (0.5-((i-22)/4))*(A[i])*(T**(-0.5-((i-22)/4)))

    sdvint = sdvint - (term*rho4/4)

    sdvint = sdvint - ((-A29/T2)*rho5/5)

    #//
    #// exponential term

    expterm = math.exp(-gamma*rho2)

    sdvint = sdvint - ((-A31/T2) - (2*A32/T3))*((expterm - 1)/(-term2))
    sdvint = sdvint - ((-A34/T2) - (2*A35/T3))*((-rho2*expterm/term2) - ((expterm-1)/(term2*gamma)))
    return udvint, sdvint

#
# PVT 6 
#
# Note that PVT 6 requires an additional Q function (Q1 or Q2) which must be specified in 
# the dictionary that is passed to the function

def pvt6(thermconst):
    """
// STANFORD PRESSURE-DENSITY-TEMPERATURE EQUATION (P-6)
//
//
//
// P =  rho R T [ 1 + rho Q + rho^2 (d Q/ d rho)|T]
//
//
    """
    rho = 1.0 / thermconst['SpVol']
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    Rgas = thermc.R_BAR/thermconst['MW']
    Q = thermconst['StanfordQ']

    if (Q == 1):
        Q = pvt6Q1(thermconst)
        dQdrho = pvt6dQ1drho(thermconst)
    elif (Q == 2):
        Q = pvt6Q2(thermconst)
        dQdrho = pvt6dQ2drho(thermconst)
        drho = 0.05*rho
        dQdrhoc = (pvt6Q2(thermconst) - Q)/drho
    else:
        print ('Invalid Q in PVT 6')

    P = rho*Rgas*T*(1 + rho*(Q + rho*dQdrho))

    return P

def pvt6Q1(thermconst):
    """
// Function StanfordPvt6Q1 will calculate the value of the
// function (Q-1) as defined in the Stanford book.  The
// function is used in conjuction with the P-v-T equation
// #6.
    """
    rho = 1.0 / thermconst['SpVol']
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    Rgas = thermc.R_BAR/thermconst['MW']

    Ta = thermconst['StanfordATa']
    tauc = thermconst['StanfordAtauc']
    Q = rhoterm = Tterm = 0

    tau = Ta/T
    taudiff = tau - tauc

    #
    # create the A array

    A = []
    load2DArray(thermconst, A, 'StanfordA', 1, 10, 1,7)


    #/////////////////////////////////////
    #// BEGIN for (i=0;i<9;i++)
    #/////////////////////////////////////
        
    for i in range(0,9):
        rhoterm = rho**i

        #/////////////////////////////////////
        #// BEGIN for (j=0;j<6;j++)
        #/////////////////////////////////////
                
        for j in range(0,6):
            Tterm = taudiff**j
            Q = Q + (A[i][j]*rhoterm*Tterm)
                
        #/////////////////////////////////////
        #// END for (j=0;j<6;j++)
        #/////////////////////////////////////

        
    #/////////////////////////////////////
    #// END for (i=0;i<9;i++)
    #/////////////////////////////////////
        

    return Q

def pvt6dQ1drho(thermconst):
    rho = 1.0 / thermconst['SpVol']
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    Rgas = thermc.R_BAR/thermconst['MW']

    Ta = thermconst['StanfordATa']
    tauc = thermconst['StanfordAtauc']
    dQdrho = rhoterm = Tterm = 0

    tau = Ta/T
    taudiff = tau - tauc

    #
    # create the A array

    A = []
    load2DArray(thermconst, A, 'StanfordA', 1, 10, 1,7)


    #/////////////////////////////////////
    #// BEGIN for (i=0;i<9;i++)
    #/////////////////////////////////////
        
    for i in range(0,9):
        rhoterm = rho**(i-1)

        #/////////////////////////////////////
        #// BEGIN for (j=0;j<6;j++)
        #/////////////////////////////////////
                
        for j in range(0,6):
            Tterm = taudiff**j
            dQdrho = dQdrho + i*(A[i][j]*rhoterm*Tterm)
                
        #/////////////////////////////////////
        #// END for (j=0;j<6;j++)
        #/////////////////////////////////////

        
    #/////////////////////////////////////
    #// END for (i=0;i<9;i++)
    #/////////////////////////////////////

    return dQdrho

def pvt6Q2(thermconst):
    rho = 1.0 / thermconst['SpVol']
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    Rgas = thermc.R_BAR/thermconst['MW']

    Ta = thermconst['StanfordATa']
    tauc = thermconst['StanfordAtauc']
    E = thermconst['StanfordAE']
    Q = rhoterm = Tterm = 0

    tau = Ta/T
    taucdiff = tau - tauc
    expterm = math.exp(-E*rho)

    #
    # create the A, tau, and rho arrays

    A = []
    load2DArray(thermconst, A, 'StanfordA', 1, 11, 1,8)

    tauA = []
    rhoA = []

    loadArray(thermconst,tauA,'StanfordTAUA',1,8)
    loadArray(thermconst,rhoA,'StanfordRHOA',1,8)

    #/////////////////////////////////////
    #// BEGIN for (j=0;j<7;j++)
    #/////////////////////////////////////
        
    for j in range (0,7):
        tauaj = tauA[j]
        tauajdiff = tau - tauaj

        Tterm = (tauajdiff**(j-1))

        #/////////////////////////////////////
        #// BEGIN for (i=0;i<8;i++)
        #/////////////////////////////////////

        rhoterm = 0

        rhoaj = rhoA[j]
        rhodiff = rho - rhoaj

        for i in range(0,8):
            k = i + 1
            l = j + 1
            Aij = A[i][j]
            rhoterm = rhoterm + Aij*(rhodiff**i)
                
        #/////////////////////////////////////
        #// END for (i=0;i<8;i++)
        #/////////////////////////////////////

        #/////////////////////////////////////
        #// BEGIN for (i=8;i<10;i++)
        #/////////////////////////////////////

        for i in range(8,10):
            k = i + 1
            l = j + 1
            Aij = A[i][j]
            rhoterm = rhoterm + expterm*Aij*(rho**(i-8))
                
        #/////////////////////////////////////
        #// END for (i=0;i<8;i++)
        #/////////////////////////////////////

        Q = Q + taucdiff*Tterm*rhoterm
        
    #/////////////////////////////////////
    #// END for (j=0;j<7;j++)
    #/////////////////////////////////////

    return Q

def pvt6dQ2drho(thermconst):
    """

    """
    rho = 1.0 / thermconst['SpVol']
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    Rgas = thermc.R_BAR/thermconst['MW']

    Ta = thermconst['StanfordATa']
    tauc = thermconst['StanfordAtauc']
    E = thermconst['StanfordAE']
    dQdrho = rhoterm = Tterm = 0

    tau = Ta/T
    taucdiff = tau - tauc
    expterm = math.exp(-E*rho)

    #
    # create the A, tau, and rho arrays

    A = []
    load2DArray(thermconst, A, 'StanfordA', 1, 11, 1,8)

    tauA = []
    rhoA = []

    loadArray(thermconst,tauA,'StanfordTAUA',1,8)
    loadArray(thermconst,rhoA,'StanfordRHOA',1,8)
        
    #/////////////////////////////////////
    #// BEGIN for (j=0;j<7;j++)
    #/////////////////////////////////////
        
    for j in range(0,7):
        tauaj = tauA[j]
        tauajdiff = tau - tauaj
                
        Tterm = (tauajdiff**(j-1))

        #/////////////////////////////////////
        #// BEGIN for (i=0;i<8;i++)
        #/////////////////////////////////////

        rhoterm = 0

        rhoaj = rhoA[j]
        rhodiff = rho - rhoaj

        for i in range(0,8):
            
            k = i + 1
            l = j + 1

            Aij = A[i][j]
            rhoterm = rhoterm + Aij*i*(rhodiff**(i-1))
                
        #/////////////////////////////////////
        #// END for (i=0;i<8;i++)
        #/////////////////////////////////////


        #/////////////////////////////////////
        #// BEGIN for (i=8;i<10;i++)
        #/////////////////////////////////////

        for i in range(8,10):
            k = i + 1
            l = j + 1

            Aij = A[i][j]
            rhoterm = rhoterm + expterm*Aij*(rho**(i-9))*((i-8) - E*rho)
                
        #/////////////////////////////////////
        #// END for (i=0;i<8;i++)
        #/////////////////////////////////////

        dQdrho = dQdrho + Tterm*rhoterm

    
        
    #/////////////////////////////////////
    #// END for (j=0;j<7;j++)
    #/////////////////////////////////////

    dQdrho = taucdiff*dQdrho
        
    return dQdrho

def dvint6(thermconst):
    """

    """
    Q = thermconst['StanfordQ']
    if (Q == 1):
        udvint, sdvint = dvint6Q1(thermconst)
    elif (Q == 2):
        udvint, sdvint = dvint6Q2(thermconst)
    else:
        print ('Invalid Q in PVT 6')

    return udvint, sdvint

def dvint6Q1(thermconst):
    """

    """
    rho = 1.0 / thermconst['SpVol']
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    Rgas = thermc.R_BAR/thermconst['MW']

    Ta = thermconst['StanfordATa']
    tauc = thermconst['StanfordAtauc']
    udvint = sdvint = rhoterm = Tterm = 0

    tau = Ta/T
    taucdiff = tau - tauc

    #
    # create the A array

    A = []
    load2DArray(thermconst, A, 'StanfordA', 1, 10, 1,7)

    #/////////////////////////////////////
    #// BEGIN for (i=0;i<9;i++)
    #/////////////////////////////////////
        
    for i in range(0,9):
        rhoterm = (rho**(i+1))

        #/////////////////////////////////////
        #// BEGIN for (j=0;j<6;j++)
        #/////////////////////////////////////
                
        for j in range(0,6):
            Tterm = (taucdiff**(j-1))
            udvint = udvint + j*rhoterm*Tterm*(A[i][j])
                
        #/////////////////////////////////////
        #// END for (j=0;j<6;j++)
        #/////////////////////////////////////

    #/////////////////////////////////////
    #// END for (i=0;i<9;i++)
    #/////////////////////////////////////

    udvint = Rgas*Ta*(udvint)
        
    #//
    #// calculate the sdv integral


    #/////////////////////////////////////
    #// BEGIN for (i=0;i<9;i++)
    #/////////////////////////////////////
        
    for i in range(0,9):
        rhoterm = (rho**(i+1))

        #/////////////////////////////////////
        #// BEGIN for (j=0;j<6;j++)
        #/////////////////////////////////////
                
        for j in range(0,6):
            Tterm = (taucdiff**(j-1))
            sdvint = sdvint + rhoterm*Tterm*(A[i][j])*(taucdiff-(j*Ta/T))
                
        #/////////////////////////////////////
        #// END for (j=0;j<6;j++)
        #/////////////////////////////////////

       
    #/////////////////////////////////////
    #// END for (i=0;i<9;i++)
    #/////////////////////////////////////

    sdvint = -Rgas*(sdvint) - Rgas*math.log(rho)

    return udvint, sdvint

def pvt6dQ2dT(thermconst):
    """

    """
    rho = 1.0 / thermconst['SpVol']
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    Rgas = thermc.R_BAR/thermconst['MW']

    Ta = thermconst['StanfordATa']
    tauc = thermconst['StanfordAtauc']
    E = thermconst['StanfordAE']
    dQdT = rhoterm = Tterm = 0

    tau = Ta/T
    taucdiff = tau - tauc
    expterm = math.exp(-E*rho)

    #
    # create the A, tau, and rho arrays

    A = []
    load2DArray(thermconst, A, 'StanfordA', 1, 11, 1,8)

    tauA = []
    rhoA = []

    loadArray(thermconst,tauA,'StanfordTAUA',1,8)
    loadArray(thermconst,rhoA,'StanfordRHOA',1,8)

    #/////////////////////////////////////
    #// BEGIN for (j=0;j<7;j++)
    #/////////////////////////////////////
        
    for j in range (0,7):
        tauaj = tauA[j]
        tauajdiff = tau - tauaj

        Tterm = (tauajdiff**(j-1))
        Tterm2 = (tauajdiff**(j-2))


        #/////////////////////////////////////
        #// BEGIN for (i=0;i<8;i++)
        #/////////////////////////////////////

        rhoterm = 0

        rhoaj = rhoA[j]
        rhodiff = rho - rhoaj

        for i in range(0,8):
            Aij = A[i][j]
            rhoterm = rhoterm + Aij*(rhodiff**i)
                
        #/////////////////////////////////////
        #// END for (i=0;i<8;i++)
        #/////////////////////////////////////

        #/////////////////////////////////////
        #// BEGIN for (i=8;i<10;i++)
        #/////////////////////////////////////

        for i in range(8,10):
            Aij = A[i][j]
            rhoterm = rhoterm + expterm*Aij*(rho**(i-8))
                
        #/////////////////////////////////////
        #// END for (i=0;i<8;i++)
        #/////////////////////////////////////

        dQdT = dQdT - rhoterm*(Ta/(T*T))*(taucdiff*(j-1)*Tterm2 + Tterm)

    #/////////////////////////////////////
    #// END for (j=0;j<7;j++)
    #/////////////////////////////////////

    return dQdT

def dvint6Q2(thermconst):
    """

    """
    rho = 1.0 / thermconst['SpVol']
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    Rgas = thermc.R_BAR/thermconst['MW']

    Q = pvt6Q2(thermconst)
    dQdT = pvt6dQ2dT(thermconst)

    udvint = -Rgas*rho*T*T*dQdT
    sdvint = -Rgas*rho*(Q + T*dQdT) - Rgas*math.log(rho)

    return udvint, sdvint


#
# pvt - 7
#
def pvt7(thermconst):
    """
//
// STANFORD PRESSURE-DENSITY-TEMPERATURE EQUATION (P-7)
//
//
//
// P =  rho R T + rho^2 B(T) + rho^3 C(T) + rho^4  D(T)
//
//              + rho^5 E(T)
//
//      B(T) = -T^2 exp(A1 + (A2 / T))
//
//      C(T) = T exp(A3 + (A4/T) + (A5/T^2))
//
//      D(T) = -T exp(A6 + (A7/T) + (A8/T^2))
//
//      E(T) = A9 T
//
    """

    rho = 1.0 / thermconst['SpVol']
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    Rgas = thermc.R_BAR/thermconst['MW']
    #//
    #// calculate powers of rho

    rho2 = rho*rho
    rho3 = rho2*rho
    rho4 = rho3*rho
    rho5 = rho4*rho
    rho6 = rho5*rho

    B = pvt7B(thermconst)
    C = pvt7C(thermconst)
    D = pvt7D(thermconst)
    E = pvt7E(thermconst)

    P = rho*Rgas*T + rho2*B + rho3*C + rho4*D + rho5*E

    return P

def pvt7B(thermconst):
    """
    Calculate the B function for P-7

    B(T) = -T^2 exp(A1 + (A2 / T))
    """
    A1 = thermconst['StanfordA1']
    A2 = thermconst['StanfordA2']
    T = thermconst['Temperature']

    expterm = math.exp(A1 + (A2/T))
    B = -T*T*expterm
    return B

def pvt7dBdT(thermconst):
    A1 = thermconst['StanfordA1']
    A2 = thermconst['StanfordA2']
    T = thermconst['Temperature']
    
    expterm = math.exp(A1 + (A2/T))
    dBdT = -2*T*expterm + A2*expterm

    return dBdT

def  pvt7C(thermconst):
    """
    C function for P-7

    C(T) = T exp(A3 + (A4/T) + (A5/T^2))
    """
    A3 = thermconst['StanfordA3']
    A4 = thermconst['StanfordA4']
    A5 = thermconst['StanfordA5']
    T = thermconst['Temperature']

    expterm = math.exp(A3 + (A4/T) + (A5/(T*T)))
    C = T*expterm
    return C

def pvt7dCdT(thermconst):
    A3 = thermconst['StanfordA3']
    A4 = thermconst['StanfordA4']
    A5 = thermconst['StanfordA5']
    T = thermconst['Temperature']
    expterm = math.exp(A3 + (A4/T) + (A5/(T*T)))

    dCdT = expterm + T*expterm*((-A4/(T*T)) - (2*A5/(T*T*T)))
    return dCdT

def pvt7D(thermconst):
    """

    D(T) = -T exp(A6 + (A7/T) + (A8/T^2))

    """
    A6 = thermconst['StanfordA6']
    A7 = thermconst['StanfordA7']
    A8 = thermconst['StanfordA8']
    T = thermconst['Temperature']
    expterm = math.exp(A6 + (A7/T) + (A8/(T*T)))

    D = -T*expterm

    return D

def pvt7dDdT(thermconst):
    A6 = thermconst['StanfordA6']
    A7 = thermconst['StanfordA7']
    A8 = thermconst['StanfordA8']
    T = thermconst['Temperature']
    expterm = math.exp(A6 + (A7/T) + (A8/(T*T)))

    dDdT = -expterm - T*expterm*((-A7/(T*T)) - (2*A8/(T*T*T)))
    return dDdT

def pvt7E(thermconst):
    """
     
    E(T) = A9 T

    """
    A9 = thermconst['StanfordA9']
    T = thermconst['Temperature']

    E = A9 * T

    return E

def pvt7dEdT(thermconst):

    return thermconst['StanfordA9']


def dvint7(thermconst):
    """

    dv integrals for P-7

    """
    rho = 1.0 / thermconst['SpVol']
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    Rgas = thermc.R_BAR/thermconst['MW']
    #//
    #// calculate powers of rho

    rho2 = rho*rho
    rho3 = rho2*rho
    rho4 = rho3*rho
    rho5 = rho4*rho
    rho6 = rho5*rho

    B = pvt7B(thermconst)
    C = pvt7C(thermconst)
    D = pvt7D(thermconst)
    E = pvt7E(thermconst)
    dBdT = pvt7dBdT(thermconst)
    dCdT = pvt7dCdT(thermconst)
    dDdT = pvt7dDdT(thermconst)
    dEdT = pvt7dEdT(thermconst)

    #//
    #// calculate u dv integral

    udvint = (B - T*dBdT)*rho + ((C - T*dCdT)*rho2/2) + ((D - T*dDdT)*rho3/3) + ((E - T*dEdT)*rho4/4)

    #//
    #// calculate s dv integral

    sdvint = -dBdT*rho - (dCdT*rho2/2) - (dDdT*rho3/3) - (dEdT*rho4/4)

    sdvint = sdvint - Rgas*math.log(rho)

    return udvint, sdvint

#
# P-8
#

def pvt8(thermconst):
    """
//
// P =  R T/(v-b) + [a / (v (v + b) T^(1/2))]
//

    Note that this is essentially the Redlich Kwong eqn
    """
    v = thermconst['SpVol']
    T = thermconst['Temperature']
    Rgas = thermc.R_BAR/thermconst['MW']
    a = thermconst['StanfordAa']
    b = thermconst['StanfordAb']

    P = (Rgas*T/(v-b)) - (a/(v*(v+b)*math.sqrt(T)))

    return P

def dvint8(thermconst):
    v = thermconst['SpVol']
    T = thermconst['Temperature']
    Rgas = thermc.R_BAR/thermconst['MW']
    a = thermconst['StanfordAa']
    b = thermconst['StanfordAb']

    term = math.log(v/(v+b))
    term2 = math.sqrt(T)
        
    #//
    #// calculate the u dv integral

    udvint = 3*a*term/(2*b*term2)

    #//
    #// calculate the s dv integral

    sdvint = a*term/(2*b*T*term2)

    sdvint = sdvint + Rgas*math.log(v-b)

    return udvint, sdvint


def pvt9(thermconst):
    """
//
// P = rho R T - rho^2 T exp(A_1 + (A_2/T) + (A_3 ln(T)))

    """
    A1 = thermconst['StanfordA1']
    A2 = thermconst['StanfordA2']
    A3 = thermconst['StanfordA3']

    rho = 1.0 / thermconst['SpVol']
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    Rgas = thermc.R_BAR/thermconst['MW']

    rho2 = rho**2

    P = rho*Rgas*T - rho2*T*math.exp(A1 + (A2/T) + (A3*math.log(T)))

    return P

def dvint9(thermconst):
    A1 = thermconst['StanfordA1']
    A2 = thermconst['StanfordA2']
    A3 = thermconst['StanfordA3']

    rho = 1.0 / thermconst['SpVol']
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    Rgas = thermc.R_BAR/thermconst['MW']

    rho2 = rho**2

    term = math.exp(A1 + (A2/T) + (A3*math.log(T)))

    udvint = term*(-A2 + (A3*T))*rho

    sdvint = term*(1 + (A2/T) - A3)*rho - Rgas*math.log(rho)

    return udvint, sdvint

def pvt10(thermconst):
    """
// P =  rho R T + rho^2 SUM[i=1,4]{A_i T^(2-i)} + rho^3 SUM[i=5,9]{A_i T^(7-i)}
//              
//              + rho^4 (A_10 T + A_11) + rho^5 (A_12 / T) + rho^6 A_13
//
//              + rho^3 (A_14 / T^2 + A_15 / T^3) exp(-gamma rho^2)
//
//              + rho^5 (A_16 / T^2 + A_17 / T^3) exp(-gamma rho^2)
//
    """
    #
    # create the A array

    A = []
    #for i in range(1, 20):
    #    key = 'StanfordA' + str(i)
    #    A.append(thermconst[key])
    loadArray(thermconst, A, 'StanfordA', 1, 18)
    # print (Atest)

    
    
    A10 = thermconst['StanfordA10']
    A11 = thermconst['StanfordA11']
    A12 = thermconst['StanfordA12']
    A13 = thermconst['StanfordA13']
    A14 = thermconst['StanfordA14']
    A15 = thermconst['StanfordA15']
    A16 = thermconst['StanfordA16']
    A17 = thermconst['StanfordA17']
    gamma = thermconst['StanfordAgamma']
        
    # print (A)

    rho = 1.0 / thermconst['SpVol']
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    Rgas = thermc.R_BAR/thermconst['MW']

    #//
    #// calcualte powers of rho

    rho2 = rho*rho
    rho3 = rho2*rho
    rho4 = rho3*rho
    rho5 = rho4*rho
    rho6 = rho5*rho


    #//
    #// powers of T

    T2 = T*T
    T3 = T2*T
    #//
    #// 1st term in equation

    P = rho*Rgas*T

    #//
    #// 2nd term in equation

    term = 0

    for i in range(0,4):
        term = term + (A[i])*(T**(2.0-(i+1)))

    term = term*rho2

    P = P + term

    #//
    #// 3rd term in equation

    term = 0

    for i in range(4,9):
        term = term + (A[i])*(T**(7.0-(i+1)))

    term = term*rho3

    P = P + term

    #//
    #// 4th term in equation

    P = P + rho4*(A10*T + A11)

    #//
    #// 5th term in equation

    P = P + (rho5*A12/T)

    #//
    #// 6th term in equation

    P = P + rho6*A13

    expterm = math.exp(-gamma*rho2)

    #//
    #// 7th term in equation

    P = P + rho3*expterm*((A14/T2) + (A15/T3))


    #//
    #// 8th term in equation

    P = P + rho5*expterm*((A16/T2) + (A17/T3))
        
    return P

def dvint10(thermconst):
    #
    # create the A array

    A = []
    #for i in range(1, 20):
    #    key = 'StanfordA' + str(i)
    #    A.append(thermconst[key])
    loadArray(thermconst, A, 'StanfordA', 1, 18)
    # print (Atest)

    
    
    A10 = thermconst['StanfordA10']
    A11 = thermconst['StanfordA11']
    A12 = thermconst['StanfordA12']
    A13 = thermconst['StanfordA13']
    A14 = thermconst['StanfordA14']
    A15 = thermconst['StanfordA15']
    A16 = thermconst['StanfordA16']
    A17 = thermconst['StanfordA17']
    gamma = thermconst['StanfordAgamma']
        
    # print (A)

    rho = 1.0 / thermconst['SpVol']
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    Rgas = thermc.R_BAR/thermconst['MW']

    #//
    #// calcualte powers of rho

    rho2 = rho*rho
    rho3 = rho2*rho
    rho4 = rho3*rho
    rho5 = rho4*rho
    rho6 = rho5*rho


    #//
    #// powers of T

    T2 = T*T
    T3 = T2*T
    T4 = T3*T

    #//
    #// calculate u dv int


    #//
    #// 2nd term in equation

    term = 0

    for i in range(0,4):
        term = term + (i*(A[i])*(T**(1-i)))

    term = term*rho

    udvint = term

    #//
    #// 3rd term in equation

    term = 0

    for i in range(4,9):
        term = term + ((i-5)*(A[i])*(T**(6-i)))

    udvint = udvint + 0.5*term*rho2
    udvint = udvint + (A11*rho3/3) + (A12*rho4/(2*T)) + (A13*rho5/5)

    expterm = math.exp(-gamma*rho2)
    term = (expterm - 1)/(2*gamma)
    term2 = ((3*A14/T2) + (4*A15/T3))
        
    udvint = udvint - term*term2

    term2 = ((3*A16/T2) + (4*A17/T3))
    term = -(rho2*expterm/(2*gamma)) - ((expterm - 1)/(2*gamma*gamma))

    udvint = udvint + term2*term

    #//
    #// calculate s dv integral

    term = 0

    for i in range(0,4):
        term = term + ((1-i)*(A[i])/(T**i))

    sdvint = -term*rho

    term = 0

    for i in range(4,9):
        term = term + ((6-i)*(A[i])*(T**(5-i)))

    sdvint = sdvint - 0.5*term*rho2
    sdvint = sdvint - (A10*rho3/3) + (A12*rho4/(4*T2))

    expterm = math.exp(-gamma*rho2)
    term = (expterm - 1)/(2*gamma)
    term2 = ((-2*A14/T3) - (3*A15/T4))
        
    sdvint = sdvint + term*term2

    term2 = ((-2*A16/T3) - (3*A17/T4))

    term = -(rho2*expterm/(2*gamma)) - ((expterm - 1)/(2*gamma*gamma))


    sdvint = sdvint + term2*term - Rgas*math.log(rho)

    return udvint, sdvint

def pvt11(thermconst):
    """
// STANFORD PRESSURE-DENSITY-TEMPERATURE EQUATION (P-11)
//
//
//
// P =  rho R T + rho^2 B1(T) + rho^3 B2(T) + rho^4 B3(T)
//
//      B1(T) = (A_1 T + A_2 T^2 + A_3 T^3 + A_4 T^4) exp(alpha/T)
//
//      B2(T) = (A_5 T + A_6 T^2 + A_7 T^3 ) exp(alpha/T)
//
//      B3(T) = (A_8 T + A_9 T^2 ) exp(alpha/T)

    """

    rho = 1.0 / thermconst['SpVol']
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    Rgas = thermc.R_BAR/thermconst['MW']
    #//
    #// calcualte powers of rho

    rho2 = rho*rho
    rho3 = rho2*rho
    rho4 = rho3*rho
    rho5 = rho4*rho
    rho6 = rho5*rho

    B1 = pvt11B1(thermconst)
    B2 = pvt11B2(thermconst)
    B3 = pvt11B3(thermconst)
    P = rho*Rgas*T + rho2*B1 + rho3*B2 + rho4*B3

    return P

def pvt11B1(thermconst):

    A1 = thermconst['StanfordA1']
    A2 = thermconst['StanfordA2']
    A3 = thermconst['StanfordA3']
    A4 = thermconst['StanfordA4']
    alpha = thermconst['StanfordAalpha']

    T = thermconst['Temperature']
    
    expterm = math.exp(alpha/T)

    B1 = (T*(A1 + T*(A2 + T*(A3 + T*A4))))*expterm

    return B1

def pvt11dB1dt(thermconst):

    A1 = thermconst['StanfordA1']
    A2 = thermconst['StanfordA2']
    A3 = thermconst['StanfordA3']
    A4 = thermconst['StanfordA4']
    alpha = thermconst['StanfordAalpha']

    T = thermconst['Temperature']
    
    expterm = math.exp(alpha/T)

    dB1dT = ((A1 + T*(2*A2 + T*(3*A3 + 4*A4*T))))*expterm - (alpha*(T*(A1 + T*(A2 + T*(A3 + T*A4))))*expterm/(T*T))

    return dB1dT


def pvt11B2(thermconst):

    A5 = thermconst['StanfordA5']
    A6 = thermconst['StanfordA6']
    A7 = thermconst['StanfordA7']
    alpha = thermconst['StanfordAalpha']

    T = thermconst['Temperature']
    
    expterm = math.exp(alpha/T)

    B2 = (T*(A5 + T*(A6 + A7*T)))*expterm

    return B2

def pvt11dB2dt(thermconst):
    A5 = thermconst['StanfordA5']
    A6 = thermconst['StanfordA6']
    A7 = thermconst['StanfordA7']
    alpha = thermconst['StanfordAalpha']

    T = thermconst['Temperature']
    
    expterm = math.exp(alpha/T)

    dB2dT = (A5 + T*(2*A6 + 3*A7*T))*expterm - (alpha*(T*(A5 + T*(A6 + A7*T))*expterm/(T*T)))

    return dB2dT

def pvt11B3(thermconst):

    A8 = thermconst['StanfordA8']
    A9 = thermconst['StanfordA9']
    alpha = thermconst['StanfordAalpha']

    T = thermconst['Temperature']
    
    expterm = math.exp(alpha/T)

    B3 = (T*(A8 + A9*T))*expterm

    return B3

def pvt11dB3dt(thermconst):

    A8 = thermconst['StanfordA8']
    A9 = thermconst['StanfordA9']
    alpha = thermconst['StanfordAalpha']

    T = thermconst['Temperature']
    
    expterm = math.exp(alpha/T)

    dB3dT = (A8 + 2*A9*T)*expterm - (alpha*(T*(A8 + A9*T))*expterm/(T*T))

    return dB3dT

def dvint11(thermconst):

    rho = 1.0 / thermconst['SpVol']
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    Rgas = thermc.R_BAR/thermconst['MW']
    #//
    #// calcualte powers of rho

    rho2 = rho*rho
    rho3 = rho2*rho
    rho4 = rho3*rho
    rho5 = rho4*rho
    rho6 = rho5*rho

    B1 = pvt11B1(thermconst)
    B2 = pvt11B2(thermconst)
    B3 = pvt11B3(thermconst)
    dB1dT = pvt11dB1dt(thermconst)
    dB2dT = pvt11dB2dt(thermconst)
    dB3dT = pvt11dB3dt(thermconst)

    #//
    #// calculate u dv integral

    udvint = (B1 - T*dB1dT)*rho + ((B2 - T*dB2dT)*rho2/2) + ((B3 - T*dB3dT)*rho3/3)

    #//
    #// calculate s dv integral

    sdvint = -dB1dT*rho - (dB2dT*rho2/2) - (dB3dT*rho3/3)

    sdvint = sdvint - Rgas*math.log(rho)

    return udvint, sdvint


#
#  SATURATION EQUATIONS
# 

def psat1(thermconst):
    """
/////////////////////////////////////////////////
// Function StanfordS1 will calculate the 
// saturation pressure given the temperature
// as based on equation (S-1) in the Stanford
// properties book

        //
        // DECLARATIONS
        //
        //      INTEGERS
        //
        //              i                       array index, counter
        //              offset          offset in consts array to 
        //                                      start of saturation pressure
        //                                      data
        //
        //      DOUBLE
        //
        //              F1                      equation parameter
        //              F2                      equation parameter
        //              F3                      equation parameter
        //              F4                      equation parameter
        //              F5                      equation parameter
        //              gamma           equation parameter
        //              Psat            saturation pressure                     [Pa]
        //              term            intermediate calculation
        //                                      term



    """
    T = thermconst['Temperature']
    F1 = thermconst['StanfordF1']
    F2 = thermconst['StanfordF2']
    F3 = thermconst['StanfordF3']
    F4 = thermconst['StanfordF4']
    F5 = thermconst['StanfordF5']
    gamma = thermconst['StanfordFgamma']

    Psat = F1 + (F2/T) + (F3 * math.log(T)) + (F4 * T) 
    if gamma != 0:
        Psat = Psat + (F5 * (gamma - T) * math.log(gamma - T) / T)

    return math.exp(Psat)

#
# dPsat/dT for S-1

def dpsatdt1(thermconst):
    """
    // Function dpsatdt1 will calculate the 
// derivative of the saturation pressure with 
// respect to temperature given the temperature
// as based on equation (S-1) in the Stanford
// properties book
    """
    T = thermconst['Temperature']
    F1 = thermconst['StanfordF1']
    F2 = thermconst['StanfordF2']
    F3 = thermconst['StanfordF3']
    F4 = thermconst['StanfordF4']
    F5 = thermconst['StanfordF5']
    gamma = thermconst['StanfordFgamma']

    Psat = F1 + (F2/T) + (F3 * math.log(T)) + (F4 * T) 
    if gamma != 0:
        Psat = Psat + (F5 * (gamma - T) * math.log(gamma - T) / T)

    Psat = math.exp(Psat)

    dPdT = (-F2/(T*T)) + (F3/T) + F4
    logterm = math.log(gamma -T)
    term = gamma -T

    if (F5 != 0):
        dPdT = dPdT - ((F5 * gamma * logterm + F5 * T)/(T*T))
    dPdT = Psat*dPdT
    return dPdT


def psat2(thermconst):
    """
 Function StanfordS2 will calculate the 
 saturation pressure given the temperature
 as based on equation (S-2) in the Stanford
 properties book

      INPUT:
      
             double          *consts         pointer to constants 
                                                              array
             double          T                       temperature                             [K]

      OUTPUT:

              NULL

      RETURN:

              double          Psat            saturation pressure             [Pa]


    """
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    Pc = thermconst['PCrit']
    Rgas = thermc.R_BAR/thermconst['MW']

    #
    # create the F array

    F = []
    ###for i in range(1, 9):
    ###    key = 'StanfordF' + str(i)
    ###    F.append(thermconst[key])

    loadArray(thermconst,F,'StanfordF',1,9)
        
    #print (F)

    Tp = thermconst['StanfordTp']

    term = 0
    TpTerm = (T/Tp) - 1
    for i in range(0,8):
        term = term + F[i] * (TpTerm)**i

    term = ((Tc/T) - 1) * term
    
    Psat = Pc * math.exp(term)

    return Psat


#
# dPsat/dT S-2

def dpsatdt2(thermconst):
    """

    """
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    Pc = thermconst['PCrit']
    Rgas = thermc.R_BAR/thermconst['MW']

    #
    # create the F array

    F = []

    loadArray(thermconst,F,'StanfordF',1,9)
        
    #print (F)

    Tp = thermconst['StanfordTp']

    term = 0
    TpTerm = (T/Tp) - 1
    for i in range(0,8):
        term = term + F[i] * (TpTerm**i)

    term2 = ((Tc/T) - 1) * term
    
    Psat = Pc * math.exp(term2)

    dPdT = Psat*term*(-Tc)/(T*T);

    term = F[1];

    for i in range(2,8):
        term2 = TpTerm**(i-1)
        term = term + (i*F[i]*term2)

    dPdT = dPdT + (Psat*term*((Tc/T) - 1)/Tp)

    return dPdT


def psat3(thermconst):
    """
    
 STANFORD SATURATION PRESSURE EQUATION (S-3)

 Ln(P/Pt) = F_1 X +  F_2 X^2 + F_3 X^3 + F_4 X (1 - X)^(alpha)

      X = (1 - (Tt/T)) / (1 - (Tt/Tc))
    """

    T = thermconst['Temperature']
    F1 = thermconst['StanfordF1']
    F2 = thermconst['StanfordF2']
    F3 = thermconst['StanfordF3']
    F4 = thermconst['StanfordF4']
    Pt = thermconst['StanfordPt']
    Tt = thermconst['StanfordTt']
    alpha = thermconst['StanfordFalpha']

    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    X = (1 - (Tt/T))/(1 - (Tt/Tc))
        
    Psat = X*(F1 + X*(F2 + X*F3))

    Psat = Psat + F4*X*((1-X)**alpha)
        
    Psat = Pt*math.exp(Psat)

    return Psat

#
# dPsat/dT S-3

def dpsatdt3(thermconst):
    """

    """
    T = thermconst['Temperature']
    F1 = thermconst['StanfordF1']
    F2 = thermconst['StanfordF2']
    F3 = thermconst['StanfordF3']
    F4 = thermconst['StanfordF4']
    Pt = thermconst['StanfordPt']
    Tt = thermconst['StanfordTt']
    alpha = thermconst['StanfordFalpha']

    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    X = (1 - (Tt/T))/(1 - (Tt/Tc))
        
    Psat = X*(F1 + X*(F2 + X*F3))

    if (alpha != 0):
            Psat = Psat + F4*X*((1-X)**alpha)
    else:
        Psat = Psat + F4*X
        
    Psat = Pt*math.exp(Psat)

    dXdT = Tt/(T*T*(1 - (Tt/Tc)))

    term1 = (1-X)**alpha
    term2 = (1-X)**(alpha - 1)
 
    term1 = term1*dXdT
    term2 = term2*alpha*X*dXdT

    dPdT = Psat*(F1*dXdT + X*(2*F2*dXdT + 3*F3*X*dXdT) + F4*(term1 - term2))
                        
    return dPdT

def psat4(thermconst):
    """
//
// STANFORD SATURATION PRESSURE EQUATION (S-4)
//
// Ln(P) = F_1/T +  F_2 + F_3 T + F_4  (T_c - T)^alpha
//                      + F_5 T^3 + F_6 T^4 + F_7 T^5 + F_8 T^6
//                      + F_9 ln (T)
//
//
    """

    T = thermconst['Temperature']
    F1 = thermconst['StanfordF1']
    F2 = thermconst['StanfordF2']
    F3 = thermconst['StanfordF3']
    F4 = thermconst['StanfordF4']
    F5 = thermconst['StanfordF5']
    F6 = thermconst['StanfordF6']
    F7 = thermconst['StanfordF7']
    F8 = thermconst['StanfordF8']
    F9 = thermconst['StanfordF9']
    alpha = thermconst['StanfordFalpha']

    #    
    # set powers of T

    T2 = T*T
    T3 = T2*T
    T4 = T3*T
    T5 = T4*T
    T6 = T5*T

    Psat = (F1/T) + F2 + F3*T + F4*((Tc-T)**alpha) + F5*T3 + F6*T4 \
                    + F7*T5 + F8*T6 + F9*math.log(T)

    Psat = math.exp(Psat)

    return Psat

#
# dPsatdT S-4

def dpsatdt4(thermconst):

    T = thermconst['Temperature']
    F1 = thermconst['StanfordF1']
    F2 = thermconst['StanfordF2']
    F3 = thermconst['StanfordF3']
    F4 = thermconst['StanfordF4']
    F5 = thermconst['StanfordF5']
    F6 = thermconst['StanfordF6']
    F7 = thermconst['StanfordF7']
    F8 = thermconst['StanfordF8']
    F9 = thermconst['StanfordF9']
    alpha = thermconst['StanfordFalpha']

    #    
    # set powers of T

    T2 = T*T
    T3 = T2*T
    T4 = T3*T
    T5 = T4*T
    T6 = T5*T

    Psat = (F1/T) + F2 + F3*T + F4*((Tc-T)**alpha) + F5*T3 + F6*T4 \
                    + F7*T5 + F8*T6 + F9*math.log(T)

    Psat = math.exp(Psat)

    dPdT = Psat*((-F1/T2) + F3 + 3*F5*T2 + 4*F6*T3 + 5*F7*T4 + 6*F8*T5 + (F9/T));
    
    dPdT = dPdT - (Psat*(alpha*F4*((Tc-T)**alpha-1)))
                        
    return dPdT


def psat5(thermconst):
    """
    //
    // STANFORD SATURATION PRESSURE EQUATION (S-5)
    //
    // Ln(P) = SUM{i=1,10} [Fi T^(2-i)]
    //
    //
    """

    T = thermconst['Temperature']
    #
    # create the F array

    F = []

    loadArray(thermconst,F,'StanfordF',1,11)

    Psat = 0

    for i in range(0,10):
        Psat = Psat + F[i] * (T**(1-i))

    Psat = math.exp(Psat)
    return Psat

def dpsatdt5(thermconst):
    #
    # create the F array

    T = thermconst['Temperature']
    F = []

    loadArray(thermconst,F,'StanfordF',1,11)

    Psat = 0

    for i in range(0,10):
        Psat = Psat + F[i] * (T**(1-i))

    Psat = math.exp(Psat)

    dPdT = 0
        
    for i in range(0,10):
        dPdT = dPdT + ((1-i) * F[i] / (T**(i)))       

    dPdT = Psat*dPdT
    return dPdT

def psat6(thermconst):
    """
    //
    // STANFORD SATURATION PRESSURE EQUATION (S-6)
    //
    // Ln(P/P_c) = (T_c/T - 1) SUM(i=1,8){ F_i [a(T - T_p )]^(i-1)
    //
    //
    """
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    Pc = thermconst['PCrit']
        

    Tp = thermconst['StanfordTp']
    a = thermconst['StanfordFa']

    F = []
    loadArray(thermconst,F,'StanfordF',1,11)
        

    term2 = a*(T - Tp)


    term = F[0]

    for i in range (1,8):
        term3 = term2**i
        term = term + F[i]*term3

    term = ((Tc/T) - 1)*term

    Psat = Pc*math.exp(term)

    return Psat

def dpsatdt6(thermconst):
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    Pc = thermconst['PCrit']
        

    Tp = thermconst['StanfordTp']
    a = thermconst['StanfordFa']

    F = []
    loadArray(thermconst,F,'StanfordF',1,11)
        

    term2 = a*(T - Tp)


    term = F[0]

    for i in range (1,8):
        term3 = term2**i
        term = term + F[i]*term3

    Psat = Pc*math.exp(((Tc/T) - 1)*term)

    dPdT = Psat*term*(-Tc)/(T*T)

    term = a*F[1]

    for i in range (2,8):
        term4 = term2**(i-1)
        term = term + (i*a*(F[i])*term4)

    dPdT = dPdT + (Psat*term*((Tc/T) - 1))

    return dPdT

def psat7(thermconst):
    """
    //
    // STANFORD SATURATION PRESSURE EQUATION (S-7)
    //
    // Ln(P) = F_1 +  F_2/T + F_3 ln (T)
    //
    //
    """

    T = thermconst['Temperature']
    F1 = thermconst['StanfordF1']
    F2 = thermconst['StanfordF2']
    F3 = thermconst['StanfordF3']

    Psat = F1 + (F2/T) + F3*math.log(T)

    Psat = math.exp(Psat)

    return Psat

def dpsatdt7(thermconst):
    T = thermconst['Temperature']
    F1 = thermconst['StanfordF1']
    F2 = thermconst['StanfordF2']
    F3 = thermconst['StanfordF3']

    Psat = F1 + (F2/T) + F3*math.log(T)

    Psat = math.exp(Psat)
    T2 = T*T
    dPdT = Psat*((-F2/T2) + (F3/T))
    return dPdT

def psat8(thermconst):
    """
// Ln(P) = F_1 +  F_2/T + F_3 ln (T)
//                      + F_4 T + F_5 T^2 + F_6 T^3

    """

    T = thermconst['Temperature']
    F1 = thermconst['StanfordF1']
    F2 = thermconst['StanfordF2']
    F3 = thermconst['StanfordF3']
    F4 = thermconst['StanfordF4']
    F5 = thermconst['StanfordF5']
    F6 = thermconst['StanfordF6']

    Psat = F1 + (F2/T) + F3*math.log(T) + T*(F4 + T*(F5 + F6*T))

    Psat = math.exp(Psat)

    return Psat

def dpsatdt8(thermconst):
    T = thermconst['Temperature']
    F1 = thermconst['StanfordF1']
    F2 = thermconst['StanfordF2']
    F3 = thermconst['StanfordF3']
    F4 = thermconst['StanfordF4']
    F5 = thermconst['StanfordF5']
    F6 = thermconst['StanfordF6']

    Psat = F1 + (F2/T) + F3*math.log(T) + T*(F4 + T*(F5 + F6*T))

    Psat = math.exp(Psat)

    dPdT = (-F2/(T*T)) + (F3/T) + F4 + T*(2*F5 + 3*F6*T)

    dPdT = Psat*dPdT

    return dPdT

def psat9(thermconst):
    """
//
// Ln(P) = F_1 +  F_2/T + F_3 ln (T) + (F4 / T^2)
//
//

    """

    T = thermconst['Temperature']
    F1 = thermconst['StanfordF1']
    F2 = thermconst['StanfordF2']
    F3 = thermconst['StanfordF3']
    F4 = thermconst['StanfordF4']
    Psat = F1 + (F2/T) + F3*math.log(T) + (F4/(T*T))

    Psat = math.exp(Psat)

    return Psat

def dpsatdt9(thermconst):

    T = thermconst['Temperature']
    F1 = thermconst['StanfordF1']
    F2 = thermconst['StanfordF2']
    F3 = thermconst['StanfordF3']
    F4 = thermconst['StanfordF4']

    T2 = T*T
    Psat = F1 + (F2/T) + F3*math.log(T) + (F4/(T*T))

    Psat = math.exp(Psat)
    dPdT = Psat*((-F2/T2) + (F3/T) - (2*F4/(T*T*T)))
    return dPdT

def psat10(thermconst):
    """
//
// Ln(P) = F_1 + F_2/T + F_3 T + F_4 T^2
//
//

    """

    T = thermconst['Temperature']
    F1 = thermconst['StanfordF1']
    F2 = thermconst['StanfordF2']
    F3 = thermconst['StanfordF3']
    F4 = thermconst['StanfordF4']

    Psat = F1 + (F2/T) + F3*T + F4*T*T
    Psat = math.exp(Psat)
    return Psat


def dpsatdt10(thermconst):
    T = thermconst['Temperature']
    F1 = thermconst['StanfordF1']
    F2 = thermconst['StanfordF2']
    F3 = thermconst['StanfordF3']
    F4 = thermconst['StanfordF4']

    Psat = F1 + (F2/T) + F3*T + F4*T*T
    Psat = math.exp(Psat)
    dPdT = Psat*((-F2/(T*T)) + F3 + 2*F4*T)
    return dPdT

#
#  CV EQUATIONS
# 

#
# C-0

def cv0(thermconst):
    """
          Cv = G1
    """
    T = thermconst['Temperature']
    T0 = thermconst['StanfordT0']
    G1 = thermconst['StanfordG1']

    cv = G1
    cvint = G1 * (T- T0)
    cvtint = G1 * math.log(T/T0)

    return cv, cvint, cvtint

def cv1(thermconst):
    """
// STANFORD SPECIFIC HEAT EQUATION (C-1)
//
//      Cv = Sum(i=1,4) {G_i T^(i-1)} + G_5/T^2


    """
    T = thermconst['Temperature']
    T0 = thermconst['StanfordT0']

    #
    # create the G array

    G = []
    #for i in range(1, N+1):
    #    key = 'StanfordG' + str(i)
    #    G.append(thermconst[key])

    loadArray(thermconst,G,'StanfordG',1,6)

    print (G)
    #
    # calculate cv

    cv = 0

    for i in range(0,4):
        cv = cv + G[i] * T**i

    cv = cv + G[4]/(T*T)

    #
    # calculate INT[T0,T] {Cv dT}

    cvint = 0

    for i in range(0,5):
        cvint = cvint + (G[i] / (i+1)) * ((T**(i+1)) - (T0**(i+1)))


    cvint = cvint - (G[4])*((1/T) - (1/T0))

    #
    # calculate INT[T0,T] {Cv/T dT}
    cvtint = G[0] * math.log(T/T0)

    for i in range(1,4):
        cvtint = cvtint + (G[i]/i) * (T**(i) - T0**(i))

    cvtint = cvtint - ((G[4]/(2*T*T)) - (G[4]/(2*T0*T0)))

    return cv, cvint, cvtint 

#
# C-2

def cv2(thermconst):
    """
#
# STANFORD SPECIFIC HEAT EQUATION (C-2)
#
#      Cv = Sum(i=1,N) {G_i T^(i-1)}
    """
    T = thermconst['Temperature']
    T0 = thermconst['StanfordT0']
    N = thermconst['StanfordcvN']

    #
    # create the G array

    G = []
    #for i in range(1, N+1):
    #    key = 'StanfordG' + str(i)
    #    G.append(thermconst[key])

    loadArray(thermconst,G,'StanfordG',1,N+1)

    print (G)
    #
    # calculate cv

    cv = 0

    for i in range(0,N):
        cv = cv + G[i] * T**i

    #
    # calculate INT[T0,T] {Cv dT}

    cvint = 0

    for i in range(0,N):
        cvint = cvint + (G[i] / (i+1)) * ((T**(i+1)) - (T0**(i+1)))

    #
    # calculate INT[T0,T] {Cv/T dT}
    cvtint = G[0] * math.log(T/T0)

    for i in range(1,N):
        cvtint = cvtint + (G[i]/i) * (T**(i) - T0**(i))

    return cv, cvint, cvtint 

#
# C-3
#

def cv3(thermconst):
    """
//
// STANFORD SPECIFIC HEAT EQUATION (C-3)
//
//      Cv = G_1 + G_2 T^(1/3) + G_3 T^(2/3) + G_4 T + G_5 exp(beta/T)  [(beta/T)/(exp(beta/T) - 1)]^2
//


/////////////////////////////////////////////////
// Function StanfordC3 will calculate the 
// specific heat, and two integrals for the 
// specific heat
//

    """

    
    T = thermconst['Temperature']
    G1 = thermconst['StanfordG1']
    G2 = thermconst['StanfordG2']
    G3 = thermconst['StanfordG3']
    G4 = thermconst['StanfordG4']
    G5 = thermconst['StanfordG5']
    T0 = thermconst['StanfordT0']

    beta = thermconst['StanfordGbeta']

    #//
    #// Set T0



    #//
    #// powers of T

    T2 = T*T
    T13 =T**(1/3)
    T23 = T**(2/3)
    T43 = T**(4/3)
    T53 = T**(5/3)
    T02 = T0*T0
    T013 = T0**(1/3)
    T023 = T0**(2/3)
    T043 = T0**(4/3)
    T053 = T0**(5/3)

    expterm = math.exp(beta/T)
    expterm0 = math.exp(beta/T0)
    expm1 = expterm - 1.0
    expm10 = expterm0 - 1.0

    #//
    #// zero out integrals

    cv = 0
    cvint = 0
    cvtint = 0

    cv = G1 + G2*T13 + G3*T23 + G4*T

    cv = cv + G5*expterm*beta*beta/(T*T*expm1*expm1)

    #//
    #// calcualte INT[T0,T] {Cv dT}

    cvint = G1*(T-T0) + (3*G2*(T43 - T043)/4) + (3*G3*(T53 - T053)/5) \
                            + 0.5*G4*(T2 - T02) + G5*((beta/expm1) - (beta/expm10))

    #//
    #// calculate INT[T0,T] {Cv/T dT}

    cvtint = G1*(math.log(T/T0)) + (3*G2*(T13 - T013)) + (3*G3*(T23 - T023)/2) + G4*(T - T0)

    term = math.log(expm1)*expterm - math.log(expm1) - math.log(expterm)*expterm
    term = term/expm1

    cvtint = cvtint - G5*term

    term = math.log(expm10)*expterm0 - math.log(expm10) - math.log(expterm0)*expterm0
    term = term/expm10
    cvtint = cvtint + G5*term
       
    return cv,cvint, cvtint;

#
# C-4

def cv4(thermconst):
    """
//
// STANFORD SPECIFIC HEAT EQUATION (C-4)
//
//      Cv = Sum(i=1,7) {G_i T^(i-4)} + G_8 exp(beta/T)  [(beta/T)/(exp(beta/T) - 1)]^2
//

    """

    T = thermconst['Temperature']

    #
    # create the G array

    G = []
    #for i in range(1, N+1):
    #    key = 'StanfordG' + str(i)
    #    G.append(thermconst[key])

    loadArray(thermconst,G,'StanfordG',1,9)
    beta = thermconst['StanfordGbeta']
    T0 = thermconst['StanfordT0']

    G1 = thermconst['StanfordG1']
    G2 = thermconst['StanfordG2']
    G3 = thermconst['StanfordG3']
    G4 = thermconst['StanfordG4']
    G5 = thermconst['StanfordG5']
    G6 = thermconst['StanfordG6']
    G7 = thermconst['StanfordG7']
    G8 = thermconst['StanfordG8']


    #//
    #// powers of T

    T2 = T*T
    T3 = T2*T
    T4 = T3*T
    T02 = T0*T0
    T03 = T02*T0
    T04 = T03*T0

    expterm = math.exp(beta/T);
    expterm0 = math.exp(beta/T0);
    expm1 = expterm - 1.0;
    expm10 = expterm0 - 1.0;

    #//
    #// zero out integrals

    cv = 0
    cvint = 0
    cvtint = 0

    cv = (G1/T3) + (G2/T2) + (G3/T) + G4 + G5*T + G6*T2 + G7*T3

    cv = cv + G8*expterm*beta*beta/(T*T*expm1*expm1)

    #//
    #// calcualte INT[T0,T] {Cv dT}

    cvint = ((-G1/2)*((1/T2) - (1/T02))) - (G2*((1/T)-(1/T0)))       \
                            + G3*(math.log(T) - math.log(T0)) + G4*(T- T0) +   \
                            + (G5*(T2 - T02)/2) + (G6*(T3 - T03)/3)  \
                            + (G7*(T4 - T04)/4) + (G8*beta/(expm1)   \
                            - (G8*beta/expm10))

    #//
    #// calculate INT[T0,T] {Cv/T dT}

    cvtint = (-G1*((1/T3) - (1/T03))/3) - (G2*((1/T2) - (1/T02))/2) - (G3*((1/T) - (1/T0)))             \
                            + G4*(math.log(T) - math.log(T0)) + G5*(T - T0) + (G6*(T2-T02)/2) + (G7*(T3-T03)/3)  

    term = (1/(beta*T)) + (1/(beta*T*(expterm-1))) - (math.log(expterm-1)/(beta*beta))   \
                    - (2*math.log(T)/(beta*beta)) + (math.log(T2)/(beta*beta))

    cvtint = cvtint + G8*beta*beta*term

    term = (1/(beta*T0)) + (1/(beta*T0*(expterm0-1))) - (math.log(expterm0-1)/(beta*beta))  \
                    - (2*math.log(T0)/(beta*beta)) + (math.log(T02)/(beta*beta))

    cvtint = cvtint - G8*beta*beta*term
        
    return cv, cvint, cvtint

#
# C-5

def cv5(thermconst):


    """
//
// STANFORD SPECIFIC HEAT EQUATION (C-5)
//
//      Cv = G1         T <= T1
//
//
//      Cv = Sum(i=1,12) {G_i [ln(T/T_1)]^(i-1) T1 <= T <= T2
//      Cv = Sum(i=13,17) {G_i [ln(T/T_1)]^(i-1)        T > T2
//


/////////////////////////////////////////////////
// Function StanfordC5 will calculate the 
// specific heat, and two integrals for the 
// specific heat
//

    """
    T = thermconst['Temperature']
    T0 = thermconst['StanfordT0']
    T1 = thermconst['StanfordT1']
    T2 = thermconst['StanfordT2']
    #
    # create the G array

    G = []
    #for i in range(1, N+1):
    #    key = 'StanfordG' + str(i)
    #    G.append(thermconst[key])

    loadArray(thermconst,G,'StanfordG',1,18)
    #print (G)

    #//
    #// zero out integrals

    cv = 0
    cvint = 0
    cvtint = 0

    #
    # Now switch based on the value for the temperature


    #//
    #// Region I 

    #//////////////////////////////
    #// BEGIN if (T <= T1)
    #//////////////////////////////

    if (T <= T1):
        cv = G[0]
        cvint = cv*(T - T0)
        cvtint = cv*math.log(T/T0)


    #//////////////////////////////
    #// END if (T <= T1)
    #//////////////////////////////

    #//////////////////////////////////
    #// BEGIN if ((T>T1) && (T<=T2)) 
    #//////////////////////////////////

    if (T>=T1):
        if T > T2:
            T = T2
        cv = G[0]
        cvint = cv*(T1 - T0)
        cvtint = cv*math.log(T1/T0)
        TT1 = T/T1
        cvint = cvint + cv*(T-T1)
        Iprev = (T-T1)

        cvtint = cvtint + cv*(math.log(T/T1))

        #//////////////////////////////
        #// BEGIN for (i=1;i<12;i++)
        #//////////////////////////////

        for i in range(1,12):
            cv = cv + G[i]*(math.log(TT1)**i)
            term = T*(math.log(TT1)**i)
            term = term - i*Iprev
            cvint = cvint + (G[i])*term
            Iprev = term
            term2 = ((math.log(TT1)**(i+1)))/(i+1)
            cvtint = cvtint + G[i]*term2

        #//////////////////////////////
        #// END for (i=1;i<12;i++)
        #//////////////////////////////
    
    if (thermconst['Temperature']>=T2):
        T = thermconst['Temperature']
        TT2 = T/T2
        cv = G[12]
        cvint = cvint + cv*(T-T2)
        Iprev = (T-T2)
        cvtint = cvtint + cv*math.log(TT2)
                
        #//////////////////////////////
        #// BEGIN for (i=13;i<17;i++)
        #//////////////////////////////

        for i in range(13,17):
            cv = cv + (G[i])*(math.log(TT2)**(i-12))
            term = T*(math.log(TT2)**(i-12))
            term = term - (i-12)*Iprev
            cvint = cvint + (G[i])*term
            Iprev = term
            term2 = (math.log(TT2)**(i-11))/(i-11)
            cvtint = cvtint + (G[i])*term2

        #//////////////////////////////
        #// END for (i=13;i<17;i++)
        #//////////////////////////////

    return cv, cvint, cvtint


#
# C-6
#

def cv6(thermconst):
    """
 STANFORD SPECIFIC HEAT EQUATION (C-6)

      Cv = Sum(i=1,6) {G_i T^(i-2)}

 Function cv6 will calculate the 
 specific heat, and two integrals for the 
 specific heat

    """
    T = thermconst['Temperature']
    T0 = thermconst['StanfordT0']

    #
    # create the G array

    G = []
    
    #for i in range(1, 7):
    #    key = 'StanfordG' + str(i)
    #    G.append(thermconst[key])
    
    loadArray(thermconst, G, 'StanfordG', 1, 7)
        
    # print (G)

    #
    # calculate cv

    cv = cvint = cvtint = 0

    for i in range(0,6):
        cv = cv + G[i]  * (T**i-1)

    #
    # calculate INT[T0,T] {Cv dT}

    cvint = G[0] * math.log(T/T0)

    for i in range(1,6):
        cvint = cvint + (G[i]/(i)) * ((T**i) - (T0**i))

    #
    # calculate INT[T0,T] {Cv/T dT}

    cvtint = -G[0] * ((1/T) - (1/T0))
    cvtint = cvtint + G[1] * math.log(T/T0)

    for i in range(2,6):
        cvtint = cvtint + (G[i]/(i-1)) * (T**(i-1) - T0**(i-1))

    return cv, cvint, cvtint
 
#
# C-7
#

def cv7(thermconst):
    """
    //
    //      STANFORD SPECIFIC HEAT EQUATION (C-7)
    //
    //      Cv = G1 + (G2 exp(-beta/T) / T^2)
    //
    """

    G1 = thermconst['StanfordG1']
    G2 = thermconst['StanfordG2']
    beta = thermconst['StanfordGbeta']

    #//
    #// Set T,T0
    T = thermconst['Temperature']
    T0 = thermconst['StanfordT0']

    #//
    #// zero out integrals

    cv = 0
    cvint = 0
    cvtint = 0


    expterm = math.exp(-beta/T)
    expterm0 = math.exp(-beta/T0)
        
    cv = G1 + (G2*expterm/(T*T))

    cvint = G1*(T- T0) + (G2*(expterm - expterm0)/beta)

    cvtint = G1*math.log(T/T0) + G2*((1/(beta*beta)) + (1/(beta*T)))*(expterm - expterm0)

    return cv, cvint, cvtint


#
#  RHOSAT EQUATIONS
# 

#
# D-1

def rhosat1(thermconst):
    """
//
// STANFORD SATURATED LIQUID DENSITY EQUATION (D-2)
//
// rho_f = SUM(i=1,6) {D_i X^((i-1)/3)
//
//      X = 1 - T/T_c
//
//
    """
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    #
    # create the D array

    D = []
    #for i in range(1, 7):
    #    key = 'StanfordD' + str(i)
    #    D.append(thermconst[key])

    loadArray(thermconst,D,'StanfordD',1,8)

    
    # print (D)
    
    X = 1 - (T/Tc)
    rhof = 0
    for i in range(0,5):
        rhof = rhof + D[i] * X**(i/3)

    rhof = rhof + D[5] * X**(1/2) + D[6] * X*X
    return rhof



#
# D-2
#

def rhosat2(thermconst):
    """
    
 STANFORD SATURATED LIQUID DENSITY EQUATION (D-2)

 rho_f = SUM(i=1,6) {D_i X^((i-1)/3)

      X = 1 - T/T_c


    """
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    #
    # create the D array

    D = []
    #for i in range(1, 7):
    #    key = 'StanfordD' + str(i)
    #    D.append(thermconst[key])

    loadArray(thermconst,D,'StanfordD',1,7)
    
    # print (D)
    
    X = 1 - (T/Tc)
    rhof = 0
    for i in range(0,6):
        rhof = rhof + D[i] * X**(i/3)

    return rhof

#
# D-3

def rhosat3(thermconst):
    """

    """

    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    rhoc = 1/thermconst['vCrit']
    D1 = thermconst['StanfordD1']
    D2 = thermconst['StanfordD2']
    D3 = thermconst['StanfordD3']
    alpha = thermconst['StanfordDalpha']
    rhot = thermconst['Stanfordrhot']
    Tt = thermconst['StanfordTt']
    #
    # now set the parameter w

    w = (Tc - T)/(Tc- Tt)

    rhof = D1*(1.0-(w**(2.0/3.0))) + D2*(1.0-(w**(4.0/3.0))) + D3*(1.0 - w*w)

    rhof = math.exp(rhof)

    rhof = rhof*(w**alpha)
    rhof = rhoc + rhof*(rhot - rhoc)
    return rhof

#
# D-4
def rhosat4(thermconst):
    """
//
// STANFORD SATURATED LIQUID DENSITY EQUATION (D-4)
//
// rho_f = rho_c + D_1 X^(alpha) + SUM(i=1,6) {D_(i+1) X^(1+((i-1)/3))
//
//      X = 1 - T/T_c
//
//


/////////////////////////////////////////////////
// Function StanfordD4 will calculate the 
// saturated liquid density given the temperature
// as based on equation (D-4) in the Stanford
// properties book
//
//      INPUT:
//      
//              double          *consts         pointer to constants 
//                                                              array
//              double          T                       temperature                             [K]
//
//      OUTPUT:
//
//              NULL
//
//      RETURN:
//
//              double          rhof            saturated liquid                [kg/m^3]
//
//
    """

    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    alpha = thermconst['StanfordDalpha']
    rhoc = 1/thermconst['vCrit']
    #
    # create the D array

    D = []
    #for i in range(1, 7):
    #    key = 'StanfordD' + str(i)
    #    D.append(thermconst[key])

    loadArray(thermconst,D,'StanfordD',1,7)

    
    # print (D)
    
    X = 1 - (T/Tc)
    rhof = 0

    for i in range(0,6):
        rhof = rhof + D[i]*(X**(1+(i/3.0)))

    rhof = rhof + (D[0]*(X**alpha)) + rhoc

    return rhof

def rhosat5(thermconst):
    """
// rho_f = rho_c [ 1 + SUM(i=1,N) {D_i X^(i/3)}]
//
//      X = 1 - T/T_c
    """
    T = thermconst['Temperature']
    Tc = thermconst['TCrit']
    rhoc = 1/thermconst['vCrit']
    #
    # create the D array

    D = []
    #for i in range(1, 7):
    #    key = 'StanfordD' + str(i)
    #    D.append(thermconst[key])

    loadArray(thermconst,D,'StanfordD',1,9)

    
    # print (D)
    
    X = 1 - (T/Tc)
    rhof = 0
    for i in range (0,8):
        rhof = rhof + D[i]*X**((i+1)/3)

    rhof = rhoc*(1 + rhof)
    return rhof

def rhosat6(thermconst):
    """
// rho_f = SUM(i=1,4) {D_i (T - T_p)^(i-1)]
    """
    T = thermconst['Temperature']
    Tp = thermconst['StanfordTp']
    rhoc = 1/thermconst['vCrit']
    #
    # create the D array

    D = []
    #for i in range(1, 7):
    #    key = 'StanfordD' + str(i)
    #    D.append(thermconst[key])

    loadArray(thermconst,D,'StanfordD',1,5)

    X = T - Tp;
    rhof = 0;

    for i in range(0,4):
        rhof = rhof + D[i]*X**i

    return rhof

def rhosat7(thermconst):
    """
    rho_f = D1 + D2 X + D3 X^alpha + D4 X^2
    """
    T = thermconst['Temperature']
    alpha = thermconst['StanfordDalpha']
    Tc = thermconst['TCrit']
    rhoc = 1/thermconst['vCrit']
    
    D1 = thermconst['StanfordD1']
    D2 = thermconst['StanfordD2']
    D3 = thermconst['StanfordD3']
    D4 = thermconst['StanfordD4']
    alpha = thermconst['StanfordDalpha']

    X = 1.0 - (T/Tc)
    rhof = D1 + D2*X + D3*X**alpha + D4*X*X

    return rhof

def rhosat8(thermconst):
    """
// rho_f = D1 + D2 + D3 X^2
//
//      X = T - Tp
    """
    T = thermconst['Temperature']
    Tp = thermconst['StanfordTp']

    D1 = thermconst['StanfordD1']
    D2 = thermconst['StanfordD2']
    D3 = thermconst['StanfordD3']

    X = T - Tp
    rhof = D1 + D2*X + D3*X*X

    return rhof



#
# UTILITY FUNCTIONS
# 

def loadArray(dict, array, arrayString, startVal, endVal):
    """
    Use this function to load values into 
    an array

    Note endval is not inclusive in the range for the indicies
    """

    for i in range(startVal, endVal):
        try:
            key = arrayString + str(i)
            array.append(dict[key])
        except:
            array.append(0.0)

    return

def load2DArray(dict, array, arrayString, startVal1, endVal1, startVal2, endVal2):
    """
    Use this function to load values into 
    a 2-D array

    Note the endvals are not inclusive in the range for the indicies
    """

    for i in range(startVal1, endVal1):
        array.append([])
        for j in range(startVal2, endVal2):
            try:
                key = arrayString + str(i) + "," + str(j)
                array[i-1].append(dict[key])
            except:
                array[i-1].append(0.0)

    return

#
# ADDITIONAL EQUATIONS
#

def vtp(thermo):
    thermo['pActual'] = thermo['Pressure']
    v0 = thermo['SpVol']
    thermtuple = 'Stanford vtp',thermo
    #print (thermtuple, type(thermtuple)
    v = scipy.optimize.newton(vftp, v0, None, thermtuple)
    return v

#
# def vftp: used to calculate the pressure
# from guessed value of v in secant/newton's method

def vftp(v,str, thermo):
    vsave = thermo['SpVol']
    thermo['SpVol'] = v
    #print ('vftp, v = ', v)
    pcalc = thermo['pvt'](thermo)
    thermo['SpVol'] = vsave
    return pcalc - thermo['pActual']

def tpv(thermo):
    thermo['pActual'] = thermo['Pressure']
    T0 = thermo['Temperature']
    thermtuple = 'Stanford tpv',thermo
    #print (thermtuple, type(thermtuple))
    T = scipy.optimize.newton(tfpv, T0, None, thermtuple)
    return T


def tfpv(T,str, thermo):
    Tsave = thermo['Temperature']
    thermo['Temperature'] = T
    #print ('Stanford in tfpv', T)
    pcalc = thermo['pvt'](thermo)
    thermo['Temperature'] = Tsave
    return pcalc - thermo['pActual']

def tsat(thermo):
    thermo['pActual'] = thermo['Pressure']
    Tsat = thermo['Temperature']
    thermtuple = 'Stanford tsat',thermo
    #print (thermtuple, type(thermtuple))
    T = scipy.optimize.newton(tsatfp, Tsat, None, thermtuple)
    return T

def tsatfp(T,str, thermo):
    Tsave = thermo['Temperature']
    thermo['Temperature'] = T
    #print ('Stanford in tsatfp', T)
    pcalc = thermo['psat'](thermo)
    thermo['Temperature'] = Tsave
    return pcalc - thermo['pActual']

def vts(thermconst):
    thermconst['sActual'] = thermconst['SpEntropy']
    v0 = thermconst['SpVol']
    thermtuple = 'Stanford vts',thermconst
    #print (thermtuple, type(thermtuple)
    v = scipy.optimize.newton(vfts, v0, None, thermtuple)
    return v

def vfts(v,str, thermo):
    vsave = thermo['SpVol']
    thermo['SpVol'] = v
    P = thermo['pvt'](thermo)
    thermo['Pressure'] = P
    (cx, cxint, cxtint) = thermo['cvfunc'](thermo)
    (udvint, sdvint) = thermo['dvint'](thermo)
    scalc = cxtint + sdvint + thermo['sRef']
    #print ('vftp, v = ', v)
    thermo['SpVol'] = vsave
    return scalc - thermo['sActual']

def vth(thermconst):
    thermconst['hActual'] = thermconst['SpEnthalpy']
    v0 = thermconst['SpVol']
    thermtuple = 'Stanford vth',thermconst
    #print (thermtuple, type(thermtuple)
    v = scipy.optimize.newton(vfth, v0, None, thermtuple)
    return v

def vfth(v,str, thermo):
    vsave = thermo['SpVol']
    thermo['SpVol'] = v
    P = thermo['pvt'](thermo)
    thermo['Pressure'] = P
    (cx, cxint, cxtint) = thermo['cvfunc'](thermo)
    (udvint, sdvint) = thermo['dvint'](thermo)
    u = cxint + udvint + thermo['xRef']
    hcalc = u + thermo['Pressure']*thermo['SpVol']
    #print ('vftp, v = ', v)
    thermo['SpVol'] = vsave
    return hcalc - thermo['hActual']

def tph(thermconst):
    thermconst['hActual'] = thermconst['SpEnthalpy']
    T0 = thermconst['Temperature']
    thermtuple = 'Stanford tph',thermconst
    #print (thermtuple, type(thermtuple)
    T = scipy.optimize.newton(tfph, T0, None, thermtuple)
    return T

def tfph(T,str, thermo):
    Tsave = thermo['Temperature']
    vsave = thermo['SpVol']
    thermo['Temperature'] = T
    thermo['SpVol'] = ideal.vtp(thermo)
    v  = vtp(thermo)
    thermo['SpVol'] = v
    (cx, cxint, cxtint) = thermo['cvfunc'](thermo)
    (udvint, sdvint) = thermo['dvint'](thermo)
    u = cxint + udvint + thermo['xRef']
    h = u + thermo['Pressure']*thermo['SpVol']
    thermo['Temperature'] = Tsave
    thermo['SpVol'] = vsave
    return h - thermo['hActual']

def tps(thermconst):
    thermconst['sActual'] = thermconst['SpEntropy']
    T0 = thermconst['Temperature']
    thermtuple = 'Stanford tps',thermconst
    #print (thermtuple, type(thermtuple)
    T = scipy.optimize.newton(tfps, T0, None, thermtuple)
    return T

def tfps(T,str, thermo):
    Tsave = thermo['Temperature']
    vsave = thermo['SpVol']
    thermo['Temperature'] = T
    thermo['SpVol'] = ideal.vtp(thermo)
    v  = vtp(thermo)
    thermo['SpVol'] = v
    (cx, cxint, cxtint) = thermo['cvfunc'](thermo)
    (udvint, sdvint) = thermo['dvint'](thermo)
    scalc = cxtint + sdvint + thermo['sRef']
    thermo['Temperature'] = Tsave
    thermo['SpVol'] = vsave
    return scalc - thermo['sActual']

def setStanfordEqns(thermconst):
    """
    This function will set the appropriate Stanford equations in the 
    the dictionary
    """

    #
    #
    # Set the equations

    #
    # PVT
    #

    thermconst['pvt'] = getattr(sys.modules['StanfordEqns'],"pvt" + str(thermconst['Stanford_pvt']))
    thermconst['dvint'] = getattr(sys.modules['StanfordEqns'],"dvint" + str(thermconst['Stanford_pvt']))
    thermconst['vtp'] = vtp
    thermconst['tpv'] = tpv
    thermconst['tsat'] = tsat

    #
    # Psat
    #
    thermconst['psat'] = getattr(sys.modules['StanfordEqns'],"psat" + str(thermconst['Stanford_psat']))
    thermconst['dpsatdt'] = getattr(sys.modules['StanfordEqns'],"dpsatdt" + str(thermconst['Stanford_psat']))

    #
    # rhosat
    #
    thermconst['rhosat'] = getattr(sys.modules['StanfordEqns'],"rhosat" + str(thermconst['Stanford_rhosat']))
    #
    # cv
    #
    thermconst['cvfunc'] = getattr(sys.modules['StanfordEqns'],"cv" + str(thermconst['Stanford_cv']))

    return


#
# calc hfg, sfg

def calchfgsfg(thermconst):
    T = thermconst['Temperature']
    dpsatdt = thermconst['dpsatdt'](thermconst)
    vg = vtp(thermconst)
    rhosat = thermconst['rhosat'](thermconst)
    vf = 1/rhosat
    vfg = vg - vf
    hfg = T * vfg * dpsatdt
    sfg = hfg / T
    return hfg, sfg

#
# Calculate u,h,s,x

def calcuhsx(thermconst):
    """
    Function calcuhsx will calculate the 
    thermodynamic properties u,h,s,x given 
    T, P, v and the region.  If the region is the 
    2-phase region then the user may supply
    one other property (h, s, or x) based on 
    user key '2phaseprop'
    """

    region = thermconst['Region']
    (cx, cxint, cxtint) = thermconst['cvfunc'](thermconst)
    (udvint, sdvint) = thermconst['dvint'](thermconst)
    u = cxint + udvint + thermconst['xRef']
    h = u + thermconst['Pressure'] * thermconst['SpVol']
    s = cxtint + sdvint + thermconst['sRef']

    if region in ('SHV','SC'):
        thermconst['SpEnergy'] = u
        thermconst['SpEnthalpy'] = h
        thermconst['SpEntropy'] = s
        thermconst['Quality'] = -1
    elif region == '2phase':

        psat = thermconst['psat'](thermconst)
        vsave = thermconst['SpVol']
        #
        # calculate vg
        thermconst['Pressure'] = psat
        thermconst['SpVol'] = ideal.vtp(thermconst)
        vg = vtp(thermconst)
        rhosat = thermconst['rhosat'](thermconst)
        vf = 1/rhosat
        
        thermconst['SpVol'] = vg
        (cx, cxint, cxtint) = thermconst['cvfunc'](thermconst)
        (udvint, sdvint) = thermconst['dvint'](thermconst)
        u = cxint + udvint + thermconst['xRef']
        h = u + thermconst['Pressure'] * thermconst['SpVol']
        s = cxtint + sdvint + thermconst['sRef']
        (hfg, sfg) = calchfgsfg(thermconst)
        thermconst['SpVol'] = vsave

        if thermconst['2phaseprop'] == 'h':
            hactual = thermconst['SpEnthalpy']
            #print (hactual, hfg)
            x = 1 + ((hactual - h)/hfg)
            v = vf + x * (vg - vf)
            thermconst['SpVol'] = v
            thermconst['SpEnergy'] = hactual - thermconst['Pressure'] * v
            thermconst['SpEntropy'] = s + (x - 1) * sfg
            thermconst['Quality'] = x
        elif thermconst['2phaseprop'] == 's':
            sactual = thermconst['SpEntropy']
            x = 1 + ((sactual - s)/sfg)
            v = vf + x * (vg - vf)
            thermconst['SpEnthalpy'] = h  + (x-1)*hfg
            thermconst['SpEnergy'] = thermconst['SpEnthalpy'] - thermconst['Pressure'] * v
            thermconst['Quality'] = x
        elif thermconst['2phaseprop'] == 'v':
            vactual = thermconst['SpVol']
            x = 1 + ((vactual - vg)/(vg - vf))
            if x > 1.0:
                pass
            thermconst['SpEnthalpy'] = h  + (x-1)*hfg
            thermconst['SpEnergy'] = thermconst['SpEnthalpy'] - thermconst['Pressure'] * vactual
            thermconst['SpEntropy'] = s + (x - 1) * sfg
            thermconst['Quality'] = x
        elif thermconst['2phaseprop'] == 'x':
            x = thermconst['Quality']
            #print (vf, vg)
            v = vf + x * (vg - vf)
            thermconst['SpVol'] = v
            thermconst['SpEnthalpy'] = h  + (x-1)*hfg
            thermconst['SpEnergy'] = thermconst['SpEnthalpy'] - thermconst['Pressure'] * v
            thermconst['SpEntropy'] = s + (x - 1) * sfg
        elif region == 'CL':
            print ('need CL search')
    return

