"""
NBS steam tables EOS functions library

This file contains the NBS steam tables

03-05-2013

sjpfu

"""

import ThermoConst as thermc
import math
import scipy.optimize
#import IdealGas as ideal


def NBSprops(thermconst):
    units(thermconst)
    job = thermconst['job']
    thermconst['NBSrho'] = 0
    thermconst['SpHeatp'] = 0
    thermconst['dpdt'] = 0 
    thermconst['NBSdpdrho'] = 0
    thermconst['NBSConvertFlag'] = "toInternal"
    convertUnits(thermconst)
    setData(thermconst)
    gascon = thermconst['NBSgascon']

    if job == 1:
    #
    # case 1, Inputs are P, T
        p = thermconst['Pressure']
        t = thermconst['Temperature']
#        p = p/thermconst['NBSfp']
#        convertTemperature(thermconst)
        dgss = p/(t*0.4)
        psat = 2e4
        rt = t*gascon
        thermconst['NBSdl'] = 0
        thermconst['NBSdv'] = 0
        
        bb(thermconst)
        dgss = p/(0.4*(t))
        if (t < thermconst['NBSTz']):
            thermconst['Pressure'] = psat
            psat = pcorr(thermconst)
            thermconst['Pressure'] = p
        if (p > psat):
                dgss = thermconst['NBSdl']
        
        thermconst['NBSd'] = dgss
        #bb(thermconst)
        rho = dfind(thermconst)
        thermconst['NBSd'] = rho
        therm(thermconst)
        

        if (p < psat):
            thermconst['Region'] = "SHV"
        elif (p > psat):
            thermconst['Region'] = "CL"
        else:
            thermconst['Region'] = "2phase"

        calcValues(thermconst)

    #
    # case 2, Inputs are T, v

    elif (job == 2):
        v = thermconst['SpVol']
        rho = 1/v
#        rho = rho*thermconst['NBSfd']
        print ('v,rho in = ',1/rho,rho)
        t = thermconst['Temperature']
#        thermconst['NBSConvertFlag'] = "toInternal"
#        convertTemperature(thermconst)
        #rt = t*gascon
        vf = 0 
        vg = 0

        if (t < thermconst['NBSTz']):
            thermconst['NBSdl'] = 0
            thermconst['NBSdv'] = 0
            psat = 2e4
            thermconst['Pressure'] = psat
            psat = pcorr(thermconst)
            vf = 1/thermconst['NBSdl']
            vg = 1/thermconst['NBSdv']
            #print ('vf,vg = ',1/thermconst['NBSdl'],1/thermconst['NBSdv'])
            #print ('rhof,rhog = ',thermconst['NBSdl'],thermconst['NBSdv'])
        
        if (vf < v) and (v < vg):
            thermconst['Region'] = "2phase"
            thermconst['2phaseprop'] = 'v'
            thermconst['vActual'] = v
            calc2Phase(thermconst)
        elif (v > vg):
            thermconst['Region'] = "SHV"
            thermconst['NBSd'] = rho
            calcValues(thermconst)
        else:
            thermconst['Region'] = "CL"
            thermconst['NBSd'] = rho
            calcValues(thermconst)

    elif (job == 3):
        #
        # case 3, Inputs are T,s
        t = thermconst['Temperature']
        s = thermconst['SpEntropy']
        thermconst['sActual'] = s
        vf = 0 
        vg = 0

        if (t < thermconst['NBSTz']):
            thermconst['NBSdl'] = 0
            thermconst['NBSdv'] = 0
            psat = 2e4
            thermconst['Pressure'] = psat
            psat = pcorr(thermconst)
            vf = 1/thermconst['NBSdl']
            vg = 1/thermconst['NBSdv']

        thermconst['Region'] = "2phase"
        thermconst['2phaseprop'] = 'x'
        thermconst['Quality'] = 1
        calc2Phase(thermconst)
        sg = thermconst['SpEntropy']
        thermconst['Quality'] = 0
        calc2Phase(thermconst)
        sf = thermconst['SpEntropy']
        thermconst['SpEntropy'] = s
        if (sf <= s) and (s <= sg):
            #
            # 2-phase region
            #print ('sf,s,sg = ',sf,s,sg)
            thermconst['2phaseprop'] = 's'
            thermconst['SpEntropy'] = s
            calc2Phase(thermconst)
            thermconst['SpEntropy'] = s
            rho = thermconst['NBSrho'] 
        elif (s > sg):
            #
            # SHV, call dts using dv as guess
            thermconst['NBSrho'] = thermconst['NBSdv']
            rho = dts(thermconst)
            thermconst['NBSrho'] = rho
            thermconst['Region'] = "SHV"
            calcValues(thermconst)
        else:
            #
            # CL, call dts using dl as guess
            thermconst['NBSrho'] = thermconst['NBSdl']
            rho = dts(thermconst)
            thermconst['NBSrho'] = rho
            thermconst['Region'] = "CL"
            calcValues(thermconst)


    elif (job == 4):
        #
        # case 4, Inputs are T,x 
        t = thermconst['Temperature']
        x = thermconst['Quality']
        psat = 2e4
        rt = t*gascon
        thermconst['NBSdl'] = 0
        thermconst['NBSdv'] = 0
        
        bb(thermconst)
        dgss = psat/(0.4*(t))
        if (t < thermconst['NBSTz']):
            thermconst['Pressure'] = psat
            psat = pcorr(thermconst)
        else:
            print ('T > Tcrit')

        thermconst['Region'] = "2phase"
        thermconst['2phaseprop'] = 'x'
        thermconst['xActual'] = x
        calc2Phase(thermconst)
        rho = thermconst['NBSrho'] 

    elif (job == 5):
        #
        # case 5, inputs are P, v
        p = thermconst['Pressure']
        v =  thermconst['SpVol']
        rho = 1/v
        thermconst['NBSdl'] = 0
        thermconst['NBSdv'] = 0
        tsat = tcorr(thermconst)
        vf = 1/thermconst['NBSdl']
        vg = 1/thermconst['NBSdv']

        if (vf <= v) and (v <= vg):
            #
            # 2 phase
            thermconst['2phaseprop'] = 'v'
            thermconst['SpVol'] = v
            thermconst['vActual'] = v
            thermconst['Temperature'] = tsat
            calc2Phase(thermconst)
            thermconst['Region'] = '2phase'
        else:
            t = tpd(thermconst)
            if (t > tsat):
                thermconst['Region'] = "SHV"
            else:
                thermconst['Region'] = "CL"
            thermconst['Temperature'] = t
            calcValues(thermconst)
        thermconst['Pressure'] = p
        thermconst['SpVol'] = v

    elif (job == 6):
        #
        # case 6, Inputs are P,h
        p = thermconst['Pressure']
        h = thermconst['SpEnthalpy']
        thermconst['hActual'] = h
        thermconst['NBSdl'] = 0
        thermconst['NBSdv'] = 0
        tsat = tcorr(thermconst)
        vf = 1/thermconst['NBSdl']
        vg = 1/thermconst['NBSdv']
        thermconst['Temperature'] = tsat
        thermconst['Region'] = "2phase"
        thermconst['2phaseprop'] = 'x'
        thermconst['Quality'] = 1
        calc2Phase(thermconst)
        hg = thermconst['SpEnthalpy']
        thermconst['Quality'] = 0
        calc2Phase(thermconst)
        hf = thermconst['SpEnthalpy']
        thermconst['SpEnthalpy'] = h
        
        if (hf <= h) and (h <= hg):
            #
            # 2-phase region
            thermconst['2phaseprop'] = 'h'
            thermconst['SpEnthalpy'] = h
            calc2Phase(thermconst)
            thermconst['SpEnthalpy'] = h
            rho = thermconst['NBSrho'] 
        elif (h > hg):
            #
            # SHV, call dts using dv as guess
            thermconst['NBSrho'] = thermconst['NBSdv']
            thermconst['Temperature'] = tsat
            thermconst['hActual'] = h
            t = tph(thermconst)
            thermconst['Temperature'] = t
            thermconst['Region'] = "SHV"
            calcValues(thermconst)
            rho = thermconst['NBSrho']
        else:
            #
            # CL, call dts using dl as guess
            thermconst['NBSrho'] = thermconst['NBSdl']
            thermconst['hActual'] = h
            thermconst['Temperature'] = tsat
            t = tph(thermconst)
            thermconst['Temperature'] = t
            thermconst['Region'] = "CL"
            calcValues(thermconst)
            rho = thermconst['NBSrho']

        thermconst['Pressure'] = p
        thermconst['SpEnthalpy'] = h

    elif (job == 7):
        #
        # case 7, Inputs are P,s
        p = thermconst['Pressure']
        s = thermconst['SpEntropy']
        thermconst['sActual'] = s
        thermconst['NBSdl'] = 0
        thermconst['NBSdv'] = 0
        tsat = tcorr(thermconst)
        vf = 1/thermconst['NBSdl']
        vg = 1/thermconst['NBSdv']
        thermconst['Temperature'] = tsat
        thermconst['Region'] = "2phase"
        thermconst['2phaseprop'] = 'x'
        thermconst['Quality'] = 1
        calc2Phase(thermconst)
        sg = thermconst['SpEntropy']
        thermconst['Quality'] = 0
        calc2Phase(thermconst)
        sf = thermconst['SpEntropy']
        thermconst['SpEntropy'] = s
        print (s,sf,sg)
        if (sf <= s) and (s <= sg):
            #
            # 2-phase region
            #print ('sf,s,sg = ',sf,s,sg)
            thermconst['2phaseprop'] = 's'
            thermconst['SpEntropy'] = s
            calc2Phase(thermconst)
            thermconst['SpEntropy'] = s
            rho = thermconst['NBSrho'] 
        elif (s > sg):
            #
            # SHV, call dts using dv as guess
            thermconst['NBSrho'] = thermconst['NBSdv']
            thermconst['Temperature'] = tsat
            thermconst['sActual'] = s
            t = tps(thermconst)
            thermconst['Temperature'] = t
            #thermconst['NBSrho'] = rho
            thermconst['Region'] = "SHV"
            calcValues(thermconst)
            rho = thermconst['NBSrho']
        else:
            #
            # CL, call dts using dl as guess
            thermconst['NBSrho'] = thermconst['NBSdl']
            thermconst['sActual'] = s
            thermconst['Temperature'] = tsat
            t = tps(thermconst)
            print ('CL: ',t)
            thermconst['Temperature'] = t
            thermconst['Region'] = "CL"
            calcValues(thermconst)
            rho = thermconst['NBSrho']

        thermconst['Pressure'] = p
        thermconst['SpEntropy'] = s


    elif (job == 8):
        #
        # case 8, Inputs are P, x
        p = thermconst['Pressure']
        x = thermconst['Quality']

        thermconst['NBSdl'] = 0
        thermconst['NBSdv'] = 0
        tsat = tcorr(thermconst)
        thermconst['Temperature'] = tsat
        bb(thermconst)

        thermconst['Region'] = "2phase"
        thermconst['2phaseprop'] = 'x'
        thermconst['xactual'] = x
        calc2Phase(thermconst)
        thermconst['Pressure'] = p
        rho = thermconst['NBSrho'] 



        #bb(thermconst)
        #qq(thermconst)
        #zdum = base(thermconst)
        #p = rt*rho*thermconst['NBSz'] + thermconst['NBSq']
        #dq = rt*(thermconst['NBSz']+thermconst['NBSy']*thermconst['NBSdz'])+thermconst['NBSq5']
        #therm(thermconst)        


    thermconst['NBSrho'] = rho
    thermconst['NBSConvertFlag'] = "toExternal"
    convertUnits(thermconst)
    #convertTemperature(thermconst)
    #thermconst['Pressure'] = p*thermconst['NBSfp']
    #thermconst['SpVol'] = 1/rho
    #thermconst['SpEnthalpy'] = thermconst['NBShd']*rt*thermconst['NBSfh']
    #thermconst['SpEnergy'] = thermconst['NBSud']*rt*thermconst['NBSfh']
    #thermconst['SpEntropy'] = gascon*thermconst['NBSsd']*thermconst['NBSfh']*thermconst['NBSft']
    #thermconst['SpHeatp'] = thermconst['NBScpd']*gascon*thermconst['NBSfh']*thermconst['NBSft']
    #thermconst['dpdt'] = thermconst['NBSdpdt']*thermconst['NBSfp']*thermconst['NBSft']
    #thermconst['NBSdpdrho'] = thermconst['NBSdpdrho']*thermconst['NBSfd']*thermconst['NBSfp']


    return


def  bb(thermconst):
    """
    Function bb will calcualte the the B's of eqns 3,4
    using the constant coefficients and the T
    b's are calculated in cm3/g
    """
    v = []
    v.append(1)
    #
    # temperature parameter
    t = thermconst['Temperature']

    #
    # aconst block
    #
    # function requires tz

    tz = thermconst['NBSTz']

    #
    # bconst block
    #
    # reqruires bp and bq 

    bp = thermconst['NBSbp']
    bq = thermconst['NBSbq']

    for i in range (1,10):
        v.append(tz*v[i-1]/t)

    b1 = bp[0] + bp[1]*(math.log(1.0/v[1]))
    b2 = bq[0]
    b1t = bp[1]*v[1]/tz
    b2t = 0.0
    b1tt = 0.0
    b2tt = 0.0

    for i in range(2,10):
        b1 = b1 + bp[i]*v[i-1]
        b2 = b2 + bq[i]*v[i-1]
        b1t = b1t - (i-1)*bp[i]*v[i-1]/t
        b2t = b2t - (i-1)*bq[i]*v[i-1]/t
        b1tt = b1tt + bp[i]*(i-1)*(i-1)*v[i-1]/(t*t)
        b2tt = b2tt + bq[i]*(i-1)*(i-1)*v[i-1]/(t*t)
 
    b1tt = b1tt - b1t/t
    b2tt = b2tt - b2t/t


    #
    # set values in the ellcon block
    thermconst['NBSb1'] = b1
    thermconst['NBSb2'] = b2
    thermconst['NBSb1t'] = b1t
    thermconst['NBSb2t'] = b2t
    thermconst['NBSb1tt'] = b1tt
    thermconst['NBSb2tt'] = b2tt
    #print ('bb: ', b1,b2,b1t,b2t,b1tt,b2tt)
    return

def base(thermconst):
    """
    Function base calcualtes Z = Pbase/rho R T and also 
    Abase, Gbase, Sbase, Ubase, Hbase, CVbase, 1/(DRT) * DP/DT
    """
    #
    # Parameters from the dictionary
    t = thermconst['Temperature']
    d = thermconst['NBSd']

    #
    # aconst block
    gascon = thermconst['NBSgascon']

    #
    # ellcon block
    #
    # uses all parameters
    #
    g1 = thermconst['NBSg1']
    g2 = thermconst['NBSg2']
    gf = thermconst['NBSgf']
    b1 = thermconst['NBSb1']
    b2 = thermconst['NBSb2']
    b1t = thermconst['NBSb1t']
    b2t = thermconst['NBSb2t']
    b1tt = thermconst['NBSb1tt']
    b2tt = thermconst['NBSb2tt']


    y = 0.25*b1*d
    b2b1gf = (b2/b1) - gf
    x = 1 - y
    z0 = (1 + g1*y + g2*y*y)/(x**3)
    z = z0 + 4*y*b2b1gf
    dz0 = (g1 + 2*(g2)*(y))/(x**3)
    dz0 = dz0 + 3*(1 + g1*y + (g2*y*y)/(x**4))
    dz = dz0 + 4*b2b1gf

    ab = -math.log(x) - ((g2 - 1)/x) + (28.16666667/(x*x)) + 4*(y)*b2b1gf + 15.166666667 + math.log(d*t*(gascon)/0.101325)
    gb = ab + z
    bb2tt = t*t*b2tt
    ub = (-t*b1t*(z - 1 - d*b2)/b1) - d*t*b2t
    hb = z + ub
    tb1tb1 = t*b1t/b1
    tb1ttb1 = t*b1tt/b1
    cvb = 2*ub + (z0 - 1)*(tb1tb1*tb1tb1 - t*tb1ttb1)
    cvb = cvb - d*(bb2tt - gf*b1tt*t*t)
    cvb = cvb - tb1tb1*tb1tb1*y*dz0
    dpdtb = (z/t) + d*((dz*b1t/4) + b2t -(b2*(b1t)/(b1)))
    sb = ub - ab

    #
    # aconst block
    thermconst['NBSz'] = z
    thermconst['NBSdz'] = dz
    thermconst['NBSy'] = y

    #
    # basef block
    thermconst['NBSab'] = ab
    thermconst['NBSgb'] = gb
    thermconst['NBSsb'] = sb
    thermconst['NBSub'] = ub
    thermconst['NBShb'] = hb
    thermconst['NBScvb'] = cvb
    thermconst['NBSdpdtb'] = dpdtb

    return z

def qq(thermconst):
    """
    This function calculates residual contributions to Pressure (Q)
    Helmholtz (ar), dp/drho (q5), and also gibbs, entropy, internal energy
    enthalpy, ischoric heat capacity, and dpdt

    Inputs are T (K) and rho (g/cm3)
    """

    t = thermconst['Temperature']
    d = thermconst['NBSd']

    #
    # aconst block
    gascon = thermconst['NBSgascon']
    aa = thermconst['NBSaa']
    tz = thermconst['NBSTz']

    #
    # nconst block
    nc = thermconst['NBSnc']
    ii = thermconst['NBSii']
    g = thermconst['NBSg']
    jj = thermconst['NBSjj']

    #
    # addcon block
    adz = thermconst['NBSadz']
    aat = thermconst['NBSaat']
    atz = thermconst['NBSatz']
    aad = thermconst['NBSaad']

    

    #
    # establish req'd lists

    qr = []
    qt = []
    qzr = []
    qzt = []

    for i in range(0,12):
        qr.append(0)
        qt.append(0)
        qzr.append(0)
        qzt.append(0)


    rt = gascon*t
    qr[0] = 0
    q5  = 0
    q = 0
    ar = 0
    gr = 0
    hr = 0
    dadt = 0
    cvr = 0
    dpdtr = 0
    e = math.exp(-aa*d)
    q10 = d*d*e
    q20 = 1 - e
    qr[1] = q10
    v = tz/t
    qt[0] = t/tz
    for i in range(1,10):
        qr[i+1] = qr[i]*q20
        qt[i] = qt[i-1]*v
    
    for i in range(0,9):
        qzr[i] = qr[i+2]
        qzt[i] = qt[i+1]

    for i in range(0,nc+1):
        k = ii[i] + 1
        kz = k;
        l = jj[i]
        kmz = l
        qp = g[i]*(aa)*qr[k]*qt[l]
        q = q + qp
        q5 = q5 + aa*qp*((2/d) - aa*(1 - e*(k-1)/q20))
        ar = ar + g[i]*qr[k+1]*qt[l]/(k*q10*rt)
        dfdt = (q20**kz)*(1-l)*qzt[l]/(tz*k)
        d2f = l*dfdt
        dpt = dfdt*q10*aa*k/q20
        dadt = dadt + g[i]*dfdt
        dpdtr = dpdtr + g[i]*dpt
        cvr = cvr + (g[i]*d2f/(gascon))

    qp = 0
    q2a = 0

    for i in range(nc+1,40):
        if (g[i] != 0):
            k = ii[i]
            kz = k
            km = jj[i]
            kmz = km
            ddz = adz[i - nc - 1]
            delta = d/ddz - 1
            if (math.fabs(delta) < 1e-10):
                delta = 1e-10
            dd = delta**2
            ex1 = -delta**kz*(aad[i - nc - 1])
            dex = math.exp(ex1)*(delta**kmz)
            att = aat[i - nc - 1]
            tx = atz[i - nc - 1]
            tau = t/tx - 1
            ex2 = -att*tau*tau
            tex = math.exp(ex2)
            q10 = dex*tex
            qm = km/delta - k*(aad[i - nc - 1])*(delta**(kz-1))
            fct = qm*(d**2)*q10/ddz
            q5t = fct*(2/d + qm/ddz) - ((d/ddz)**2)*q10*(km/(delta**2) + k*(k-1)*(aad[i - nc - 1])*(delta**(kz-2)))
            q5 = q5 + q5t*(g[i])
            qp = qp + (g[i])*fct
            dadt = dadt - 2*(g[i])*att*tau*q10/tx
            dpdtr = dpdtr - 2*(g[i])*att*tau*fct/tx
            q2a = q2a + t*(g[i])*(4*att*ex2 + 2*att)*q10/(tx*tx)
            ar = ar + q10*(g[i])/rt

    sr = -dadt/(gascon)
    ur = ar + sr
    cvr = cvr + q2a/(gascon)
    q = q + qp

    #
    # set values in the resf block
    thermconst['NBSar'] = ar
    thermconst['NBSgr'] = gr
    thermconst['NBSsr'] = sr
    thermconst['NBSur'] = ur
    thermconst['NBShr'] = hr
    thermconst['NBScvr'] = cvr
    thermconst['NBSdpdtr'] = dpdtr
    #
    # set values in the qqqq block
    thermconst['NBSq'] = q
    thermconst['NBSq5'] = q5

    return

def dfind(thermconst):
    """
    Function to find the density corresponding to the 
    input pressure p (MPa) and t(K) using an initial
    guess density (g/cm3). the output density is in g/cm3
    also dp/drho is calculated
    """

    #
    # input parameters
    d = thermconst['NBSd']
    dd = d
    t = thermconst['Temperature']
    p = thermconst['Pressure']

    #
    # aconst block
    gascon = thermconst['NBSgascon']

    rt = gascon*t
    if (dd <= 0):
        dd = 1e-8
    elif (dd > 1.9):
        dd = 1.9
    l = 1

    #print ('entering dfind ',p,dd)
    while (l<=30):
        if (dd <= 0):
            dd = 1e-8
        elif (dd > 1.9):
                dd = 1.9

        thermconst['NBSd'] = dd
        qq(thermconst)
        #
        # get the new values from the qqqq block
        q = thermconst['NBSq']
        q5 = thermconst['NBSq5']

        z = base(thermconst)
        #
        #get the new values from teh aconst block
        y = thermconst['NBSy']
        dz = thermconst['NBSdz']

        #print ('q,q5,z,y,dz ',dd,q,q5,z,y,dz)

        pp = rt*dd*z + q
        dpd = rt*(z + y*dz) + q5
        skip_flag = 1
        if (dpd < 0):
            if (d >= 0.2967):
                dd = dd*1.02
            elif (d < 0.2967):
                dd = dd*0.98
            if (l <= 10):
                skip_flag = 0
        if ((dpd > 0) and (skip_flag == 1)):
            dpdx = (dpd)*1.1
            if (dpdx < 0.1):
                dpdx = 0.1
            dp = math.fabs(1 - pp/p)
            if (dp > 1e-8):
                if ((d < 0.3) or (dp > 1e-7)):
                    if ((d < 0.7) or (dp > 1e-6)):
                        x = (p - pp)/dpdx
                        if (math.fabs(x) > .1):
                            x = x*0.1/math.fabs(x)
                        dd = dd + x
                    else:
                        l = 31
                else:
                    l = 31
            else:
                l = 31
        l = l + 1

    thermconst['NBSdpdrho'] = dpd
    thermconst['NBSrho'] = dd
    #print ('exiting dfind ',dd,dpd)

    return dd

def therm(thermconst):
    """
    Function therm calculates the thermodynamic
    functions in dimensionless units (ad = a/rt,
    gd = g/rt, sd=s/r, ud=u/rt, hd=h/rt, cvd/cv/r,
    and cpd=cp/r)
    """
    #
    # set the temperature
    t = thermconst['Temperature']
    d = thermconst['NBSd']

    #
    # aconst block
    gascon = thermconst['NBSgascon']
    zb = thermconst['NBSz']
    dzb = thermconst['NBSdz']
    y = thermconst['NBSy']
    uref = thermconst['NBSuRef']
    sref = thermconst['NBSsRef']

    #
    # qqqq block
    qp = thermconst['NBSq']
    qdp = thermconst['NBSq5']
    

    #
    # call ideal function to set idf block

    ideal(thermconst)

    ai = thermconst['NBSai']
    gi = thermconst['NBSgi']
    si = thermconst['NBSsi']
    ui = thermconst['NBSui']
    hi = thermconst['NBShi']
    cvi = thermconst['NBScvi']
    cpi = thermconst['NBScpi']

    #
    # basef block
    ab = thermconst['NBSab']
    gb = thermconst['NBSgb']
    sb = thermconst['NBSsb']
    ub = thermconst['NBSub']
    hb = thermconst['NBShb']
    cvb = thermconst['NBScvb']
    dpdtb = thermconst['NBSdpdtb']

    #
    # resf block
    ar = thermconst['NBSar']
    gr = thermconst['NBSgr']
    sr = thermconst['NBSsr']
    ur = thermconst['NBSur']
    hr = thermconst['NBShr']
    cvr = thermconst['NBScvr']
    dpdtr = thermconst['NBSdpdtr']

    rt = gascon * t
    z = zb + qp/(rt*d)
    dpdd = rt*(zb + y*dzb) + qdp

    ad =  ab + ar + ai - (uref/t) + sref
    gd = ad + z
    ud = ub + ur + ui - (uref/t)
    dpdt = rt*d*dpdtb + dpdtr
    cvd = cvb + cvr + cvi
    cpd = cvd + (t*(dpdt**2)/(d*d*dpdd*gascon))
    hd = ud + z
    sd = sb + sr + si - sref
    dvdt = dpdt/(dpdd*d*d)
    cjtt = 1/(d - (t*dvdt))
    cjth = cjtt/(cpd*gascon)

    #
    # fcts block
    thermconst['NBSad'] = ad
    thermconst['NBSgd'] = gd
    thermconst['NBSsd'] = sd
    thermconst['NBSud'] = ud
    thermconst['NBShd'] = hd
    thermconst['NBScvd'] = cvd
    thermconst['NBScpd'] = cpd
    thermconst['NBSdpdt'] = dpdt
    thermconst['NBSdvdt'] = dvdt
    thermconst['NBScjtt'] = cjtt
    thermconst['NBScjth'] = cjth

    return

def ps(thermconsts):
    """
    Function ps is an approixmation to the
    vapor pressure ps as a function of temperature.
    The capor pressure agrees with the vapor pressure
    predicted by the surface to within 0.2% to within a 
    degree or so of the critical tempeature and can be
    used as an initial guess for further refinement by
    imposing the condition that gl=gv
    """

    t = thermconsts['Temperature']
    #/* set the array a */

    a = []

    a.append(-7.8889166e0)
    a.append(2.5514255e0)
    a.append(-6.716169e0)
    a.append(33.239495e0)
    a.append(-105.38479e0)
    a.append(174.35319e0)
    a.append(-148.39348e0)
    a.append(48.631602e0)

    if (t < 314.0):
        pl = 6.3573118 - 8858.843/t + 607.56335*(t**(-0.6))
        ps = 0.1*math.exp(pl)
    else:
        v = t/647.25
        w = math.fabs(1 - v)
        b = 0
        for i in range(0,8):
            z = i + 1
            b = b + a[i]*(w**((z+1)/2))
        q = b/v
        ps = 22.093*math.exp(q)

    return ps

def tsat(thermconsts):
    """
    Function tsat calcualtes teh saturation temperature
    for a given pressure by an ierative process using
    psat and tdpsdt
    """

    p = thermconsts['Pressure']

    tsat = 0

    if (p < 22.05):
        k = 0
        pl = 2.302585 + math.log(p)
        tg = 372.83 + pl*(27.7589 + pl*(2.3819 + pl*(0.24834 + pl*0.0193855)))
        while (k < 8):
                if (tg < 273.15):
                    tg = 273.15
                if (tg > 647):
                    tg = 647
                k = k + 1
                thermconsts['Temperature'] = tg
                pp = ps(thermconsts)
                dp = tdpsdt(thermconsts)
                if (math.fabs(1 - pp/p) < 1e-5):
                    k = 8
                else:
                    tg = tg*(1 + (p - pp)/dp)
    tsat = tg
    return tsat

def tdpsdt(thermconsts):
    """
    Function tdpsdt calculates 
    T * (dPs/dT)
    """

    t = thermconsts['Temperature']

    # /* set the a array */
    a = []

    a.append(-7.8889166e0)
    a.append(2.5514255e0)
    a.append(-6.716169e0)
    a.append(33.239495e0)
    a.append(-105.38479e0)
    a.append(174.35319e0)
    a.append(-148.39348e0)
    a.append(48.631602e0)
    v = t/647.25
    w = 1 - v
    b = 0
    c = 0
    for i in range (0,8):
        z = i+1
        y = a[i]*(w**((z+1)/2))
        c = c + y*(0.5 - 0.5*z - 1/v)/w
        b = b + y
    q = b/v
    tdpsdt = 22.093*math.exp(q)*c
    return tdpsdt

def ideal(thermconsts):
    """
    Function ideal calculates the thermodynamic 
    properties for water in the ideal gas state 
    from the function of H W Wooley
    """

    t = thermconsts['Temperature']
    c= []
    c.append(.19730271018e2)
    c.append(.209662681977e2)
    c.append(-.483429455355e0)
    c.append(.605743189245e1)
    c.append(22.56023885e0)
    c.append(-9.87532442e0)
    c.append( -.43135538513e1)
    c.append(.458155781e0)
    c.append(-.47754901883e-1)
    c.append(.41238460633e-2)
    c.append(-.27929052852e-3)
    c.append(.14481695261e-4)
    c.append(-.56473658748e-6)
    c.append(.16200446e-7)
    c.append(-.3303822796e-9)
    c.append(.451916067368e-11)
    c.append(-.370734122708e-13)
    c.append(.137546068238e-15)
    tt = t/1e2
    tl = math.log(tt)

    gi = -(c[0]/tt + c[1])*tl
    hi = (c[1] + c[0]*(1 - tl)/tt)
    cpi = c[1] - c[0]/tt
    for i in range(2,18):
        id = i - 5
        gi = gi - c[i]*(tt**id)
        hi = hi + c[i]*(i-5)*(tt**id)
        cpi = cpi + c[i]*(i-5)*(i-4)*(tt**id)

    ai = gi - 1
    ui = hi - 1
    cvi = cpi - 1
    si = ui - ai

    #
    # set the idf block
    thermconsts['NBSai'] = ai
    thermconsts['NBSgi'] = gi
    thermconsts['NBSsi'] = si
    thermconsts['NBSui'] = ui
    thermconsts['NBShi'] = hi
    thermconsts['NBScvi'] = cvi
    thermconsts['NBScpi'] = cpi
    return

def corr(thermconst):
    """
    Function corr will calcualte the liquid and 
    vapor densities for an input of t and p.  It also 
    calculates delg = (gl-gv)/rt for use in calculating
    the correction to the capor pressure for delg
    """

    #
    # Inputs are t and p
    t = thermconst['Temperature']
    p = thermconst['Pressure']

    #
    # aconst block
    gascon = thermconst['NBSgascon']
    zb = thermconst['NBSz']
    dzb = thermconst['NBSdz']
    y = thermconst['NBSy']
    uref = thermconst['NBSuRef']
    sref = thermconst['NBSsRef']

    #
    # dl and dv are guesses for liquid and vapor densities
    dl = thermconst['NBSdl']
    dv = thermconst['NBSdv']

    rt = gascon*t

    #delg = 0

    #print ('corr: T,P = ',t,p,dl)

    if (t <= 646.3):
        dliq = dl
        thermconst['NBSd'] = dliq
        if (dliq <= 0):
            dliq = 1.11 - 0.0004*t

        bb(thermconst)
        thermconst['NBSd'] = dliq
        dl = dfind(thermconst)
        thermconst['NBSd'] = dl                        
        therm(thermconst)
        gl = thermconst['NBSgd']
        dvap = dv
        if (dvap <= 0):
            dvap = p/rt
            thermconst['NBSdv'] = dvap
        thermconst['NBSd'] = dvap
        dv = dfind(thermconst)
        #print ('corr ',1/dl,1/dv)
        if (dv < 5e-7):
            dv = 5e-7
        gv = gl
        delg = gl - gv
    else:
        p = 0
        if (t < 647.126):
            thermconst['NBSdelg'] = 0;
            bb(thermconst)
            tau = 0.657128*((1 - t/647.126)**0.325)
            dl = 0.322 + tau
            dv = 0.322 - tau
            thermconst['NBSd'] = dv
            zdum = base(thermconst)
            qq(thermconst)
            p = rt*(dv)*(thermconst['NBSz']) + thermconst['NBSq']

    thermconst['NBSdl'] = dl
    thermconst['NBSdv'] = dv
    thermconst['NBSdelg'] = delg
    return

def pcorr(thermconst):
    """
    Function pcorr will calculate the vapor pressure, p
    and the vapor densities corresponing to the input
    temperature t
    """

    t = thermconst['Temperature']    
    #psave = thermconst['Pressure']
    gascon = thermconst['NBSgascon']
    delg = 1
    p = ps(thermconst)
    thermconst['Pressure'] = p
    #print ('pcorr: psat = ',p)
    while (delg > 1e-4):
        corr(thermconst)
        delg = thermconst['NBSdelg']
        dv = thermconst['NBSdv']
        dl = thermconst['NBSdl']
        #print ('pcorr dl, dv = ', 1/dl,1/dv)
        dp = delg*(gascon)*(t)/(1/(dv) - 1/(dl))
        p = p + dp
    
        #thermconst['Pressure'] = psave
    return p


def tcorr(thermconst):
    """
    Function tcorr will calculate the vapor temperature, t
    and the vapor densities corresponing to the input
    pressure p
    """
    p = thermconst['Pressure']
    t = tsat(thermconst)
    gascon = thermconst['NBSgascon']
    #print ('tsat = ',t)
    thermconst['Temperature'] = t
    delg = 1
    if (t > 0):
        while (delg > 1e-4):
            thermconst['NBSdelg'] = delg
            p = ps(thermconst)
            corr(thermconst)
            delg = thermconst['NBSdelg']
            dv = thermconst['NBSdv']
            dl = thermconst['NBSdl']
            dp = delg*(gascon)*(t)/(1/(dv) - 1/(dl))
            t = (t)*(1 - dp/tdpsdt(thermconst))

    return t

def setData(thermconst):
    """
    Set all fixed parameters in the NBS EOS
    """



    #/* set the values in the common block aconst
    #        aconst[0] = wm
    #        aconst[1] = gascon
    #        aconst[2] = tz
    #        aconst[3] = aa
    #        aconst[4] = zb
    #        aconst[5] = dzb
    #        aconst[6] = yb
    #        aconst[7] = uref
    #        aconst[8] = sref
    #*/
                        



    thermconst['NBSWM'] = 18.0152
    thermconst['NBSgascon'] = 0.461522
    thermconst['NBSTz'] = 647.073
    thermconst['NBSaa'] = 1.0
    thermconst['NBSz'] = 0
    thermconst['NBSdz'] = 0
    thermconst['NBSy'] = 0
    thermconst['NBSnc'] = 35
    thermconst['NBSuRef'] = -4328.455039
    thermconst['NBSsRef'] = 7.6180802

    #/*

    #        set the values of the array g

    #*/

    g = []

    g.append(-.53062968529023e3)
    g.append(.22744901424408e4)
    g.append(.78779333020687e3)
    g.append(-.69830527374994e2)
    g.append(.17863832875422e5)
    g.append(-.39514731563338e5)
    g.append(.33803884280753e5)
    g.append(-.13855050202703e5)
    g.append(-.25637436613260e6)
    g.append(.48212575981415e6)
    g.append(-.34183016969660e6)
    g.append(.12223156417448e6)
    g.append(.11797433655832e7)
    g.append(-.21734810110373e7)
    g.append(.10829952168620e7)
    g.append(-.25441998064049e6)
    g.append(-.31377774947767e7)
    g.append(.52911910757704e7)
    g.append(-.13802577177877e7)
    g.append(-.25109914369001e6)
    g.append(.46561826115608e7)
    g.append(-.72752773275387e7)
    g.append(.41774246148294e6)
    g.append(.14016358244614e7)
    g.append(-.31555231392127e7)
    g.append(.47929666384584e7)
    g.append(.40912664781209e6)
    g.append(-.13626369388386e7)
    g.append(.69625220862664e6)
    g.append(-.10834900096447e7)
    g.append(-.22722827401688e6)
    g.append(.38365486000660e6)
    g.append(.68833257944332e4)
    g.append(.21757245522644e5)
    g.append(-.26627944829770e4)
    g.append(-.70730418082074e5)
    g.append(-.225e0)
    g.append(-1.68e0)
    g.append(.055e0)
    g.append(-93.0e0)

    thermconst['NBSg'] = g

    #/*
        
    #        set the values of the array ii

    #*/
    ii= []

    ii.append(0)
    ii.append(0)
    ii.append(0)
    ii.append(0)
    ii.append(1)
    ii.append(1)
    ii.append(1)
    ii.append(1)
    ii.append(2)
    ii.append(2)
    ii.append(2)
    ii.append(2)
    ii.append(3)
    ii.append(3)
    ii.append(3)
    ii.append(3)
    ii.append(4)
    ii.append(4)
    ii.append(4)
    ii.append(4)
    ii.append(5)
    ii.append(5)
    ii.append(5)
    ii.append(5)
    ii.append(6)
    ii.append(6)
    ii.append(6)
    ii.append(6)
    ii.append(8)
    ii.append(8)
    ii.append(8)
    ii.append(8)
    ii.append(2)
    ii.append(2)
    ii.append(0)
    ii.append(4)
    ii.append(2)
    ii.append(2)
    ii.append(2)
    ii.append(4)

    thermconst['NBSii'] = ii

    #/*

    #        set the jj array

    #*/
    jj = []

    jj.append(2)
    jj.append(3)
    jj.append(5)
    jj.append(7)
    jj.append(2)
    jj.append(3)
    jj.append(5)
    jj.append(7)
    jj.append(2)
    jj.append(3)
    jj.append(5)
    jj.append(7)
    jj.append(2)
    jj.append(3)
    jj.append(5)
    jj.append(7)
    jj.append(2)
    jj.append(3)
    jj.append(5)
    jj.append(7)
    jj.append(2)
    jj.append(3)
    jj.append(5)
    jj.append(7)
    jj.append(2)
    jj.append(3)
    jj.append(5)
    jj.append(7)
    jj.append(2)
    jj.append(3)
    jj.append(5)
    jj.append(7)
    jj.append(1)
    jj.append(4)
    jj.append(4)
    jj.append(4)
    jj.append(0)
    jj.append(2)
    jj.append(0)
    jj.append(0)

    thermconst['NBSjj'] = jj

    #/*
    #        set the values in the common block ellcon

    #                ellcon[0] = g1
    #                ellcon[1] = g2
    #                ellcon[2] = gf
    #                ellcon[3] = b1
    #                ellcon[4] = b2
    #                ellcon[5] = b1t
    #                ellcon[6] = b2t
    #                ellcon[7] = b1tt
    #                ellcon[8] = b2tt
    #*/

    thermconst['NBSg1'] =  11.0
    thermconst['NBSg2'] = 44.333333333333
    thermconst['NBSgf'] = 3.5


    #/*
        
    #        set the values of the array bp

    #*/

    bp = []

    bp.append(.7478629e0)
    bp.append(-.3540782e0)
    bp.append(0.0)
    bp.append(0.0)
    bp.append(.7159876e-2)
    bp.append(0.e0)
    bp.append(-.3528426e-2)
    bp.append(0.0)
    bp.append(0.0)
    bp.append(0.0)

    thermconst['NBSbp'] = bp

    #/*

    #        set the values of the array bq

    #*/

    bq = []

    bq.append(1.1278334e0)
    bq.append(0)
    bq.append(-.5944001e0)
    bq.append(-5.010996e0)
    bq.append(0.e0)
    bq.append(.63684256e0)
    bq.append(0)
    bq.append(0)
    bq.append(0)
    bq.append(0)

    thermconst['NBSbq'] = bq

    #/*

    #        set the values in the array atz

    #*/
    atz = []
    atz.append(64.0e1)
    atz.append(64.0e1)
    atz.append(641.6e0)
    atz.append(27.e1)

    thermconst['NBSatz'] = atz


    #/*
    #        set the values in the array adz
    #*/

    adz = []

    adz.append(0.319)
    adz.append(0.319)
    adz.append(0.319)
    adz.append(1.55e0)

    thermconst['NBSadz'] = adz

    #/*
    #        set the values in the array aat
    #*/

    aat = []

    aat.append(2.0e4)
    aat.append(2.0e4)
    aat.append(4.e4)
    aat.append(25.e0)
    thermconst['NBSaat'] = aat

    #/*
    #        set the values in the array aad
    #*/

    aad = []
    aad.append(34.e0)
    aad.append(4.e1)
    aad.append(3.e1)
    aad.append(1.05e3)

    thermconst['NBSaad'] = aad
    return

def  units(thermconst):
    """
    Function units will set the unit conversion factors
    for the NBS functions
    """

    unitsys = thermconst['NBSUnitSys']

    if unitsys == 1:


        #P --> Mpa
        #T --> K
        #D --> kg/m^3
        #H --> kJ/kg
        thermconst['NBSfp'] = 1
        thermconst['NBSft'] = 1
        thermconst['NBSfd'] = 1e-3
        thermconst['NBSfh'] = 1

    elif unitsys == 2:
        #P --> Mpa
        #T --> deg C
        #D --> kg/m^3
        #H --> kJ/kg

        thermconst['NBSfp'] = 1
        thermconst['NBSft'] = 1
        thermconst['NBSfd'] = 1e-3
        thermconst['NBSfh'] = 1

    elif unitsys == 3:
        #P --> psia
        #T --> deg F
        #D --> lb/ft^3
        #H --> btu/lb

        thermconst['NBSfp'] = 145.038
        thermconst['NBSft'] = 5.0/9.0
        thermconst['NBSfd'] = 0.016018
        thermconst['NBSfh'] = 0.4299226

    elif unitsys == 4:
        #P --> Pa
        #T --> K
        #D --> kg/m^3
        #H --> J/kg

        thermconst['NBSfp'] = 1e6
        thermconst['NBSft'] = 1
        thermconst['NBSfd'] = 1e-3
        thermconst['NBSfh'] = 1e3

    else:
         print("Invalid unit system for NBS")

    return

def convertTemperature(thermconst):
    """
    This will convert the temperature based 
    on the unit system
    """

    unitsys = thermconst['NBSUnitSys']
    convertflag = thermconst['NBSConvertFlag']

    if convertflag == "toInternal":
        if unitsys == 2:
            #
            # convert deg C to K
            thermconst['Temperature'] = thermconst['Temperature'] + 273.15
        elif unitsys == 3:
            #
            # convert deg F to K
            thermconst['Temperature'] = 5*(thermconst['Temperature'] + 459.67)/9
    elif convertflag == "toExternal":
        if unitsys == 2:
            #
            # convert K to deg C
            thermconst['Temperature'] = thermconst['Temperature'] - 273.15
        elif unitsys == 3:
            #
            # convert deg F to K
            thermconst['Temperature'] = (9*thermconst['Temperature']/5) - 459.67

    return

def calcValues(thermconst):
    """
    This function will set the values of
    u,h,s,cv,cp given T and rho
    """
    t = thermconst['Temperature']
    rho = thermconst['NBSd']

    gascon = thermconst['NBSgascon']

    rt = t*gascon
    bb(thermconst)
    qq(thermconst)
    zdum = base(thermconst)
    p = rt*rho*thermconst['NBSz'] + thermconst['NBSq']
    therm(thermconst)
    thermconst['Pressure'] = p
    thermconst['SpVol'] = 1/rho
    thermconst['SpEnthalpy'] = thermconst['NBShd']*rt
    thermconst['SpEnergy'] = thermconst['NBSud']*rt
    thermconst['SpEntropy'] = gascon*thermconst['NBSsd']
    thermconst['SpHeatp'] = thermconst['NBScpd']*gascon
    thermconst['dpdt'] = thermconst['NBSdpdt']
    thermconst['NBSdpdrho'] = thermconst['NBSdpdrho']

    return

def calc2Phase(thermconsts):
    """
    Function calc2Phase will calculate the 
    thermodynamic properties in the 2-phase 
    The user has to supply one other property 
    (v, h, s, or x) based on user key '2phaseprop'
    The function assumes the user has already calculated
    dl and dv (rhof and rhog) and they are stored in the 
    dictionary
    """

    rhof = thermconsts['NBSdl']
    thermconsts['NBSd'] = rhof
    calcValues(thermconsts)
    hf = thermconsts['SpEnthalpy']
    sf = thermconsts['SpEntropy']
    uf = thermconsts['SpEnergy']
    vf = 1/rhof
    rhog = thermconsts['NBSdv']
    thermconsts['NBSd'] = rhog
    calcValues(thermconsts)
    hg = thermconsts['SpEnthalpy']
    sg = thermconsts['SpEntropy']
    ug = thermconsts['SpEnergy']
    vg = 1/rhog
    hfg = hg - hf
    sfg = sg - sf
    vfg = vg - vf
    ufg = ug - uf

    if thermconsts['2phaseprop'] == 'h':
        hactual = thermconsts['hActual']
        #print (hactual, hfg)
        x = (hactual - hf)/hfg
        v = vf + x * (vg - vf)
        thermconsts['SpVol'] = v
        thermconsts['NBSrho'] = 1/v
        thermconsts['SpEnergy'] = uf + x * ufg
        thermconsts['SpEntropy'] = sf + x * sfg
        thermconsts['Quality'] = x
    elif thermconsts['2phaseprop'] == 's':
        sactual = thermconsts['sActual']
        x = ((sactual - sf)/sfg)
        v = vf + x * (vg - vf)
        thermconsts['SpVol'] = v
        thermconsts['NBSrho'] = 1/v
        thermconsts['SpEnthalpy'] = hf  + x*hfg
        thermconsts['SpEnergy'] = uf + x * ufg
        thermconsts['Quality'] = x
    elif thermconsts['2phaseprop'] == 'v':
        vactual = thermconsts['vActual']
        #print ('vactual =',vactual,vg,vf)
        x = 1 + ((vactual - vg)/(vg - vf))
        thermconsts['SpEnthalpy'] = hf  + x * hfg
        thermconsts['SpEnergy'] = uf + x * ufg
        thermconsts['SpEntropy'] = sf + x * sfg
        thermconsts['SpVol'] = vactual
        thermconsts['NBSrho'] = 1/vactual
        thermconsts['Quality'] = x
    elif thermconsts['2phaseprop'] == 'x':
        x = thermconsts['Quality']
        #print (vf, vg)
        v = vf + x * (vg - vf)
        thermconsts['SpVol'] = v
        thermconsts['NBSrho'] = 1/v
        thermconsts['SpEnthalpy'] = hf  + x*hfg
        thermconsts['SpEnergy'] = uf + x * ufg
        thermconsts['SpEntropy'] = sf + x * sfg
    return

def convertUnits(thermconst):
    """
    This will convert the temperature based 
    on the unit system
    """

    unitsys = thermconst['NBSUnitSys']
    convertflag = thermconst['NBSConvertFlag']

    if convertflag == "toInternal":
        thermconst['Pressure'] = thermconst['Pressure']/thermconst['NBSfp']
        thermconst['NBSrho'] = thermconst['NBSrho']*thermconst['NBSfd']
        thermconst['SpVol'] = thermconst['SpVol']/thermconst['NBSfd']
        thermconst['SpEnthalpy'] = thermconst['SpEnthalpy']/thermconst['NBSfh']
        thermconst['SpEnergy'] = thermconst['SpEnergy']/thermconst['NBSfh']
        thermconst['SpEntropy'] = thermconst['SpEntropy']/(thermconst['NBSfh']*thermconst['NBSft'])
        thermconst['SpHeatp'] = thermconst['SpHeatp']/(thermconst['NBSfh']*thermconst['NBSft'])
        thermconst['dpdt'] = thermconst['dpdt']/(thermconst['NBSfp']*thermconst['NBSft'])
        thermconst['NBSdpdrho'] = thermconst['NBSdpdrho']/(thermconst['NBSfd']*thermconst['NBSfp'])

        if unitsys == 2:
            #
            # convert deg C to K
            thermconst['Temperature'] = thermconst['Temperature'] + 273.15
        elif unitsys == 3:
            #
            # convert deg F to K
            thermconst['Temperature'] = 5*(thermconst['Temperature'] + 459.67)/9
    elif convertflag == "toExternal":
        thermconst['Pressure'] = thermconst['Pressure']*thermconst['NBSfp']
        thermconst['NBSrho'] = thermconst['NBSrho']/thermconst['NBSfd']
        thermconst['SpVol'] = thermconst['SpVol']*thermconst['NBSfd']
        thermconst['SpEnthalpy'] = thermconst['SpEnthalpy']*thermconst['NBSfh']
        thermconst['SpEnergy'] = thermconst['SpEnergy']*thermconst['NBSfh']
        thermconst['SpEntropy'] = thermconst['SpEntropy']*thermconst['NBSfh']*thermconst['NBSft']
        thermconst['SpHeatp'] = thermconst['SpHeatp']*thermconst['NBSfh']*thermconst['NBSft']
        thermconst['dpdt'] = thermconst['dpdt']*thermconst['NBSfp']*thermconst['NBSft']
        thermconst['NBSdpdrho'] = thermconst['NBSdpdrho']*thermconst['NBSfd']*thermconst['NBSfp']
        if unitsys == 2:
            #
            # convert K to deg C
            thermconst['Temperature'] = thermconst['Temperature'] - 273.15
        elif unitsys == 3:
            #
            # convert deg F to K
            thermconst['Temperature'] = (9*thermconst['Temperature']/5) - 459.67

    return

def dts(thermconsts):
    thermconsts['sActual'] = thermconsts['SpEntropy']
    d0 = thermconsts['NBSrho']
    thermtuple = 'NBS dts',thermconsts
    #print (thermtuple, type(thermtuple)
    d = scipy.optimize.newton(dfts, d0, None, thermtuple)
    return d

def dfts(d,str, thermo):
    dsave = thermo['NBSrho']
    thermo['NBSrho'] = d
    thermconst['NBSd'] = d
    calcValues(thermconst)
    #print ('dfts, v = ', v)
    thermo['NBSrho'] = dsave
    return thermconst['SpEntropy'] - thermo['sActual']

def tpd(thermconsts):
    thermconsts['pActual'] = thermconsts['Pressure']
    p0 = thermconsts['Pressure']
    thermtuple = 'NBS tpd',thermconsts
    #print (thermtuple, type(thermtuple)
    t = scipy.optimize.newton(tfpd, t0, None, thermtuple)
    return d

def tfpd(d,str, thermo):
    calcValues(thermconst)
    #print ('tfpd, v = ', v)
    return thermconst['Pressure'] - thermo['pActual']



def tps(thermconsts):
    thermconsts['sActual'] = thermconsts['SpEntropy']
    t0 = thermconsts['Temperature']
    thermtuple = 'NBS tps',thermconsts
    #print (thermtuple, type(thermtuple)
    t = scipy.optimize.newton(tfps, t0, None, thermtuple)
    return t

def tfps(t,str, thermconst):

    p = thermconst['Pressure']
    gascon = thermconst['NBSgascon']
    thermconst['Temperature'] = t

    dgss = p/(t*0.4)
    psat = 2e4
    rt = t*gascon
    thermconst['NBSdl'] = 0
    thermconst['NBSdv'] = 0
        
    bb(thermconst)
    dgss = p/(0.4*(t))
    if (t < thermconst['NBSTz']):
        thermconst['Pressure'] = psat
        psat = pcorr(thermconst)
        thermconst['Pressure'] = p
    if (p > psat):
            dgss = thermconst['NBSdl']
        
    thermconst['NBSd'] = dgss
    rho = dfind(thermconst)
    
    thermconst['NBSd'] = rho
    therm(thermconst)
    thermconst['SpEntropy'] = gascon*thermconst['NBSsd']
    err = thermconst['SpEntropy'] - thermconst['sActual']
    #print ('P,T, rho = ',p,t,rho,err,thermconst['SpEntropy'],thermconst['sActual'])
    #print ('tfpd, v = ', v)
    return thermconst['SpEntropy'] - thermconst['sActual']

def tph(thermconsts):
    thermconsts['hActual'] = thermconsts['SpEnthalpy']
    t0 = thermconsts['Temperature']
    thermtuple = 'NBS tph',thermconsts
    t = scipy.optimize.newton(tfph, t0, None, thermtuple)
    return t

def tfph(t,str, thermconst):

    p = thermconst['Pressure']
    gascon = thermconst['NBSgascon']
    thermconst['Temperature'] = t

    dgss = p/(t*0.4)
    psat = 2e4
    rt = t*gascon
    thermconst['NBSdl'] = 0
    thermconst['NBSdv'] = 0
        
    bb(thermconst)
    dgss = p/(0.4*(t))
    if (t < thermconst['NBSTz']):
        thermconst['Pressure'] = psat
        psat = pcorr(thermconst)
        thermconst['Pressure'] = p
    if (p > psat):
            dgss = thermconst['NBSdl']
        
    thermconst['NBSd'] = dgss
    rho = dfind(thermconst)
    
    thermconst['NBSd'] = rho
    therm(thermconst)
    thermconst['SpEnthalpy'] = thermconst['NBShd']*rt
    return thermconst['SpEnthalpy'] - thermconst['hActual']

def c_trans(thermconst):
    """
    Function c_tran will calculate the
    transport properties of water.  The fluid 
    properties must be calclated first

    """

    t_in = thermconst['Temperature']
    c_p = thermconst['SpHeatp']
    rho_in = 1/thermconst['SpVol']
    dpdt = thermconst['dpdt']
    drhodp = 1/thermconst['NBSdpdrho']

    t_star = 647.27e0
    if (t_in < t_star):
        #/* set constants */
        r_star = 317.763e0
        p_star = 22.115e0
        cap_c = 3.7711e-8
        omega = 0.4678e0
        cap_a = 18.66e0
        cap_b = 1.00e0
        a1 = []
        for i in range(0,4):
            a1.append(0.0181583e0)
            a1.append(0.0177624e0)
            a1.append(0.0105287e0)
            a1.append(-0.0036744e0)

        b1 = []

        for i in range(0,6):
            b1.append([])
            for j in range(0,5):                     
                b1.append(0.501938e0)
                b1.append(0.235622e0)
                b1.append(-0.274637e0)
                b1.append(0.145831e0)
                b1.append(-0.0270448e0)

                b1.append(0.162888e0)
                b1.append(0.789393e0)
                b1.append(-0.743539e0)
                b1.append(0.263129e0)
                b1.append(-0.0253093e0)

                b1.append(-0.130356e0)
                b1.append(0.673665e0)
                b1.append(-0.959456e0)
                b1.append(0.347247e0)
                b1.append(-0.0267758e0)

                b1.append(0.907919e0)
                b1.append(1.207552e0)
                b1.append(-0.687343e0)
                b1.append(0.213486e0)
                b1.append(-0.0822904e0)

                b1.append(-0.551119e0)
                b1.append(0.0670665e0)
                b1.append(-0.497089e0)
                b1.append(0.100754e0)
                b1.append(0.0602253e0)

                b1.append(0.146543e0)
                b1.append(-0.0843370e0)
                b1.append(0.195286e0)
                b1.append(-0.032932e0)
                b1.append(-0.0202595e0)

        a2 = []
        a2.append(2.022223e0)
        a2.append(14.11166e0)
        a2.append(5.25597e0)
        a2.append(-2.01870e0)

        for i in range(0,5):
            for j in range(0,6):
                b2.append(1.3293046e0)
                b2.append(-0.40452437e0)
                b2.append(0.24409490e0)
                b2.append(0.018660751e0)
                b2.append(-0.12961068e0)
                b2.append(0.044809953e0)

                b2.append(1.7018363e0)
                b2.append(-2.2156845e0)
                b2.append(1.6511057e0)
                b2.append(-0.76736002e0)
                b2.append(0.37283344e0)
                b2.append(-0.11203160e0)

                b2.append(5.2246158e0)
                b2.append(-10.124111e0)
                b2.append(4.9874687e0)
                b2.append(-0.27297694e0)
                b2.append(-0.43083393e0)
                b2.append(0.13333849e0)

                b2.append(8.7127675e0)
                b2.append(-9.5000611e0)
                b2.append(4.3786606e0)
                b2.append(-0.91783782e0)
                b2.append(0.0e0)
                b2.append(0.0e0)

                b2.append(-1.8525999e0)
                b2.append(0.9340469e0)
                b2.append(0.0e0)
                b2.append(0.0e0)
                b2.append(0.0e0)
                b2.append(0.0e0)

    #/*----------------------
    #    evaluate viscosity mu
    #    ---------------------*/

        sum = 0
        for i in (0,4):
            sum = sum + a1[i]*((t_star/t_in)**i)

        eta_0 = (1.e-6)*math.sqrt((t_in/t_star))/sum
        sum = 0
        for i in range(0,6):
            for j in range(0,5):
                sum = sum + b1[i][j]*(((t_star/t_in)-1)**i)*(((rho_in/r_star)-1)**j)

        #/* set viscosoty according to eq'n C.1 */
        mu = eta_0*math.exp(rho_in*sum/r_star)
        thermconst['NBSViscosity'] = mu

        #/*-----------------------------------
        #evaluate thermal conductivity k
        #-----------------------------------*/


        #/* calculate chi_t */

        chi_t = rho_in*drhodp*p_star/(r_star*r_star)

        #/* calculate  lambda  */

        val = math.exp(-cap_a*(((t_star/t_in) - 1)**2)-cap_b*(((rho_in/r_star)-1)**4))
        xlambda = (cap_c/(mu))*((t_in*r_star/(t_star*rho_in))**2)
        xlambda = xlambda*((t_star*dpdt/p_star)**2)*(chi_t**omega)
        xlambda = xlambda*math.sqrt(rho_in/r_star)*val
      
        lambda_0 = math.sqrt(t_in/t_star)
        
        sum = 0
        for i in range(0,4):
            sum = sum + a2[i]*((t_star/t_in)**i)

        lambda_0 = lambda_0/sum
        sum = 0
        for i in range(0,5):
            for j in range(0,6):
                sum = sum + b2[i][j]*((t_star/t_in - 1)**i)*((rho_in/r_star-1)**j)
        sum = math.exp(rho_in*sum/r_star)
        k = lambda_0*sum + xlambda

        thermconst['NBSThermalConductivity'] = k


        #/*---------------------------
        #evaluate Prandtl number Pr
        #---------------------------*/

        #/* use equation C.3 from appendix C  */

        Pr = 1000*(mu)*c_p/(k)
        thermconst['NBSPrandtl'] = Pr

        #    evaluate surface tension sigma
        
   
        # use equation C.5 from NBS appendix C  

        sigma = 0.2358*(((647.15 - t_in)/647.15**1.256)*(1+(-0.625)*((647.15-t_in)/647.15)))
        thermconst['NBSSurfaceTension'] = sigma

        #/* evaluate beta */

        beta = rho_in*dvdt/1000
        thermconst['NBSCompressibility'] = beta

        #/* evaluate isothermal compressibility, kappa */
   
        kappa = rho_in*dvdt/(1000000000*dpdt)
        thermconst['NBSIsothermalCompressibility'] = kappa
        
        
    return