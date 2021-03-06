import IdealGas as ideal
import RKGas as RK
import LKGas as LK
import StanfordEqns as staneqs
import NBS as nbs
import ThermoConst as thermc
import math
import sys
import getopt
import numpy as np
import matplotlib.pyplot as plt


class ThermoProps(object):
    """A class to calculate thermo props for a specified EOS"""

    def __init__(self):
        self.Thermo = {}
        #
        # set required keys
        self.Thermo['Pressure'] = 0
        self.Thermo['Temperature'] = 0
        self.Thermo['SpVol'] = 0
        self.Thermo['SpEnthalpy'] = 0
        self.Thermo['SpEnergy'] = 0
        self.Thermo['SpEntropy'] = 0
        self.Thermo['Quality'] = 0
        self.Thermo['Region'] = ""
        self.Thermo['Density'] = 0 
        self.Thermo['chartType'] = 'none'
        self.Thermo['chartIso'] = 'none'
        self.Thermo['chartNumLines'] = 0
        self.Thermo['chartTmin'] = 0.0
        self.Thermo['chartTmax'] = 0.0
        self.Thermo['chartPmin'] = 0.0
        self.Thermo['chartPmax'] = 0.0

    #
    # define sepcial functions to get fluid properties

    def getTCrit(self):
        return self.Thermo['TCrit']

    def getPCrit(self):
        return self.Thermo['PCrit']

    def getvCrit(self):
        return self.Thermo['vCrit']

    def getFluidName(self):
        return self.Thermo['FluidName']

    def getFluidFormula(self):
        return self.Thermo['FluidFormula']

    def getMW(self):
        return self.Thermo['MW']

    def getStanfordT0(self):
        if self.Thermo['Stanfordeos'].lower() == "true":
            return self.Thermo['StanfordT0']
        else:
            print ("Stanford props not defined")
            return -1

    #
    # Define get/set functions for thermo properties
    # T, p, v, u, h, s, x

    def getPressure(self):
        return self.Thermo['Pressure']

    def setPressure(self,p):
        self.Thermo['Pressure'] = p
        return

    def getTemperature(self):
        return self.Thermo['Temperature'] 

    def setTemperature(self,t):
        self.Thermo['Temperature'] = t
        return

    def getSpecificVolume(self):
        return self.Thermo['SpVol']

    def setSpecificVolume(self,v):
        self.Thermo['SpVol'] = v
        return

    def getSpecificEnergy(self):
        return self.Thermo['SpEnergy']

    def setSpecificEnergy(self,u):
        self.Thermo['SpEnergy'] = u
        return

    def getSpecificEnthalpy(self):
        return self.Thermo['SpEnthalpy']

    def setSpecificEnthalpy(self,h):
        self.Thermo['SpEnthalpy'] = h
        return

    def getSpecificEntropy(self):
        return self.Thermo['SpEntropy']

    def setSpecificEntropy(self,s):
        self.Thermo['SpEntropy'] = s
        return
    def getQuality(self):
        return self.Thermo['Quality']

    def setQuality(self,x):
        self.Thermo['Quality'] = x
        return

    def setJobID(self,job):
        """
        set the job id for the props
        """
        self.Thermo['job'] = job
        return

    def getUnitSys(self):
        return self.Thermo['NBSUnitSys'] 

    def getRegion(self):
        return self.Thermo['Region']

    def setUnitSys(self, u):
        self.Thermo['NBSUnitSys'] = u
        return

    def setJob(self,strjob):
        jobl = strjob.lower()
        if jobl == 'tp' or jobl == 'pt':
            setJobID(1)
        elif jobl == 'tv' or jobl == 'vt':
            setJobID(2)
        elif jobl == 'ts' or jobl == 'st':
            setJobID(3)
        elif jobl == 'tx' or jobl == 'xt':
            setJobID(4)
        elif jobl == 'pv' or jobl == 'vp':
            setJobID(5)
        elif jobl == 'ph' or jobl == 'hp':
            setJobID(6)
        elif jobl == 'ps' or jobl == 'sp':
            setJobID(7)
        elif jobl == 'px' or jobl == 'xp':
            setJobID(8)
        else:
            print ('Invalid job type')
        return

    def getEOS(self):
        return self.Thermo['eos']

    def setEOS(self,eos):
        self.Thermo['eos'] = eos
        return


    def setFluid(self, FluidName, clear='True'):
        """Use this function to read fluid properties into the hash table"""
        if clear == 'True':
            self.Thermo.clear()
        Filename = FluidName + ".pfd"
        linenum = 0
        for line in open(Filename,"r"):
            # print (line)
            tokens = line.split()
            linenum = linenum + 1
            # print (tokens)
            if len(tokens) > 1:
                if tokens[0] == "#":
                    # comment line ignore
                    # print ('Comment ',line)
                    pass
                elif len(tokens) > 2:
                    if tokens[1] == "=":
                        # print ('data:   ',line)
                        # line contains data
                        if tokens[2].lower() == 'float':
                           # print (tokens[0], tokens[2], tokens[3])
                            try:
                                self.Thermo[tokens[0]] = float(tokens[3])
                            except:
                                print ('Error parsing file ',Filename," at line number ",linenum)
                                sys.exit(2)


                        elif tokens[2].lower() == 'int':
                            self.Thermo[tokens[0]] = int(tokens[3])
                        else:
                            self.Thermo[tokens[0]] = tokens[3].strip()
        return

    def setchartType(self,type):
        """
        set teh chart type for the properties.
        there are 3 types, T-s, P-h and P-v
        """
        typel =type.lower()
        if typel in ['ts','ph','pv']:
            self.Thermo['chartType'] = type
        else:
            print ('Invalid chart type')
        return

    def getchartType(self):
        """
        get the chart type from the 
        dictionary
        """
        return self.Thermo['chartType']

    def setchartIsoLines(self,iso):
        """
        This will set the iso lines to plot 
        on the graph.  Valid values are:
        for T-s chart: p,h,x,v
        for P-h chart: t,x,v,s
        for h-s chart: t,p,v,x

        since this is dependent on the chart 
        type which may or may not be set, the check will be done 
        when generating the chart
        """

        self.Thermo['chartIso'] = iso.lower()
        return

    def getchartIsoLines(self):
        return self.Thermo['chartIso']

    def setchartNumIsoLines(self,numlines):
        self.Thermo['chartNumLines'] = numlines
        return

    def getchartNumLines(self):
        return self.Thermo['chartNumLines']

    def setchartTmin(self, tmin):
        self.Thermo['chartTmin'] = tmin
        return

    def getchartTmin(self):
        return self.Thermo['chartTmin']


    def setchartTmax(self, tmax):
        self.Thermo['chartTmax'] = tmax
        return

    def getchartTmax(self):
        return self.Thermo['chartTmax']

    def setchartPmin(self, pmin):
        self.Thermo['chartPmin'] = pmin
        return

    def getchartPmin(self):
        return self.Thermo['chartPmin']


    def setchartPmax(self, pmax):
        self.Thermo['chartPmax'] = pmax
        return

    def getchartPmax(self):
        return self.Thermo['chartPmax']

    def calcProps(self):
        eos = self.Thermo['eos']
        self.Thermo['cvfunc'] = 'nothing'
        #
        # Now check if the ieal gas law is valid, if so add pointers to functions
        self.Thermo['vtp'] = 'False'
        if eos == 'Ideal':
            eosCheck = self.Thermo['Idealeos']
            if eosCheck == 'True':
                self.Thermo['pvt'] = ideal.pvt
                self.Thermo['tpv'] = ideal.tpv
                self.Thermo['vtp'] = ideal.vtp
                self.Thermo['dptdv'] = ideal.dptdv
                self.Thermo['dpvdt'] = ideal.dpvdt
                self.Thermo['dvint'] = ideal.dvint
                self.Thermo['vts'] = ideal.vts
                self.Thermo['tph'] = ideal.tph
                self.Thermo['tps'] = ideal.tps
                if self.Thermo['Stanfordeos'] == 'True':
                    self.setStanfordcv()
                self.gasProps()
      
        elif eos == "LK":
            #
            #simple fluid params
            #
            self.Thermo['LKb1simple'] = 0.1181193 # simple fluid b1
            self.Thermo['LKb2simple'] = 0.265728 # simple fluid b2
            self.Thermo['LKb3simple'] = 0.154790 # simple fluid b3
            self.Thermo['LKb4simple'] = 0.030323 # simple fluid b4
            self.Thermo['LKc1simple'] = 0.0236744 # simple fluid c1
            self.Thermo['LKc2simple'] = 0.0186984 # simple fluid c2
            self.Thermo['LKc3simple'] = 0.0      # simple fluid c3
            self.Thermo['LKd1simple'] = 0.155488e-4 # simple fluid d1
            self.Thermo['LKd2simple'] = 0.623689e-4 # simple fluid d2
            self.Thermo['LKc4simple'] = 0.042724 # simple fluid c4
            self.Thermo['LKbetasimple'] = 0.65392 # simple fluid beta
            self.Thermo['LKgammasimple'] = 0.060167 # simple fluid gamma
            #
            # reference fluid (octane) parameters
            #
            self.Thermo['LKb1ref'] = 0.2026579 # ref fluid b1
            self.Thermo['LKb2ref'] = 0.331511 # ref fluid b2
            self.Thermo['LKb3ref'] = 0.027655 # ref fluid b3
            self.Thermo['LKb4ref'] = 0.203488 # ref fluid b4
            self.Thermo['LKc1ref'] = 0.0313385 # ref fluid c1
            self.Thermo['LKc2ref'] = 0.0503618 # ref fluid c2
            self.Thermo['LKc3ref'] = 0.16901      # ref fluid c3
            self.Thermo['LKd1ref'] = 0.48736e-4 # ref fluid d1
            self.Thermo['LKd2ref'] = 0.0740336e-4 # ref fluid d2
            self.Thermo['LKc4ref'] = 0.041577 # ref fluid c4
            self.Thermo['LKbetaref'] = 1.226 # ref fluid beta
            self.Thermo['LKgammaref'] = 0.03754 # ref fluid gamma
            self.Thermo['LKRgasref'] = 72.78171 # ref Rgas
            self.Thermo['LKTcritref'] = 568.8 # ref fluid Tcrit
            self.Thermo['LKPcritref'] = 2.49e6 # ref fluid Pcrit
            self.Thermo['LKvcritref'] = 0.492 # ref fluid vcrit
            self.Thermo['LKMWref'] = 114.232 # ref fluid MW
            self.Thermo['LKomegaref'] = 0.3978 # ref fluid accentric
            self.Thermo['pvt'] = LK.pvt
            #self.Thermo['tpv'] = ideal.tpv
            self.Thermo['vtp'] = LK.vtp
            self.Thermo['dvint'] = LK.dvint
            
            self.Thermo['vts'] = LK.vts
            self.Thermo['tpv'] = LK.tpv
            self.Thermo['tph'] = LK.tph
            self.Thermo['tps'] = LK.tps
            if self.Thermo['Stanfordeos'] == 'True':
                self.setStanfordcv()
            eosCheck = self.Thermo['LKeos']
            if eosCheck == 'True':
                self.gasProps()
            else:
                print ('Invalid eos')
        elif eos == "RK":
            Tcrit = self.Thermo['TCrit']
            Pcrit = self.Thermo['PCrit']
            Rgas = thermc.R_BAR/self.Thermo['MW']
            self.Thermo['RKa'] = 0.427480*Rgas*Rgas*(Tcrit**2.5)/(Pcrit)
            self.Thermo['RKb'] = 0.086640*Rgas*Tcrit/(Pcrit)
            self.Thermo['pvt'] = RK.pvt
            #self.Thermo['tpv'] = ideal.tpv
            self.Thermo['vtp'] = RK.vtp
            self.Thermo['dvint'] = RK.dvint
            self.Thermo['vts'] = RK.vts
            self.Thermo['tpv'] = RK.tpv
            self.Thermo['tph'] = RK.tph
            self.Thermo['tps'] = RK.tps
            if self.Thermo['Stanfordeos'] == 'True':
                self.setStanfordcv()
            eosCheck = self.Thermo['RKeos']
            if eosCheck == 'True':
                self.gasProps()
            else:
                print ('Invalid eos')

        elif eos == 'Stanford':
            eosCheck = self.Thermo['Stanfordeos']
            if eosCheck == 'True':
                self.StanProps()
            else:
                print ('Invalid eos')
            
        elif eos == "NBS":
            eosCheck = self.Thermo['NBSeos']
            if eosCheck == 'True':
                self.NBS()
            else:
                print ('Invalid eos')


    def setStanfordcv(self):
        """
        Used by Ideal, Rk, Lk to set a better
        cv function if the Stanford eqn set 
        exists
        """
        cv = self.Thermo['Stanford_cv']
        if cv == 0:
            self.Thermo['cvfunc'] = staneqs.cv0
        elif cv == 1:
            self.Thermo['cvfunc'] = staneqs.cv1
        elif cv == 2:
            self.Thermo['cvfunc'] = staneqs.cv2
        elif cv == 3:
            self.Thermo['cvfunc'] = staneqs.cv3
        elif cv == 4:
            self.Thermo['cvfunc'] = staneqs.cv4
        elif cv == 5:
            self.Thermo['cvfunc'] = staneqs.cv5
        elif cv == 6:
            self.Thermo['cvfunc'] = staneqs.cv6
        else:
            pass
        return
            
            
#Function GasProps will calculate the
#thermodynamic properties of a fluid based on the
#job number.  The job number will indicate the two
#properties passed to the function.  The equation
#of state and any pertinent functions will be passed
#via the function pointer FunctionHead.  This structure
#is defined in the jrmprops.h file.  Any constants
#associated with the fluid will be in the array consts
        

#UNITS:
#        Pressure                p                       Pa
#        Temperature             t                       K
#        Specific volume         v                       m^3/kg
#        Specific Energy         u                       J/kg
#        Specific Enthalpy       h                       J/kg
#        Specific Entropy        s                       J/kg K


#INPUT:
#        INT             job     indicates properties passed to routine

#                                =       1   T,P are inputs
#                                =       2   T,v are inputs
#                                =       3   T,s are inputs
#                                =       4   T,u are inputs
#                                =       5   P,v are inputs
#                                =       6   P,h are inputs
#                                =       7   P,s are inputs
#                                =       8   P,u are inputs
#                                =       9   T,h are inputs
#                                =       10  h,s are inputs
    def gasProps(self):
        
        job = self.Thermo['job']
        self.Thermo['Region'] = 'SHV'

        if job == 1:
            # 
            # User input T,p
            if self.Thermo['vtp'] == 'False':
                print ("No vtp function")

            else:
                vtp = self.Thermo['vtp']
                v  = vtp(self.Thermo)
                if v > 0: 
                    self.Thermo['SpVol'] = v
                    if self.Thermo['cvfunc'] != 'nothing':
                        (cx, cxint, cxtint) = self.Thermo['cvfunc'](self.Thermo)
                    elif self.Thermo['cxCurveFit'] == 'True':
                        (cx, cxint, cxtint) = thermc.Cubiccv(self.Thermo)
                    (udvint, sdvint) = self.Thermo['dvint'](self.Thermo)
                    self.Thermo['SpEnergy'] = cxint +udvint + self.Thermo['xRef']
                    self.Thermo['SpEntropy'] = cxtint + sdvint + self.Thermo['sRef'];
                    # print (cxtint, sdvint, self.Thermo['sRef'])
                    self.Thermo['SpEnthalpy'] = self.Thermo['SpEnergy'] + self.Thermo['Pressure']*self.Thermo['SpVol']

        if job == 2:
            #
            # user inputs T, v
            self.Thermo['Pressure'] = self.Thermo['pvt'](self.Thermo)
            if self.Thermo['cvfunc'] != 'nothing':
                (cx, cxint, cxtint) = self.Thermo['cvfunc'](self.Thermo)
            elif self.Thermo['cxCurveFit'] == 'True':
                (cx, cxint, cxtint) = thermc.Cubiccv(self.Thermo)
            (udvint, sdvint) = self.Thermo['dvint'](self.Thermo)
            self.Thermo['SpEnergy'] = cxint +udvint + self.Thermo['xRef']
            self.Thermo['SpEntropy'] = cxtint + sdvint + self.Thermo['sRef'];
            # print (cxtint, sdvint, self.Thermo['sRef'])
            self.Thermo['SpEnthalpy'] = self.Thermo['SpEnergy'] + self.Thermo['Pressure']*self.Thermo['SpVol']

        if job == 3:
            # 
            # User input T,s
            if self.Thermo['vts'] == 'False':
                print ("No vts function")

            else:
                vts = self.Thermo['vts']
                v  = vts(self.Thermo)
                pvt = self.Thermo['pvt']
                self.Thermo['Pressure'] = pvt(self.Thermo)
                if v > 0: 
                    self.Thermo['SpVol'] = v
                    if self.Thermo['cvfunc'] != 'nothing':
                        (cx, cxint, cxtint) = self.Thermo['cvfunc'](self.Thermo)
                    elif self.Thermo['cxCurveFit'] == 'True':
                        (cx, cxint, cxtint) = thermc.Cubiccv(self.Thermo)
                    (udvint, sdvint) = self.Thermo['dvint'](self.Thermo)
                    self.Thermo['SpEnergy'] = cxint +udvint + self.Thermo['xRef']
                    #self.Thermo['SpEntropy'] = cxtint + sdvint + self.Thermo['sRef'];
                    # print (cxtint, sdvint, self.Thermo['sRef'])
                    self.Thermo['SpEnthalpy'] = self.Thermo['SpEnergy'] + self.Thermo['Pressure']*self.Thermo['SpVol']


        if job == 5:
            # 
            # User input P,v
            if self.Thermo['tpv'] == 'False':
                print ("No tpv function")

            else:
                tpv = self.Thermo['tpv']
                T  = tpv(self.Thermo)
                # print ('T = ',T)
                if T > 0: 
                    self.Thermo['Temperature'] = T
                    if self.Thermo['cvfunc'] != 'nothing':
                        (cx, cxint, cxtint) = self.Thermo['cvfunc'](self.Thermo)
                    elif self.Thermo['cxCurveFit'] == 'True':
                        (cx, cxint, cxtint) = thermc.Cubiccv(self.Thermo)
                    (udvint, sdvint) = self.Thermo['dvint'](self.Thermo)
                    self.Thermo['SpEnergy'] = cxint +udvint + self.Thermo['xRef']
                    #self.Thermo['SpEntropy'] = cxtint + sdvint + self.Thermo['sRef'];
                    # print (cxtint, sdvint, self.Thermo['sRef'])
                    self.Thermo['SpEnthalpy'] = self.Thermo['SpEnergy'] + self.Thermo['Pressure']*self.Thermo['SpVol']

        if job == 6:
            # 
            # User input P,h
            if self.Thermo['tph'] == 'False':
                print ("No tph function")

            else:
                tph = self.Thermo['tph']
                T  = tph(self.Thermo)
                # print ('T = ',T)
                if T > 0: 
                    self.Thermo['Temperature'] = T    
                    vtp = self.Thermo['vtp']
                    v = vtp(self.Thermo)
                    self.Thermo['SpVol'] = v

                    if self.Thermo['cvfunc'] != 'nothing':
                        (cx, cxint, cxtint) = self.Thermo['cvfunc'](self.Thermo)
                    elif self.Thermo['cxCurveFit'] == 'True':
                        (cx, cxint, cxtint) = thermc.Cubiccv(self.Thermo)
                    (udvint, sdvint) = self.Thermo['dvint'](self.Thermo)
                    self.Thermo['SpEnergy'] = cxint +udvint + self.Thermo['xRef']
                    #self.Thermo['SpEntropy'] = cxtint + sdvint + self.Thermo['sRef'];
                    # print (cxtint, sdvint, self.Thermo['sRef'])
                    self.Thermo['SpEnthalpy'] = self.Thermo['SpEnergy'] + self.Thermo['Pressure']*self.Thermo['SpVol']
        if job == 7:
            # 
            # User input P,s
            if self.Thermo['tph'] == 'False':
                print ("No tps function")

            else:
                tps = self.Thermo['tps']
                T  = tps(self.Thermo)
                # print ('T = ',T)
                if T > 0: 
                    self.Thermo['Temperature'] = T    
                    vtp = self.Thermo['vtp']
                    v = vtp(self.Thermo)
                    self.Thermo['SpVol'] = v

                    if self.Thermo['cvfunc'] != 'nothing':
                        (cx, cxint, cxtint) = self.Thermo['cvfunc'](self.Thermo)
                    elif self.Thermo['cxCurveFit'] == 'True':
                        (cx, cxint, cxtint) = thermc.Cubiccv(self.Thermo)
                    (udvint, sdvint) = self.Thermo['dvint'](self.Thermo)
                    self.Thermo['SpEnergy'] = cxint +udvint + self.Thermo['xRef']
                    #self.Thermo['SpEntropy'] = cxtint + sdvint + self.Thermo['sRef'];
                    # print (cxtint, sdvint, self.Thermo['sRef'])
                    self.Thermo['SpEnthalpy'] = self.Thermo['SpEnergy'] + self.Thermo['Pressure']*self.Thermo['SpVol']

    def StanProps(self):
        staneqs.Stanford(self.Thermo)

    def NBS(self):
        nbs.NBSprops(self.Thermo)


    def makeChart(self):
        """
        This function will create a chart (T-s, P-v, P-h, h-s)
        depending on user input.  The chart will show the 
        saturation dome and a number of iso lines chosen by the user
        (one set of iso lines at this time).  the code will start at 
        5K above T0 and stop 5K below Tcrit (for T-s, same for others)

        The function is only valid for Stanford and NBS
        """

        eos = self.Thermo['eos']
        useLogx = False
        useLogy = False

        if eos in ['NBS', 'Stanford']:
            print ('OK')
            chartData = []
            if self.Thermo['chartType'] == 'ts':
                #
                # make a T-s chart
               self.tsChart(chartData)
               yLabel = 'T [K]'
               xLabel = 's [J/kg K]'
            elif self.Thermo['chartType'] == 'ph':
                self.phChart(chartData)
                useLogy = True
                xLabel = 'h [J/kg]'
                yLabel = 'P [Pa]'
            else:
                #
                # right now pass
                pass

            xvals = []
            yvals = []
            zvals = chartData[0]
            n = 0
            dxmax = 0.0
            dymax = 0.0

            for i in range(2,len(chartData)):
                (xvals,yvals) = chartData[i]
                xmin = min(xvals)
                xmax = max(xvals)
                ymin = min(yvals)
                ymax = max(yvals)
                dx = math.fabs(xmax - xmin)
                dy = math.fabs(ymax - ymin)
                if dx > dxmax:
                    dxmax = dx
                if dy > dymax:
                    dymax = dy

            colors = ['black','blue','green','red','yellow']
            colorindex = 0
            numlines = self.getchartNumLines()
            for i in range(1,len(chartData)):
                (xvals,yvals) = chartData[i]
                #line = plt.Line2D(xvals,yvals,label="line "+str(n))
                line =  plt.plot(xvals,yvals)
                plt.setp(line,color=colors[colorindex])
                #
                # add text label to data
                if (i > 1) and len(xvals) > 3:
                    if (i-2) % numlines == 0:
                        colorindex = colorindex + 1
                    plt.setp(line,color=colors[colorindex])
                    midindex = int(len(xvals)/2)
                    x1 = xvals[midindex]
                    y1 = yvals[midindex]
                    x2 = xvals[midindex + 1]
                    y2 = yvals[midindex + 1]
                    dx = (x2 - x1)/dxmax
                    dy = (y2 - y1)/dymax
                    if math.fabs(dx) <= 1e-2:
                        #
                        # approx zero
                        rotate = 90
                    elif math.fabs(dy) < 1e-2:
                        #
                        # approx 0
                        rotate = 0
                    else:
                        rotate = math.degrees(math.atan(dy/dx))

                    #print (zvals[i-2],'rotate = ',rotate, 180-math.fabs(rotate))
                    if rotate < 0:
                        rotate = 180-math.fabs(rotate)
                    plt.annotate(zvals[i-2],xy=(xvals[midindex],yvals[midindex]), xytext=(xvals[midindex],yvals[midindex]),xycoords='data',ha="center", va="center", rotation=(rotate))

                #plt.contour(xvals,yvals,zvals,zvals)
            if useLogx == True:
                plt.xscale('log')
            if useLogy == True:
                plt.yscale('log')

            plt.xlabel(xLabel)
            plt.ylabel(yLabel)

            
            plt.show()

            
        else:
            print ('Invalid eos for charts')



            return

    def tsChart(self, chartData):
        """
        this function will generate the data for a 
        T-s diagram for the specified fluid
        and store the data in the chartData list
        """

        print ("Creating saturation dome")

        #
        # get the bounds values for the iso lines
        Tcrit = self.getTCrit()
        T0 = self.getStanfordT0()
        alpha = 0.99

        Tmax = self.getchartTmax()
        if Tmax == 0:
            Tmax = Tcrit

        Tmin = self.getchartTmin()

        if Tmin < T0:
            Tmin = T0

        numlines = self.Thermo['chartNumLines']
        numSatPoints = 80
        
        #
        # set the number of sat dome points 
        # and deltaT for saturation calcs

        satDomes = []
        satDomeT = []
        satDomex = []

        deltaT = (alpha*Tcrit - T0)/(numSatPoints + 1)
        Tcurr = T0
        #
        # start by going up the sat liq side of the dome
        self.setQuality(0.0)
        self.setJobID(4)
        for i in range(numSatPoints):
            Tcurr = Tcurr + deltaT
            self.setTemperature(Tcurr)
            self.calcProps()
            satDomeT.append(self.getTemperature())
            satDomes.append(self.getSpecificEntropy())
            #print (self.getRegion(), Tcurr, self.getQuality(), self.getSpecificEntropy())
#            satDomex.append(0.0)

        #
        # now add the sat vap side of the dome
        self.setQuality(1.0)
        Tcurr = alpha*Tcrit
        for i in range(numSatPoints):
            Tcurr = Tcurr - deltaT
            self.setTemperature(Tcurr)
            self.calcProps()
            satDomeT.append(self.getTemperature())
            satDomes.append(self.getSpecificEntropy())
            #print (self.getRegion(), Tcurr, self.getQuality(), self.getSpecificEntropy())
 #           satDomex.append(1.0)

        # make room for the labels
        chartData.append([])
        chartData.append((satDomes, satDomeT))

        print ("    Finished saturation dome")

        #
        # get the chart iso line variable
        # and set the number of points

        chartIso = self.Thermo['chartIso']
        numPoints = 30
        deltaT = (Tcrit - T0)/(numPoints + 1)

        isoT = []
        isos = []
        isoz = []
        
        for isovar in chartIso:

            if isovar == 'x':
                print ("Creating iso-quality lines")
                #
                # iso lines for quality
                dx = 1.0/(numlines + 1.0)
                x = 0.0
                #
                # set delta T based on T0 and Tcrit
                deltaT = ((Tcrit-(0.1*Tcrit)) - T0)/(numPoints + 1)

                numSets = len(isoT)
                for n in range(numlines):
                    setNum = n + numSets
                    isoT.append([])
                    isos.append([])
                    Tcurr = T0
                    x = x + dx
                    self.setQuality(x)
                    self.setJobID(4)
                    isoz.append('x = %.2f' % (x))
                    for m in range(numPoints):
                        Tcurr = Tcurr + deltaT
                        self.setTemperature(Tcurr)
                        self.calcProps()
                        isoT[setNum].append(self.getTemperature())
                        isos[setNum].append(self.getSpecificEntropy())
                        
                        #print (self.getRegion(), Tcurr, self.getQuality(), self.getSpecificEntropy())
                print ("    Finished iso-quality lines")

            elif isovar == 'p':
                #
                # iso lines for pressure
                # 

                # note that the user can specify 
                # a max > Tcrit, so we need to 
                # have a check on that for the iso lines



                #
                # get pcrit and Psat(T0)
                print ("Creating isobaric lines")
                pcrit = self.getPCrit()
                self.setTemperature(self.getStanfordT0())
                self.setQuality(0.0)
                self.setJobID(4)
                self.calcProps()
                psat0 = self.getPressure()

                #
                # get Pmin and Pmax

                Pmax = self.getchartPmax()
                if Pmax == 0 or Pmax > pcrit:
                    Pmax = pcrit

                Pmin = self.getchartPmin()
                if Pmin == 0 or Pmin > pcrit: 
                    Pmin = psat0

                #
                # now check the ratio of Pmax/Pmin
                # if it's large use ln to move thru
                # values instead of linear
                
                if Pmax/Pmin > 10:
                    dp = math.log(Pmax/Pmin)/(numlines + 1.0)
                    useLn = True
                else:
                    dp = (Pmax - Pmin)/(numlines + 1.0)
                    useLn = False
                numSets = len(isoT)
                psat = Pmin

                for n in range(numlines):
                    if useLn:
                        psat = math.exp(math.log(psat) + dp)
                    else:
                        psat = psat + dp
                    setNum = n + numSets
                    #
                    # create a new set of data
                    isoT.append([])
                    isos.append([])
                    #
                    # add the saturation points for the 
                    # isobaric lines

                    self.setPressure(psat)
                    self.setQuality(0.0)
                    self.setJobID(8)
                    self.calcProps()
                    tsat = self.getTemperature()
                    isoT[setNum].append(tsat)
                    isos[setNum].append(self.getSpecificEntropy())
                    self.setQuality(1.0)
                    self.setJobID(8)
                    self.calcProps()
                    isoT[setNum].append(tsat)
                    isos[setNum].append(self.getSpecificEntropy())
                    self.setJobID(1)
                    Tcurr = tsat

                    deltaT = (Tmax - Tcurr)/(numPoints + 1)

                    isoz.append('P = %1.2e' % (psat))
                    for m in range(numPoints):
                        Tcurr = Tcurr + deltaT
                        self.setTemperature(Tcurr)
                        self.calcProps()
                        isoT[setNum].append(self.getTemperature())
                        isos[setNum].append(self.getSpecificEntropy())
                print ("    Finished isobaric lines")

            elif isovar == 'v':
                #
                # iso lines for specific volume
                # 

                # note that the user can specify 
                # a max > Tcrit, so we need to 
                # have a check on that for the iso lines

                #
                # get pcrit and Psat(T0)
                print ("Creating iso-volume lines")
                vcrit = self.getvCrit()

                #
                # set vmax = vg(Tmin,x = 1)
                self.setTemperature(Tmin)
                self.setQuality(1.0)
                self.setJobID(4)
                self.calcProps()
                vmax = self.getSpecificVolume()

                #
                # set vmin = v(Tmin, x = 0)
                self.setTemperature(Tmin)
                self.setQuality(0.0)
                self.setJobID(4)
                self.calcProps()
                vmin = self.getSpecificVolume()   

                #
                # now check the ratio of vmax/vmin
                # if it's large use ln to move thru
                # values instead of linear

                if vmax/vmin > 10:
                    dv = math.log(vmax/vmin)/(numlines + 1.0)
                    useLn = True
                else:                    
                    dv = (vmax - vmin)/(numlines + 1.0)
                    useLn = False

                numSets = len(isoT)
                vcurr = vmin
                
                for n in range(numlines):
                    if useLn:
                        vcurr = math.exp(math.log(vcurr) + dv)
                    else:
                        vcurr = vcurr + dv
                    self.setJobID(2)
                    self.setSpecificVolume(vcurr)
                    setNum = n + numSets
                    #
                    # create a new set of data
                    isoT.append([])
                    isos.append([])
                    #
                    # add the saturation points for the 
                    # iso volume lines

                    Tcurr = Tmin

                    deltaT = (Tmax - Tmin)/(numPoints + 1)
                    isoz.append('v = %1.2e' % (vcurr))
                    for m in range(numPoints):
                        Tcurr = Tcurr + deltaT
                        self.setTemperature(Tcurr)
                        self.calcProps()
                        isoT[setNum].append(self.getTemperature())
                        isos[setNum].append(self.getSpecificEntropy())
                        #print (self.getRegion(), Tcurr, vcurr, self.getQuality(), self.getSpecificEntropy())

                print ("    Finished iso-volume lines")
            elif isovar == 'h':
                #
                # iso lines for specific enthalpy
                # 

                # note that the user can specify 
                # a max > Tcrit, so we need to 
                # have a check on that for the iso lines

                #
                print ("Creating iso-enthalpy lines")

                #
                # set hmax = hg(Tmax,x = 1)
                if Tmax < Tcrit:
                    self.setTemperature(Tmax)
                else:
                    self.setTemperature(0.95*Tcrit)
                self.setQuality(1.0)
                self.setJobID(4)
                self.calcProps()
                hmax = self.getSpecificEnthalpy()

                #
                # set hmin = hf(Tmax, x = 0)
                self.setTemperature(Tmax)
                self.setQuality(0.0)
                self.setJobID(4)
                self.calcProps()
                hmin = self.getSpecificEnthalpy()

                #
                # now check the ratio of hmax/hmin
                # if it's large use ln to move thru
                # values instead of linear

                if hmax/hmin > 10:
                    dh = math.log(hmax/hmin)/(numlines + 1.0)
                    useLn = True
                else:                    
                    dh = (hmax - hmin)/(numlines + 1.0)
                    useLn = False

                numSets = len(isoT)
                hcurr = hmin
                
                for n in range(numlines):
                    if useLn:
                        hcurr = math.exp(math.log(hcurr) + dh)
                    else:
                        hcurr = hcurr + dh
                    self.setJobID(9)
                    self.setSpecificEnthalpy(hcurr)
                    setNum = n + numSets
                    #
                    # create a new set of data
                    isoT.append([])
                    isos.append([])
                    #
                    # add the saturation points for the 
                    # iso volume lines

                    Tcurr = Tmin

                    deltaT = (Tmax - Tmin)/(numPoints + 1)

                    isoz.append('h = %1.2e'% (hcurr))
                    for m in range(numPoints):
                        Tcurr = Tcurr + deltaT
                        self.setTemperature(Tcurr)
                        self.calcProps()
                        isoT[setNum].append(self.getTemperature())
                        isos[setNum].append(self.getSpecificEntropy())
                        #print (self.getRegion(), Tcurr, hcurr, self.getQuality(), self.getSpecificEntropy())

                print ("    Finished iso-enthalpy lines")

            elif isovar in ('t','s'):
                print ("Can't use t or s as iso vars for Ts chart")

       
        print ("Assembling data")
        chartData[0] = isoz
        for n in range(len(isoT)):            
            chartData.append((isos[n],isoT[n]))
        print ("    Data assembly complete")

        return 


    def phChart(self, chartData):
        """
        this function will generate the data for a 
        P-h diagram for the specified fluid
        and store the data in the chartData list
        """

        print ("Creating saturation dome")

        #
        # get the bounds values for the iso lines
        Pcrit = self.getPCrit()
        Tcrit = self.getTCrit()
        T0 = self.getStanfordT0()
        self.setTemperature(T0)
        self.setQuality(0.0)
        self.setJobID(4)
        self.calcProps()
        P0 = self.getPressure()
        alpha = 0.99



        Pmax = self.getchartPmax()
        if Pmax == 0:
            Pmax = alpha*Pcrit

        Pmin = self.getchartPmin()

        if Pmin < P0:
            Pmin = P0

        Tmax = self.getchartTmax()
        if Tmax == 0:
            Tmax = Tcrit

        Tmin = self.getchartTmin()

        if Tmin < T0:
            Tmin = T0

        numlines = self.Thermo['chartNumLines']
        numSatPoints = 80
        
        #
        # set the number of sat dome points 
        # and deltaT for saturation calcs

        satDomeh = []
        satDomeP = []
        satDomex = []

        deltaT = (alpha*Tcrit - T0)/(numSatPoints + 1)
        Tcurr = T0
        #
        # start by going up the sat liq side of the dome
        self.setQuality(0.0)
        self.setJobID(4)
        for i in range(numSatPoints):
            Tcurr = Tcurr + deltaT
            self.setTemperature(Tcurr)
            self.calcProps()
            satDomeP.append(self.getPressure())
            satDomeh.append(self.getSpecificEnthalpy())
            #print (self.getRegion(), Tcurr, self.getQuality(), self.getSpecificEntropy())
#            satDomex.append(0.0)

        #
        # now add the sat vap side of the dome
        self.setQuality(1.0)
        Tcurr = alpha*Tcrit
        for i in range(numSatPoints):
            Tcurr = Tcurr - deltaT
            self.setTemperature(Tcurr)
            self.calcProps()
            satDomeP.append(self.getPressure())
            satDomeh.append(self.getSpecificEnthalpy())
            #print (self.getRegion(), Tcurr, self.getQuality(), self.getSpecificEntropy())
 #           satDomex.append(1.0)

        # make room for the labels
        chartData.append([])
        chartData.append((satDomeh, satDomeP))

        print ("    Finished saturation dome")

        #
        # get the chart iso line variable
        # and set the number of points
        
        chartIso = self.Thermo['chartIso']
        numPoints = 30
        deltaT = (Tcrit - T0)/(numPoints + 1)

        isoP = []
        isoh = []
        isoz = []
        
        for isovar in chartIso:

            if isovar == 'x':
                print ("Creating iso-quality lines")
                #
                # iso lines for quality
                dx = 1.0/(numlines + 1.0)
                x = 0.0
                #
                # set delta T based on T0 and Tcrit
                deltaP = ((Pcrit-(0.1*Pcrit)) - P0)/(numPoints + 1)

                numSets = len(isoP)
                for n in range(numlines):
                    setNum = n + numSets
                    isoP.append([])
                    isoh.append([])
                    Pcurr = P0
                    x = x + dx
                    self.setQuality(x)
                    self.setJobID(8)
                    isoz.append('x = %.2f' % (x))
                    for m in range(numPoints):
                        Pcurr = Pcurr + deltaP
                        self.setPressure(Pcurr)
                        self.calcProps()
                        isoP[setNum].append(self.getPressure())
                        isoh[setNum].append(self.getSpecificEnthalpy())
                        
                        #print (self.getRegion(), Tcurr, self.getQuality(), self.getSpecificEntropy())
                print ("    Finished iso-quality lines")
                         
            elif isovar == 't':
                #
                # iso lines for temperature
                # 

                # note that the user can specify 
                # a max > Tcrit, so we need to 
                # have a check on that for the iso lines
                print ("Creating isothermal lines")

                dt = (Tmax - Tmin)/(numlines + 1)

                numSets = len(isoP)
                tsat = Tmax

                for n in range(numlines):
                    tsat = tsat - dt
                    setNum = n + numSets
                    #
                    # create a new set of data
                    isoP.append([])
                    isoh.append([])
                    #
                    # add the saturation points for the 
                    # isothermal lines

                    self.setTemperature(tsat)
                    self.setQuality(0.0)
                    self.setJobID(4)
                    self.calcProps()
                    psat = self.getPressure()
                    isoP[setNum].append(psat)
                    isoh[setNum].append(self.getSpecificEnthalpy())
                    self.setQuality(1.0)
                    self.setJobID(4)
                    self.calcProps()
                    psat = self.getPressure()
                    isoP[setNum].append(psat)
                    isoh[setNum].append(self.getSpecificEnthalpy())
                    self.setJobID(1)
                    Pcurr = psat

                    deltaP = (Pcurr - Pmin)/(numPoints + 1)

                    isoz.append('T = %1.2e' % (tsat))
                    print ('---Saturation points set')
                    for m in range(numPoints+1):
                        Pcurr = Pcurr - deltaP
                        self.setPressure(Pcurr)
                        self.calcProps()
                        isoP[setNum].append(Pcurr)
                        isoh[setNum].append(self.getSpecificEnthalpy())
                        #print (self.getRegion(), Pcurr, tsat, self.getQuality(), self.getSpecificEnthalpy())
                print ("    Finished isothermal lines")
                
           
            elif isovar == 'v':
                #
                # iso lines for specific volume
                # 

                # note that the user can specify 
                # a max > Tcrit, so we need to 
                # have a check on that for the iso lines

                #
                # get pcrit and Psat(T0)
                print ("Creating iso-volume lines")
                vcrit = self.getvCrit()

                #
                # set vmax = vg(Tmin,x = 1)
                self.setPressure(Pmin)
                self.setQuality(0.1)
                self.setJobID(8)
                self.calcProps()
                vmax = self.getSpecificVolume()

                #
                # set vmin = v(Tmin, x = 0)
                self.setPressure(Pmax)
                self.setTemperature(Tcrit)
                self.setJobID(1)
                self.calcProps()
                vmin = self.getSpecificVolume()   

                #
                # now check the ratio of vmax/vmin
                # if it's large use ln to move thru
                # values instead of linear

                if vmax/vmin > 10:
                    dv = math.log(vmax/vmin)/(numlines + 1.0)
                    useLn = True
                else:                    
                    dv = (vmax - vmin)/(numlines + 1.0)
                    useLn = False

                numSets = len(isoP)
                vcurr = vmin
                
                for n in range(numlines):
                    if useLn:
                        vcurr = math.exp(math.log(vcurr) + dv)
                    else:
                        vcurr = vcurr + dv
                    #
                    # Now choose Pstart = min(Pmax, P(Tcrit,vcurr)
                    self.setTemperature(Tcrit)
                    self.setSpecificVolume(vcurr)
                    self.setJobID(2)
                    self.calcProps()
                    Pstart = min(Pmax,self.getPressure())
                    deltaP = (Pstart - Pmin)/(numPoints + 1)
                    self.setJobID(5)
                    
                    setNum = n + numSets
                    #
                    # create a new set of data
                    isoP.append([])
                    isoh.append([])
                    #
                    # add the saturation points for the 
                    # iso volume lines

                    Pcurr = Pstart
                    isoz.append('v = %1.2e' % (vcurr))
                    for m in range(numPoints):
                        Pcurr = Pcurr - deltaP
                        self.setPressure(Pcurr)
                        if Pcurr < Pcrit:
                            #
                            # Find v = v(Pcurr,x = 0.1)
                            self.setQuality(0.1)
                            self.setJobID(8)
                            self.calcProps()
                            vstop = self.getSpecificVolume()
                        else:
                            vstop = 0
                        if vstop < vcurr:
                            self.setSpecificVolume(vcurr)
                            self.setJobID(5)
                            self.calcProps()
                            isoP[setNum].append(self.getPressure())
                            isoh[setNum].append(self.getSpecificEnthalpy())
                            #print (self.getRegion(), Pcurr, vcurr, vstop, self.getQuality(), self.getSpecificEnthalpy())
                        else:
                            break

                print ("    Finished iso-volume lines")
                
            
            elif isovar == 's':
                #
                # iso lines for specific entropy
                # 

                # note that the user can specify 
                # a max > Tcrit, so we need to 
                # have a check on that for the iso lines

                #
                print ("Creating isentropic lines")

                #
                # set smax = sg(Pmin,x = 1)
                self.setPressure(Pmin)
                self.setQuality(1.0)
                self.setJobID(8)
                self.calcProps()
                smax = self.getSpecificEntropy()

                #
                # set smin = sf(Pmax, x = 0)
                self.setPressure(Pmin)
                self.setQuality(0.1)
                self.setJobID(8)
                self.calcProps()
                smin = self.getSpecificEntropy()

                #
                # now check the ratio of hmax/hmin
                # if it's large use ln to move thru
                # values instead of linear

                if smax/smin > 10:
                    ds = math.log(smax/smin)/(numlines + 1.0)
                    useLn = True
                else:                    
                    ds = (smax - smin)/(numlines + 1.0)
                    useLn = False

                numSets = len(isoP)
                scurr = smin
                
                for n in range(numlines):
                    if useLn:
                        scurr = math.exp(math.log(scurr) + ds)
                    else:
                        scurr = scurr + ds
                    self.setJobID(7)
                    self.setSpecificEntropy(scurr)
                    setNum = n + numSets
                    #
                    # create a new set of data
                    isoP.append([])
                    isoh.append([])
                    #
                    # add the saturation points for the 
                    # isentropic lines

                    Pcurr = Pmin

                    deltaP = (Pmax - Pmin)/(numPoints + 1)

                    isoz.append('s = %1.2e'% (scurr))
                    for m in range(numPoints):
                        Pcurr = Pcurr + deltaP
                        self.setPressure(Pcurr)
                        self.calcProps()
                        isoP[setNum].append(Pcurr)
                        isoh[setNum].append(self.getSpecificEnthalpy())
                        #
                        # put a few checks on the output to stop collecting data
                        # if x < 0.1 OR T > Tcrit, stop collecting points
                        if (self.getQuality() < 0.1) or (self.getTemperature() > 2.0*Tcrit):
                            break
                        #print (self.getRegion(), Pcurr, scurr, self.getQuality(), self.getSpecificEnthalpy())

                print ("    Finished isentropic lines")
                
            

            elif isovar in ('p','h'):
                print ("Can't use t or s as iso vars for Ph chart")

       
        print ("Assembling data")
        chartData[0] = isoz
        for n in range(len(isoP)):            
            chartData.append((isoh[n],isoP[n]))
        print ("    Data assembly complete")

        return 


        
def main():
    # parse command line options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h:cf:j:p:t:e:u:x:v:s:", ["help","chart=","iso=","numlines=","tmin=","tmax=","pmin=","pmax="])
        
    except getopt.error as msg:
        print (msg)
        print ("for help use --help")
        sys.exit(2)
    # process options
    props = ThermoProps()
    for o, arg in opts:
        if o == "--help":
            print ("to be added")
            sys.exit(0)
        if o == "-c":
            pass
        if o == "-p":
            props.setPressure(float(arg))
        if o == "-t":
            props.setTemperature(float(arg))
        if o == "-e":
            props.setEOS(arg)
        if o == "-f":
            props.setFluid(arg,'False')
        if o == "-v":
            props.setSpecificVolume(float(arg))
        if o == "-s":
            props.setSpecificEntropy(float(arg))
        if o == "-h":
            props.setSpecificEnthalpy(float(arg))
        if o == "-x":
            props.setQuality(float(arg))
        if o == "-j":
            props.setJobID(int(arg,10))
        if o == "-u":
            props.setUnitSys(int(arg,10))
        if o == "--chart":
            props.setchartType(arg)
        if o == "--iso":
            props.setchartIsoLines(arg)
        if o == "--numlines":
            props.setchartNumIsoLines(int(arg,10))
        if o == "--tmin":
            props.setchartTmin(float(arg))
        if o == "--tmax":
            props.setchartTmax(float(arg))
        if o =="--pmin":
            props.setchartPmin(float(arg))
        if o == "--pmax":
            props.setchartPmax(float(arg))




    if props.getchartType() == 'none':
        props.calcProps()
    else:
        props.makeChart()

    print (props.Thermo['FluidName'],props.Thermo['eos'], props.Thermo['MW'], props.Thermo['job'])
    print ('P = ', props.Thermo['Pressure'], 'T = ', props.Thermo['Temperature'], 'v = ', props.Thermo['SpVol'])
    print ('u = ',props.Thermo['SpEnergy'], 'h = ',props.Thermo['SpEnthalpy'],'s = ', props.Thermo['SpEntropy'])
    print ('Region = ',props.Thermo['Region'], 'x = ', props.Thermo['Quality'])
    # process arguments
    for arg in args:
        process(arg) # process() is defined elsewhere

if __name__ == "__main__":
    main()






