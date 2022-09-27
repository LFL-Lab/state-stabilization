 ## Python driver for Keithley 2400
##  Developed by USC Levenson-Falk Lab, driver convention established by James Farmer
##  Contributors:
##                  William S. Sager : williamsager7@gmail.com

'''
This is a driver for the Keithley 2400 SourceMeter.
'''

# Libraries
import pyvisa as visa
from visa import ResourceManager


#class with all necessary
class K2400SM():

    def __init__(self, address, reset=True):
        '''Initialization sequence for the SourceMeter'''
        # Open the K2400SM
        self.inst = ResourceManager().open_resource(address)
        #Shorthand
        self.w = self.inst.write
        self.r = self.inst.read
        self.q = self.inst.query
        #Print K2400SM string for confirmation
        #print(self.q('*IDN?'))
        if reset:
            self.w('*RST')

    ## The following are more basic source-measure commands
    ## They can be found in table 3-5 of the Keithley 2400 SM manual
    def id(self):
        '''Queries instrument identification info'''
        return self.q('*IDN?')

    def toggleDisplay(self, choice ):
        ''' Enable or disable front panel display.
        Arguments: 'on', 'off'.
        '''
        displayOff = ':DISPlay:ENABle OFF'
        displayOn = ':DISPlay:ENABle ON'
        if choice == 'on':
            #toggle display on
            self.w(displayOn)
        elif choice == 'off':
            #toggle display off
            self.w(displayOff)

    def sourceFunction(self, choice):
        ''' Select source function(voltage or current). Arguments: 'v', 'i'. '''
        currentMode = 'SOUR:FUNC CURRent'
        voltageMode = 'SOUR:FUNC VOLTage'
        if choice == 'i':
            self.w(currentMode)
        if choice == 'v':
            self.w(voltageMode)

    def fixedSourceMode(self,choice):
        '''
        Selects fixed sourcing mode for voltage or current.
        Choose between voltage and current.
        Arguments: 'v','i'
        '''
        currentMode = ':SOURce:CURRent:MODE FIXed'
        voltageMode = ':SOURce:VOLTage:MODE FIXed'
        if choice =='i':
            self.w(currentMode)
        elif choice == 'v':
            self.w(voltageMode)

    def setSourceRange(self, choice, range):
        '''
        Set the range for the current or voltage source.
        Arguments:'i', 'v'; n
        '''
        currentRange = ':SOURce:CURRent:RANGe '+str(range)
        voltageRange = ':SOURce:VOLTage:RANGe '+str(range)
        if choice == 'i':
            self.w(currentRange)
        elif choice == 'v':
            self.w(voltageRange)

    def setLevel(self, choice, level):
        '''
        Set the amplitude for the current or voltage source
        Arguments: 'i', 'v'; n
        '''
        currentLevel = ':SOURce:CURRent:LEVel '+str(level)
        voltageLevel = ':SOURce:VOLTage:LEVel '+str(level)
        if choice == 'i':
            self.w(currentLevel)
        elif choice == 'v':
            self.w(voltageLevel)

    def selectMeasureFN(self, choice):
    ### Throwing a data type error on instrument, error code 104 -JF(?) ###
    ### I think this error was likely due to the arguments being CURR 
    ### And VOLT instead of CURRent and VOLTage - WSS ###   
    
    ### Problem still exists. SCPI strings are case independent, but are documented like that
    ### for shorthand typing. You can abbreviate to only the uppercase letters.
    ### For Example, SENSe = SENS = sense = sens. These are all valid and equivalent.
    ### Also, you often see something like "[SENSe:]FUNCtion <argument>" in a manual.
    ### In this case, you can omit anything inside brackets and it will still work. - JF ###

        '''
        Select the measure function
        Arguments: 'i', 'v'
        '''
        if choice == 'i':
            self.w("SENS:FUNC 'CURR'")
        elif choice == 'v':
            self.w("SENS:FUNC 'VOLT'")

    def setComplianceLevel(self, choice, level):
    ### This is the one that is indiscriminate about current and voltage.
    ### It sets the compliance current even if choice = 'v'.
    ### Likely the string in setVoltage definition, changed as below - JF ###
        '''
        This function sets the current or voltage compliance.
        Arguments: 'i','v'; n
        '''
        setCurrent = 'SENSe:CURRent:PROTection '+str(level)
        setVoltage = 'SENSe:VOLTage:PROTection '+str(level) # Changed from SENS:CURR:PROT to SENS:VOLT:PROT -JF
        if choice == 'i': #Set current Level
            self.w(setCurrent)
        elif choice == 'v': #Set voltage level
            self.w(setVoltage)

    def setMeasurementRange(self, choice, level):
    ### indiscriminant on 'i' vs 'v'. In both cases, it sets the compliance current. -JF(?) ###
    ### I think this was because the arguments were 'SENSe:...' instead of 'SOURce...'- WSS ###
    ### Sorry, this was meant for the function above (setComplianceLevel), changing back to SENSe - JF ###

        '''
        This function sets the measurement range for current or voltage.
        Arguments: 'i','v'; n
        '''
        currentRange = 'SENSe:CURRent:RANGe '+str(level)
        voltageRange = 'SENSe:VOLTage:RANGe '+str(level)
        if choice == 'i': #set current range
            self.w(currentRange)
        elif choice == 'v':#set voltage range
            self.w(voltageRange)

    def outputState(self, choice):
        '''
        This function sets the output state.
        Arguments: 'on', 'off'
        '''
        if choice == 'on':#Set output state to on
            self.w(':OUTPut ON')
        elif choice == 'off':#Set output state to off
            self.w(':OUTPut OFF')
    
    def read(self):         # sets the instrument to read, but does not acquire data to computer -JF
        '''This function triggers and acquires reading.'''
        self.q(':READ?')

    ## End of functions marked in table 3-5
    ## Start of functions marked in Table 4-2, Remote Ohms Commands

    ### Considering consolodating these along with functions from
    ### Table 3-5 etc. Will wait until all functions have been added
    ### To the driver before deciding on organization though. -WSS

    def selectOhmsFN(self):
        ''' Selects ohms function'''
        self.w(":SENSe:FUNCtion RESistance")

    def setOhmsRange(self, level):
        '''
        Selects the ohms range.
        Argument: n
        '''
        ##integrating this function into selectMeasurementRange function -WSS
        ohmsRange =':SENSe:RESistance:RANGe '+ str(level)
        self.w(ohmsRange)

    def selectOhmsMode(self, mode):
        '''
        Selects the ohms mode
        Arguments: 'm' <-- manual, 'a' <-- auto
        '''
        if mode == 'a':
            self.w(':SENSe:RESistance:MODE MANual')
        elif mode == 'm':
            self.w(':SENSe:RESistance:MODE AUTO')

    def offsetCompensation(self, choice):
        '''
        This function enables/disables the offset compensation
        Arguments: 'on', 'off'
        '''
        if choice == 'on':
            self.w(':SENSe:RESistance:OCOMpensated ON')
        elif choice == 'off':
            self.w(':SENSe:RESistance:OCOMpensated OFF')

    def setOhmsCompliance(self, mode, level):
        '''
        This sets the current or voltage compliance
        Arguments: 'i', 'v'; n
        '''
        voltageMode = ':SENSe:CURRent:PROTection ' + str(level)
        currentMode = ':SENSe:VOLTage:PROTection ' + str(level)
        if mode == 'i':
            self.w(currentMode)
        elif mode == 'v':
            self.w(voltageMode)

    def SensingState(self, mode):
        '''
        This function selects between 2 and 4 wire
        sensing states.
        Arguments: '2', '4'
        '''
        if mode == '2':
            self.w(":SYSTem:RSENse OFF")
        elif mode == '4':
            self.w(":SYSTem:RSENse ON")

    ## End of remote ohms commands from table 4-3
    ## (Did not duplicate functions also found in table 3-5)

    ## Start of remote range and digits programming commands
    ## Found in Table 7-1

    def selectMeasurementRange(self, mode, choice, level, enable):
        '''
        This function selects the current or voltage measurement range
        If auto is selected the range will be enabled/disabled
        Arguments: 'i', 'v', 'R'; 'a' <-- auto, 'm' <-- manual; n; 'on', 'off'
        '''
        manCurrentMode = ':SENSe:CURRent:RANGe ' + str(level)
        manVoltageMode = ':SENSe:VOLTage:RANGe ' + str(level)
        autCurrentMode = ':SENSe:CURRent:RANGe:AUTO ' + str(enable)
        autVoltageMode = ':SENSe:VOLTage:RANGe:AUTO ' + str(enable)
        autResistanceMode = ':SENSe:RESistance:RANGe:AUTO ' + str(enable)
        manResistanceMode = ':SENSe:RESistance:RANGe ' + str(level)
        if choice == 'a':
            if mode == 'i':
                self.w(autCurrentMode)
            elif mode == 'v':
                self.w(autVoltageMode)
            elif mode == 'R':
                self.w(autResistanceMode)

        elif choice == 'm':
            if mode == 'i':
                self.w(manCurrentMode)
            elif mode == 'v':
                self.w(manVoltageMode)
            elif mode == 'R':
                self.w(manResistanceMode)

    def displayDigits(self, level):
        '''
        This function sets the display digits
        Arguments: n (n = 4, 5, 6, 7)
        '''
        nDigits = ':DISPlay:DIGits ' + str(level)
        self.w(nDigits)

    ## End functions found in table 7-1

    ## Start functions found in table 7-3, Speed commands
    def setSpeed(self, mode, level):
        '''
        This function sets the speed for current, voltage
        or resistance.
        Arguments: 'i', 'v', 'R'; n
        '''
        currentMode = ':SENSe:CURRent:NPLCycles ' + str(level)
        voltageMode = ':SENSe:VOLTage:NPLCycles ' + str(level)
        resistanceMode = ':SENSe:RESistance:NPLCycles ' + str(level)
        if mode == 'i':
            self.w(currentMode)
        elif mode == 'v':
            self.w(voltageMode)
        elif mode == 'R':
            self.w(resistanceMode)
    ## End functions found in table 7-3

    ## Start functions found in table 7-4
    def filterControl(self, fType, level, enable):
        '''
        This function selects the filter type (repeat or moving)
        sets the filter count and enables/disables the filter state
        Arguments: 'repeat', 'moving'; n; 'on', 'off'
        '''
        filterType = ':SENSe:AVERage:TCONtrol ' + str(fType)
        self.w(filterType)
        filterCount = ':SENSe:AVERage:COUNt ' + str(level)
        self.w(filterCount)
        if enable == 'on':
            self.w(":SENSe:AVERage ON")
        elif enable == 'off':
            self.w(":SENSe:AVERage OFF")
    ## End functions found in table 7-4

    ## Start functions found in table 8-1: Rel commands
    def setNullValue(set, level):
        '''
        This function defines a null value.
        Arguments: n
        '''
        self.w(':CALCulate2:NULL:OFFSet ' + str(level))
    def enRelState(self, enable):
        '''
        This function enables/disables the rel state.
        Arguments: 'on', 'off'
        '''
        if choice == 'on':
            self.w(':CALCulate:NULL:STATe ON')
        elif choice == 'off':
            self.w(':CALCulate:NULL:STATe OFF')

    def nullAcquire(self):
        '''This function automatically acquires the rel value'''
        self.w(':CALCulate2:NULL:ACQUire')

    ## End function found in table 8-1

    ## Start functions found in table 8-3: Math commands
    def defineMathExpressionName(self, name):
        '''
        This function sets the expression name. There are
        the following pre-defined names: "POWER", "OFFCOMPOHM"
        "VOLTCOEF", "VARALPHA"
        '''
        expressionName = ':CALCulate:MATH:NAME ' + str(name)
        self.w(expressionName)

    def mathDataQuery(self):
        '''This function queries math data'''
        self.w("CALCulate:DATA?")

    ## End functions found in table 8-3

    ## Start functions found in table 8-4: User defined math functions
    def UserDefinedFNUnits(self, units):
        '''
        This function sets the units for a user defined math function.
        Arguments: name <-- must be 3 ASCII characters
        '''
        unitString = ':CALCulate:MATH:UNITs ' + str(name)
        self.w(unitString)

    def defineMathName(self, name):
        '''
        This function defines a math name.
        Arguments: name (name = "user-name")
        '''
        nameString = ':CALCulate:MATH:NAME ' + str(name)
        self.w(nameString)

    def defineMathFormula(self, formula):
        '''
        This function defines a valid math operator.
        Arguments: send a formula as a string to this function for it to be
        evaluated. Terms are as follows:
                NAME 'VOLT', 'CURR', 'RES', 'TIME';
                OPERATOR '+', '-', '*', '/', '^', 'log', 'sin', 'cos', 'tan', 'exp'
        Example: formula = "(((RES - 10e3) / 10e3) * 100)"
        '''
        formulaString = ':CALCulate:MATH:EXPR ' + str(formula)
        self.w(formulaString)

    def enMathState(self, enable):
        '''
        This function enables/disables the math state
        Arguments: 'on', 'off'
        '''
        if enable == 'on':
            self.w(':CALCulate:STATe ON')
        elif enable == 'off':
            self.w(':CALCulate:STATe OFF')

    ## End functions from table 8-4

    ## Start functions from table 9-1: Data store commands

    def readBufferContents(self):
        '''This function reads the content of buffer'''
        self.w(':TRACe"DATA?')

    def clearBuffer(self):
        '''This function clears the buffer.'''
        self.w(':TRACe:CLEar')

    def bufferMemoryStatus(self):
        '''This function reads the memory status'''
        self.w(':TRACe:FREE?')

    def setBufferSize(self, size):
        '''
        This function sets the buffer size.
        Arguments: n
        '''
        sizeString = ':TRACe:POINts '+ str(n)
        self.w(sizeString)

    def nStoredReadings(self):
        '''This function queries the number of stored readings'''
        self.w(':TRACe:POINts:ACTual?')

    def specifyReadingSource(self, name):
        '''
        This function specifies the reading source.
        Arguments: 'sense','calc1','calc2'
        '''
        if name == 'sense':
            self.w('TRACe:FEED SENSe')
        elif name == 'calc1':
            self.w(':TRACe:FEED CALCulate1')
        elif name == 'calc2':
            self.w(':TRACe:FEED CALCulate2')

    def bufferStartStop(self, name):
        '''
        This function starts or stops the buffer. The argument 'next'
        fills the buffer and stops, the argument 'never' disables the buffer
        Arguments: 'next','never'
        '''
        if name == 'next':
            self.w(':TRACe:FEED:CONTrol NEXT')
        elif name == 'never':
            self.w(':TRACe:FEED:CONTrol NEVer')

    def timeStampFormat(self, name):
        '''
        This function selects the timestamp format. The argument 'first'
        refers to the first buffer reading, the argument 'delta' refers to
        the time between buffer readings.
        Arguments: 'first', 'delta'
        '''
        if name == 'first':
            self.w(':TRACe:TSTamp:FORMat ABSolute')
        elif name == 'delta':
            self.w(':TRACe:TSTamp:FORMat DELTa')

    def bufferStatisticName(self, name):
        '''
        This function selects the buffer statistic name.
        Arguments: 'mean', 'sdev' <-- standard deviation, 'max', 'min',
                    'pk2pk' <-- peak to peak
        '''
        if name == 'mean':
            self.w(':CALCulate3:FORMat MEAN')
        elif name == 'sdev':
            self.w(':CALCulate3:FORMat SDEViation')
        elif name == 'max':
            self.w(':CALCulate3:FORMat MAXimum')
        elif name == 'min':
            self.w(':CALCulate3:FORMat MINimum')
        elif name == 'pk2pk':
            self.w(':CALCulate3:FORMat PKPK')

    def readBufferStatData(self):
        '''
        This function reads the buffer statistic data.
        Note: if TRACe:FEED is set to :SENSe[1] this command will return
        one V, I, Resistance and Math result.
        '''
        self.w(':CALCulate3:DATA?')

    ## End functions found in table 9-1

    ## Start functions found in Table 10-2: Source memory saved configuration

    ###Less sure of what these functions do, will find better definitions
    ###And add better descriptions.-WSS
    def iIntegrationRate(self):
        '''This is a function for the current integration rate'''
        self.w(':SENSe:CURRent:NPLCycles')

    def RIntegrationRate(self):
        '''This is a function for the resistance integration rate'''
        self.w(':SENSe:RESistance:NPLCycles')

    def vIntegrationRate(self):
        '''This is a functino for the voltage integration rate'''
        self.w(':SENSe:VOLTage:NPLCycles')

    def concurrentFN(self):
        '''This is a function that sets concurrent functions'''
        self.w(':SENSe:FUNCtion:CONCurrent')

    def enableFN(self, en):
        '''
        This is a function that enables or disables functions
        Arguments: 'enable', 'disable'
        '''
        if en == 'enable':
            self.w(':SENSe:FUNCtion:ON')
        elif en == 'disable':
            self.w(':SENSe:FUNCtion:OFF')

    ### Section 10-1 10-2 not finished

    ## Section 10-3:Linear and log staircase sweep commands

    def currentSweep(self, iStart, iStop, iStep, iCenter, iSpan):
        '''
        This function handles current sweep modes.
        Arguments: 'iStart'(n); 'iStop'(n); 'iStep'(n); 'iCenter'(n); 'iSpan'(n);
        '''
        startString = ':SOURce:CURRent:STARt ' + str(iStart)
        stopString = ':SOURce:CURRent:STOP ' + str(iStop)
        stepString = ':SOURce:CURRent:STEP ' + str(iStep)
        centerString = ':SOURce:CURRent:CENTer ' + str(iCenter)
        spanString = ':SOURce:CURRent:SPAN ' + str(iSpan)
        self.w(':SOURce:CURRent:MODE SWEep')
        self.w(startString)
        self.w(stopString)
        self.w(stepString)
        self.w(centerString)
        self.w(spanString)

    def voltageSweep(self, vStart, vStop, vStep, vCenter, vSpan):
        '''
        This function handles voltage sweep modes.
        Arguments: 'vStart'(n); 'vStop'(n); 'vStep'(n); 'vCenter'(n); 'vSpan'(n);
        '''
        startString = ':SOURce:VOLTage:STARt ' + str(vStart)
        stopString = ':SOURce:VOLTage:STOP ' + str(vStop)
        stepString = ':SOURce:VOLTage:STEP ' + str(vStep)
        centerString = ':SOURce:VOLTage:CENTer ' + str(vCenter)
        spanString = ':SOURce:VOLTage:SPAN ' + str(vSpan)
        self.w(':SOURce:VOLTage:MODE SWEep')
        self.w(startString)
        self.w(stopString)
        self.w(stepString)
        self.w(centerString)
        self.w(spanString)

    def sweepRange(self, ranging):
        '''
        This function selects source  ranging.
        Arguments: 'best', 'auto', 'fixed'
        '''
        if ranging == 'best':
            self.w(':SOURce:SWEep:RANGing BEST')
        elif ranging == 'auto':
            self.w(':SOURce:SWEep:RANGing AUTO')
        elif ranging == 'fixed':
            self.w(':SOURce:SWEep:RANGing FIXed')

    def sweepScale(self, scale):
        '''
        This function selects sweep scale.
        Arguments: 'linear','logarithmic'
        '''
        if scale == 'linear':
            self.w(':SOURce:SWEep:SPACing LINear')
        elif scale == 'logarithmic':
            self.w(':SOURce:SWEep:SPACing LOGarithmic')

    def nSweepPoints(self, nPoints):
        '''
        This function sets the number of sweep points.
        Arguments: n
        '''
        pointsString = ':SOURce:SWEep:POINts ' + str(nPoints)
        self.w(pointsString)

    def sweepDirection(self, direction):
        '''
        This function sets the sweep direction (stop to start).
        Arguments: 'up','down'
        '''
        if direction == 'up':
            self.w(':SOURce:SWEep:DIREction UP')
        elif direction == 'down':
            self.w(':SOURce:SWEep:DIREction DOWn')

    ## End functions found in table 10-3

    ## Start functions found in table 10-7
    def memorySweep(self, mPoints, mLocation):
        '''
        This function selects memory sweep mode, sets the number of sweep
        points and selects memory start location.
        Arguments: n (number of points); n (source memory start location);
        '''
        self.w(':SOURce:FUNCtion MEM')
        pointsString = ':SOURce:MEMory:POINts ' + str(mPoints)
        self.w(pointsString)
        locationString = ':SOURce:MEMory:STARt ' + str(mLocation)

    def memorySaveRecall(self, SorR, location):
        '''
        This function either saves, or recalls a setup.
        Arguments: 'save', 'recall'; n (n = memory location)
        '''
        if SorR == 'save':
            sString = ':SOURce:MEMory:RECall ' + str(location)
            self.w(sString)
        elif SorR == 'recall':
            rString = ':SOURce:SAVE ' + str(location)
            self.w(rString)

    ## End functions found in table  10-7

    ## Start functions found in table 11-1: Remote trigger commands
    def initSM(self):
        '''This function takes the source meter out of idle state.'''
        self.w(':INITiate')

    def abortSM(self):
        '''
        This function aborts the operation, and returns the source meter
        to idle state.
        '''
        self.w(':ABORt')

    def setArmCount(self, armCount):
        '''
        This function sets the arm count.
        Arguments: n (arm count)
        '''
        armString = ':ARM:COUNt ' + str(armCount)
        self.w(armString)

    def armControlSource(self, name):
        '''
        This function specifies the arm control source.
        Arguments: 'immediate','tlink', 'timer','manual', 'bus','nstest',
                    'pstest', 'bstest'
        '''
        controlString = ':ARM:SOURce ' + name
        self.w(controlString)

    def armTimer(self, interval):
        '''
        This function sets the arm layer timer.
        Arguments: n (interval)
        '''
        intervalString = ':ARM:TIMer ' + str(interval)
        self.w(intervalString)

    def armBypassControl(self, name):
        '''
        This function controls the arm bypass.
        Arguments: 'source', 'acceptor'
        '''
        nameString = ':ARM:DIRection ' + name
        self.w(nameString)

    def armInLine(self, lineNumber):
        '''
        This function selects the arm layer input line.
        Arguments: n (input line #)
        '''
        lineString = ':ARM:ILINe ' + str(lineNumber)
        self.w(lineString)

    def armOutLine(self, lineNumber):
        '''
        This function selects the arm layer output line.
        Arguments: n (output line #)
        '''
        lineString = ':ARM:OLINe ' + str(lineNumber)

    def armOutEvents(self, event):
        '''
        This function selects the arm layer output events.
        Arguments: 'tenter', 'texit', 'none'
        '''
        tentString = ':ARM:OUTPut TENTer'
        texString = ':ARM:OUTPut TEXit'
        noneString = ':ARM:OUTPut NONE'
        if event == 'tenter':
            self.w('tentString')
        elif event == 'texit':
            self.w(texString)
        elif event == 'none':
            self.w(noneString)

    def clearTrig(self):
        '''This function clears any pending triggers immediately'''
        self.w(':TRIGger:CLEar')

    def setTrigCount(self, n):
        '''
        This function sets the trigger count
        Arguments: n
        '''
        trigString = ':TRIGger:COUNt ' + str(n)
        self.w(trigString)

    def setTrigDelay(self, n):
        '''
        This function sets the trigger delay.
        Arguments: n
        '''
        trigString = ':TRIGger:DELay ' + str(n)

    def specTrigCntrlSrc(self, source):
        '''
        This function specifies the triggers control source.
        Arguments: 'immediate', 'tlink'
        '''
        imString = ':TRIGer:SOURce IMMediate'
        tString = ':TRIGger:SOURce TLINk'
        if source == 'immediate':
            self.w(imString)
        elif source == 'tlink':
            self.w(tString)

    def cntrlTrigBypass(self, name):
        '''
        This function controls the trigger bypass
        Arguments: 'source', 'acceptor'
        '''
        sString = ':TRIGger:DIRection SOURce'
        aString = ':TRIGger:DIRection ACCeptor'
        if name == 'source':
            self.w(sString)
        elif name == 'acceptor':
            self.w(aString)

    def trigIOLine(self, IO, line):
        '''
        This function selects the triggers input or output line.
        Arguments: 'in','out'; line
        '''
        if IO == 'in':
            self.w(':TRIGger:ILINe ' + str(line))
        elif IO == 'out':
            self.w(':TRIGger:OLINe ' + str(line))

    def trigIOEvents(self, IO, events):
        '''
        This function selects the input or output layer events.
        Arguments: 'in', 'out'; 'source','delay','sense','none'
        '''
        iSrcString = ':TRIGger:INPut SOURce'
        oSrcString = ':TRIGger:OUTPut SOURce'
        iDstring = ':TRIGger:INPut DELay'
        oDstring = ':TRIGger:OUTPut DELay'
        iSenString = ':TRIGger:INPut SENSe'
        oSenString = ':TRIGger:OUTPut SENSe'
        iNstring = ':TRIGger:INPut NONE'
        oNstring = ':TRIGger:OUTPut NONE'

        if IO == 'in':
            if events == 'source':
                self.w(iSrcString)
            elif events == 'delay':
                self.w(iDstring)
            elif events == 'sense':
                self.w(iSenString)
            elif events == 'none':
                self.w(iNstring)
        elif IO == 'out':
            if events == 'source':
                self.w(oSrcString)
            elif events == 'delay':
                self.w(oDstring)
            elif events == 'sense':
                self.w(oSenString)
            elif events == 'none':
                self.w(oNstring)

    def trigSM(self):
        ''' This function is defined as "Trigger Sourcemeter (if BUS source selected)"'''
        self.w(*TRG)
