'''
A class to control the Agilent N5222A.
Contributor: Haimeng Zhang
'''
# start a measurement and read data
from visa import ResourceManager
from time import sleep

class VNA():
    
    def __init__(self,address):
        '''
        create a new VNA instance
        Parameters:
        address: str
            A string containing the VNA address of the device
        '''
        self.inst=ResourceManager().open_resource(address)#a visa object
        self.inst.write_termination = '\n'
        self.inst.read_termination = '\n'
        #self.inst.query_termination = '\n'
        self.w=self.inst.write
        self.r=self.inst.read
        self.q=self.inst.query
        '''
        if reset:
            self.w('*RST')
        '''
    def ID(self):
        return self.q('*idn?')
    
    def reset(self):
        '''resets the AWG, restoring all parameters to defaults, etc.
           -JF'''
        self.w('*rst')
    def InitSetup(self):
        '''
        perform initial setup
        Delete previous parameters, put trigger on hold, turn display on.
        '''
        self.w('CALCulate:PARameter:DELete:ALL')#Deletes all measurements on the PNA.
        self.w('SENSe:SWEep:MODE HOLD')
        self.w('DISPlay:ENABLE ON')
    def disconnect(self):
        '''
        Turns output off and disconnect
        '''
        #self.outputOff()
        self.inst.close()
        
    def clearwindow(self,wnum):
        '''
        wnum:int
            number of window to clear
        '''
        traces=self.q('DISPlay:WINDow{}:CATalog?'.format(wnum))#Returns the trace numbers for the specified window.
        
    def setup(self,measName,npoints=None,startFreq=None,stopFreq=None,centFreq=None,spanFreq=None):
        '''
        Measurement setup
        Parameters:
        centFreq:
        spanFreq:
        startFreq: str
        stopFreq: str
        nPoints:
        '''
        #self.w('CALCulate:PARameter:DELete:ALL')
        '''
        for n in portNums: 
            self.w('DISPlay:WINDow:ENABle 1'.format(n))# disable or enable analyzer display information
            #self.clearwindow(n)
        '''
        
        self.w("CALCulate:PARameter:SELect \'{}\'".format(measName))#select the measurement
        self.w("SENSe1:SWEep:TYPE LIN")
        if startFreq and stopFreq: 
            self.w('SENSe1:FREQuency:STARt {}'.format(startFreq)) 
            self.w('SENSe1:FREQuency:STOP {}'.format(stopFreq))
        if centFreq and spanFreq: 
            self.w('SENSe1:FREQuency:CENTer {}'.format(centFreq))
            self.w('SENSe1:FREQuency:SPAN {}'.format(spanFreq))
        if npoints:
            self.w('SENSe1:SWEep:POINts '+str(npoints))#Sets the number of data points for the measurement
            
    def setupQubitSpectroscopy(self,npoints = 801,freq = 5.59129,power = -30):
        '''Prepares the VNA for a spectroscopic measurement of qubit frequency using 
           another signal generator - JF'''
        measName,Sparam = str(self.q("CALC:PAR:CAT?")).strip('\"').split(',')
        self.w("CALC:PAR:SEL \'{}\'".format(measName))#select the measurement
        self.w("TRIG:SOUR MAN") # sets up for a manual trigger
        self.w("SENS:SWE:MODE CONT")
        #self.w("CALCulate1:CORRection:EDELay 49.536744ns") # set the electrical delay
        self.w("SOUR1:POW1 {}".format(power)) # sets the RF power level
        #self.w("CALCulate1:PARameter:DELete \'{}\'".format(measName)) #delete existing measurement
        #self.w("DISP:WIND1:TRAC1:DEL") #delete current trace
        #self.w("CALCulate1:PARameter:DEFine:EXT \'{}\',S21".format(measName)) # define the measurement
        #self.w("DISP:WIND1:TRAC1:FEED \'{}\'".format(measName)) # Show the trace on VNA
        self.w("CALC:FORM SMIT") # sets the measurement format to SMITH
        self.w("SENS:SWE:TYPE CW") # change to CW time sweep
        self.w("SENS:FREQ {}GHz".format(freq)) # set the output frequency
        self.w('SENS:SWE:POIN {}'.format(npoints)) # set the number of points in the sweep
        self.w("TRIG:CHAN:AUX1:DUR 1E-6") # output trigger duration
        self.w("TRIG:CHAN:AUX1:INT POIN") # send trigger for each point in sweep
        self.w("TRIG:CHAN:AUX1 ON") # enables output trigger on AUX 1
        self.w("SENS:BWID 50") # sets the IF bandwidth to 10 Hz
    
    def qubitSpectroscopyMeas(self):
        '''Initiates the measurement and pulls data - JF'''
        #self.w("DISPlay:WINDow1:STATE ON") # display window
        #self.w("SENS:SWE:MODE SING") # initiate a single measurement sweep
        self.w("INIT:IMM") # Force trigger of sweep
        delay = self.q("SENS:SWE:TIME?") # Ask how long it will take
        sleep(float(delay))                        # halt program for the measurement duration
        self.q('*OPC?')                     # verify VNA is finished
        self.w("FORM ASCii,0")            # specify data format
        results = self.q("CALC:DATA? FDATA") # pull data
        return results
       
        
    def sMeas(self,Ports,testname,setparms=None):
        '''
        Perform and save an s-parameter measurement

        Parameters:
        sPorts:str
            number of ports separated by comma
        savedir: str
        testname: str
        '''
        sParms=[]
        sPorts=Ports.split(",")

        if len(sPorts)<1 or len(sPorts)>2:
            raise ValueError('Please Specify a number of ports between 1 and 2. '
                             'Currently, {} ports are specified.'.format(str(len(sPorts))))
        
        
        '''
        for i in sPorts:
            for j in sPorts:
                s='S{}{}'.format(i,j)
                sParms.append(s)
                measName='meas'+s
                print(measName)
                self.w("CALCulate:PARameter:DEFine:EXTended \'{}\',\'{}\'".format(measName,s))
                # Creates a measurement but does NOT display it.
                self.w("DISPlay:WINDow{}:TRACe{}:FEED \'{}\'".format(i,j,measName))
                #Creates a new trace <tnum> and associates (feeds) a measurement <name> to the specified window<wnum>.
        '''
        self.w("DISPlay:WINDow1:STATE ON")
        s='S11'
        measName='Ch1_'+s
        self.w("CALCulate:PARameter:DEFine:EXTended \'{}\',\'{}\'".format(measName,s))
        self.w("DISPlay:WINDow1:TRACe1:FEED \'{}\'".format(measName))

        if setparms:
            self.setup(measName,**setparms)
        else:
            self.setup(measName)
        self.inst.timeout=450000
        self.w("SENSe1:SWEep:MODE SINGle")#Sets the number of trigger signals the specified channel will ACCEPT.
        self.q('*OPC?')
        filename='{}.s1p'.format(testname)
        print(filename)
        #Save Data
        print("CALCulate:DATA:SNP:PORTs:SAVE \'1\',\'{}\'".format(filename))

        self.w("CALCulate:PARameter:SELect \'{}\'".format(measName))
        #select the measurement
        self.w("CALCulate:DATA:SNP:PORTs:SAVE \'1\',\'C:\Documents\{}\'".format(filename))
        #Saves SNP data from the selected measurement for the specified ports.
        self.q('*OPC?')
        self.inst.timeout = 2000
        self.savescreen(measName,filename)
        self.savecsv(measName,filename)
        self.outputOff()
        
    def savescreen(self,measName,direc,filename):
        #select measurement
        self.w("CALCulate:PARameter:SELect \'{}\'".format(measName))
        self.w("MMEMory:STORe:SSCReen \'{}\{}.bmp\'".format(direc,filename))
        
    def savecsv(self,measName,filename):
        self.w("CALCulate:PARameter:SELect \'{}\'".format(measName))
        self.w("MMEMory:STORe:DATA \"C:\Documents\\test.csv\",\"CSV Formatted Data\",\"Displayed\",\"RI\",-1")
        
    def plotcsv(self):
        self.q("MMEMory:TRANsfer? \'C:\Documents\\test.csv\'")
        
    def outputOff(self):
        '''
        Turn output off and put trigger on hold
        '''
        self.w('SENSe1:SWEep:MODE HOLD')
    def MeasSetupExam(self):
        '''
        This example sets up:
        S11 traces on two channels
        10 data points
        Sweep time of 2 seconds - this is slow enough to allow us to watch as each trace is triggered.
        '''
        self.w("SYST:FPReset")
        #Deletes all traces, measurements, and windows.
        #Resets the analyzer to factory defined default settings and creates a S11 measurement named "CH1_S11_1".
        self.w("DISPlay:WINDow1:STATE ON")#Create and turn on window/channel 1
        self.w("CALCulate1:PARameter:DEFine:EXT 'MyMeas1',S11")#Define a measurement name, parameter
        self.w("DISPlay:WINDow1:TRACe1:FEED 'MyMeas1'")#Associate ("FEED") the measurement name ('MyMeas') to WINDow (1)
        self.w("SENS1:SWE:TIME 2")#Set slow sweep so we can see
        self.w("SENS1:SWE:POIN 10")#set number of points to 10
        self.w("SENS1:SWE:MODE HOLD")

    def current_status(self):
        meas_1=self.q("CALCulate1:PARameter:CATalog?")
        meas_2=self.q("CALCulate2:PARameter:CATalog?")
        window=self.q("DISPlay:CATalog?")
        trace1=self.q("DISPlay:WINDow1:CATalog?")
        trace2=self.q("DISPlay:WINDow2:CATalog?")
        print("current measurement in channel 1: {}".format(meas_1))
        print("current measurement in channel 2: {}".format(meas_2))
        print("current window: {}".format(window))
        print("trace in window 1:".format(trace1))
        print("trace in window 2:".format(trace2))
    def readdata(self,measName):
        self.w("CALCulate1:PARameter:SELect '{}'".format(measName))
        self.w("FORMat ASCii,0")
        self.q("CALCulate1:DATA? RDATA")
    #def S21meas(self,power,start=4,stop=8):
    def S21meas(self):
        measName,Sparam = str(self.q("CALC:PAR:CAT?")).strip('\"').split(',')
        self.w("CALC:PAR:SEL \'{}\'".format(measName))#select the measurement
        #self.w("CALCulate1:PARameter:DELete 'TestMeas'")
        #self.w("DISPlay:WINDow2:STATE ON")
        #self.w("CALCulate1:PARameter:DEFine:EXT 'TestMeas',S21")
        #self.w("DISPlay:WINDow2:TRACe2:FEED 'TestMeas'")
        #self.w('SENSe1:FREQuency:STARt {}ghz'.format(start))
        #self.w('SENSe1:FREQuency:STOP {}ghz'.format(stop))
        #self.w("SOUR1:POW1 {}".format(power))

        
        
        
        self.w("CALCulate1:FORMat SMITh")
        
        #self.w("SENS1:SWE:TIME 185.326e-3")#Set slow sweep so we can see
        #self.w("SENS1:SWE:POIN 201")#set number of points to 10
        self.w("SENS1:SWE:MODE CONT")
        #self.w("CALCulate1:CORRection:EDELay 51.733003NS")
        sleep(2)
        self.w("FORMat ASCii,0")
        results=self.q("CALCulate1:DATA? FDATA")
        return str(results)
    def s21memory(self):
        '''pull an existing memory trace'''
        self.w("CALC:PAR:SEL 'CH1_S21_1'")
        self.w("CALCulate1:FORMat SMITh")
        self.w("FORMat ASCii,0")
        mem=self.q("CALCulate1:DATA? FMEM")
        return str(mem)

        
    def getTrace(self):
        '''Pull data from the current VNA trace. - JF '''
        measName,Sparam = str(self.q("CALC:PAR:CAT?")).strip('\"').split(',')
        self.w("CALC:PAR:SEL \'{}\'".format(measName))#select the measurement
        self.w("SENS:SWE:MODE SINGl; ABOR")
        #self.w("SENS1:SWE:MODE CONT")
        self.w("CALCulate1:FORMat SMITh")
        #self.w("FORMat REAL,32")
        self.w("FORMat ASCii,0")
        results = self.q("CALCulate1:DATA? FDATA")
        #self.inst.write_raw("CALCulate1:DATA? FDATA")
        #results = self.inst.read_raw()
        #results = self.inst.query_binary_values("CALCulate1:DATA? FDATA")
        return results
       