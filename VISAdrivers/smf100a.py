'''
This is a driver for the Rohde & Schwarz SMF100A signal generator.

'''
from pyvisa import ResourceManager
class SMF():
    
    def __init__(self,address,reset=False):
        '''creates instance of the class, shorthand for commands, and optionally
           resets the device. - JF'''
        self.inst = ResourceManager().open_resource(address)
        self.w = self.inst.write
        self.r = self.inst.read
        self.q = self.inst.query
        if reset:
            self.w("*RST")
        #self.w("DISP:PSAV ON")
        #self.w("DISP:UPD OFF")

    def screen_save(self, state):
        '''
        changes status of screen saver option.
        argument should be either "ON" or "OFF"
        -JF'''
        comm = "DISP:PSAV %s" % state
        self.w(comm)

    def id(self):
        '''queries instrument identification info -JF'''
        return self.q("*IDN?")

    def rst(self):
        '''resets instrument state to defaults -JF'''
        self.w("*RST")

    def set_freq(self, value):
        '''sets the RF frequency in GHz.
        Example: set_freq(10) sets the signal generator to 10 GHz.
        -JF'''
        comm = "FREQ %s GHz" % value
        self.w(comm)

    def RF_ON(self):
        '''enables RF output -JF'''
        self.w("OUTP ON")

    def RF_OFF(self):
        '''disables RF output -JF'''
        self.w("OUTP OFF")

    def screenshot(self):
        '''streams a hardcopy of the display to computer as JPG -JF'''
        self.w("HCOP:DEV:LANG JPG")
        return self.q("HCOP:DATA?")

    def save_state(self, num):
        '''saves the parameters and state of the signal generator to recall later.
        Example: save_state(4) saves the state to be recalled with the same # 
        -JF'''
        comm = "*SAV %s" % num
        self.w(comm)

    def recall_state(self, num):
        '''loads the settings which were associated to a number
        by the save_state() command. Example: recall_state(4) -JF'''
        comm = "*RCL %s" % num
        self.w(comm)

    def configure_AM(self, channel, source, depth):
        '''sets the modulation depth and source for specified channel.
        channel: which output channel do you want to configure? 1 or 2.
        source: specify the modulation source signal. Choose from {LF1, LF2,
        EXT1, EXT2, NOISe}
        depth: specify modulation depth in percent. 0 to 100 in 0.1 increment
        Example: configure_AM(1,EXT2,44.1) - JF'''
        comm1 = "AM%s:DEPT %s" % (channel, depth)
        comm2 = "AM%s:SOUR %s" % (channel, source)
        self.w(comm1)
        self.w(comm2)

    def enable_AM(self, channel):
        ''' turns on Amplitude Modulation for specified channel.
        Example: enable_AM(2) -JF'''
        comm = "AM%s:STAT ON" % channel
        self.w(comm)

    def disable_AM(self, channel):
        ''' turns off Amplitude Modulation for specified channel.
        Example: disable_AM(2) -JF'''
        comm = "AM%s:STAT OFF" % channel
        self.w(comm)

    def configure_freq_sweep(self, start, stop, step, stepUnits):
        '''sets the parameters for a frequency sweep.
           start: starting frequency in GHz. minimum 1 GHz.
           stop: end frequency in GHz. maximum 22 GHz.
           step: increment frequency. specify units in stepUnits.
           stepUnits: string specifying increment units.
           dwell: time per increment. specify units with dwellUnits.
           dwellUnits: string specifying units of dwell time.
           Example: configure_freq_sweep(5, 15, 500, "MHz", 12, "ms")
           sets up a sweep from 5 to 15 GHz with increments of 500 MHz
           and dwell time of 12 ms -JF'''
        #self.w("FREQ:MODE SWE")
        comm1 = "FREQ:STAR %s GHz" % (start)
        comm2 = "FREQ:STOP %s GHz" % (stop)
        comm3 = "SWE:STEP %s %s" % (step, stepUnits)
        #comm4 = "SWE:DWEL %s %s" % (dwell, dwellUnits)
        self.w(comm1);self.w(comm2);self.w(comm3);#self.w(comm4)
        self.w("FREQ:MODE SWE")
    
    def configure_trig_freq_sweep(self, source, mode):
        '''source: string in {IMMediate, BUS, EXTernal, EAUTo}
           mode: string in {AUTO, MANual, STEP}
           -JF'''
        comm1 = "TRIG:FSW:SOUR %s" % source
        comm2 = "SWE:MODE %s" % mode
        self.w(comm1);self.w(comm2)
    
    def configure_trig_power_sweep(self, source, mode):
        '''source: string in {IMMediate, BUS, EXTernal, EAUTo}
           mode: string in {AUTO, MANual, STEP}
           -JF'''
        comm1 = "TRIG:PSW:SOUR %s" % source
        comm2 = "SWE:POW:MODE %s" % mode
        self.w(comm)

    def set_level(self, level):
        '''sets the RF output level in dBm.
        Example: set_level(-30) sets the RF out to -30 dBm. -JF'''
        # note that this is subject to an offset which can be defined with SCPI.
        # the offset is zero by default. To set the level at RF out regardless of
        # offset, use "POW:POW %s" % level
        comm = "POW %s" % level
        self.w(comm)

    def configure_power_sweep(self, start, stop, step, dwell, dwellUnits):
        '''confiture_power_sweep(start, stop, step)
        configures the start, stop, and step levels for a power sweep. All values
        in units dBm.
        Example: configure_power_sweep(-35.5,-29.2,0.1)
        -JF'''
        #self.w("POW:MODE SWE")
        comm1 = "POW:STAR %s" % start
        comm2 = "POW:STOP %s" % stop
        comm3 = "SWE:POW:STEP %s" % step
        comm4 = "SWE:POW:DWEL %s %s" % (dwell, dwellUnits)
        self.w(comm1);self.w(comm2);self.w(comm3);self.w(comm4)

    def execute_freq_sweep(self):
        ''' begins the frequency sweep.
        To set parameters, use configure_freq_sweep. -JF'''
        #self.w("TRIG:FSW:SOUR SING")
        #self.w("SOUR:SWE:FREQ:MODE AUTO")
        #self.w("FREQ:MODE SWE")
        self.w("SWE:FREQ:EXEC")

    def execute_power_sweep(self):
        ''' begins the power level sweep.
        To set parameters, use configure_power_sweep. -JF'''
        self.w("SWE:POW:MODE STEP")
        self.w("SWE:POW:EXEC")

    def reset_all_sweeps(self):
        '''Resets all active sweeps to the starting point -JF'''
        self.w("SWE:RES")

    def error_check(self):
        ''' checks for any error codes in the instrument error queue -JF'''
        return self.q("SYST:ERR:ALL?")
        
    def display_update(self,boolean):
        '''sets whether or not to update the display, disable for fast sweeps -JF'''
        self.w("SYST:DISP:UPD {}".format(boolean))

