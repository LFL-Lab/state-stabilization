'''
A class for controlling the Tektronix AWG520. 
Contributors: James Farmer
'''

#test

# TODO:
# add data validation for all input parameters.


from pyvisa import ResourceManager
class AWG():

    def __init__(self,address,reset=False):
        '''This function is automatically called on initialization of an instance
           of the AWG class. 
           -JF'''
        self.inst = ResourceManager().open_resource(address)
        self.inst.write_termination = '\n'
        self.inst.read_termination = '\n'
        # Some shorthand used throughout.
        self.w = self.inst.write
        self.r = self.inst.read
        self.q = self.inst.query
        if reset:
            self.w('*RST')

    def ID(self):
        '''returns a string identifying the AWG.
           -JF'''
        return self.q('*idn?')

    def whatFiles(self):
        '''returns the contents of the AWG file storage.
           -JF'''
        return self.q('mmem:cat?')
    
    def reset(self):
        '''resets the AWG, restoring all parameters to defaults, etc.
           -JF'''
        self.w('*rst')
        
    def send_waveform(self,filename):
        '''Sends a premade waveform file (.wfm) to the AWG520.
           Example: send_waveform("E:/waveforms/example.wfm")
           -JF'''
        from lflPython.tools.wfmtools import wfm2datablock
        data, awgFile = wfm2datablock(filename)
        s1 = 'MMEM:DATA '.encode() + '"'.encode() + awgFile.encode() + '",#'.encode()
        print(s1)
        lenlen = str(len(str(len(data)))).encode('utf-8')
        strlen = str(len(data)).encode('utf-8')
        mes = s1 + lenlen + strlen + data
        print(mes)
        self.inst.write_raw(mes)
        return awgFile
    
    def load_waveform(self, channel, filename):
        '''Loads the specified file (which should be on AWG) to
           the given channel.
           -JF'''
        mes = 'SOUR%s:FUNC:USER "%s"' % (channel,filename)
        self.w(mes)
        
    def set_mode(self, mode):
        ''' sets the run mode of AWG, such as continuous or triggered, etc.
           mode can be one of {"CONT","TRIG","GAT","ENH"}.
           CONTinuous: continuously outputs the waveform. trigger has no effect.
           TRIGgered: one waveform cycle per trigger event.
           GATed: continuously outputs waveform while trigger is high.
           ENHanced: follows the run mode specified in the loaded sequence file.
           -JF'''
        if mode in ["CONT","TRIG","GAT","ENH"]:
            mes = 'AWGC:RMOD %s' % mode
            self.w(mes)
        else: return ValueError
    
    def clear_waveforms(self):
        '''clears the waveform on both channels.
           -JF'''
        self.w('sour1:func:user ""')
        self.w('sour2:func:user ""')
        
    def marker_level(self,channel,m1_low,m1_high,m2_low,m2_high):
        '''sets high and low levels of markers 1 and 2 for the
           given channel.
           -JF'''
        mes1 = 'SOUR%s:MARK1:HIGH %s' % (channnel, m1_high)
        mes2 = 'SOUR%s:MARK2:HIGH %s' % (channnel, m2_high)
        mes3 = 'SOUR%s:MARK1:LOW %s' % (channnel, m1_low)
        mes4 = 'SOUR%s:MARK2:LOW %s' % (channnel, m2_low)
        mes = mes1 + ';' + mes2 + ';' + mes3 + ';' + mes4
        self.w(mes)
        
    def signal_level(self, channel, volts, offset=False):
        '''sets the peak to peak voltage of specified channel to
           given level in volts. Optionally, a voltage offset can
           be specified as well.
           -JF'''
        mes = 'SOUR%s:VOLT %s' % (channel,volts)
        self.w(mes)
        if offset:
            mes2 = 'SOUR%s:VOLT:OFFS %s' % (channel,offset)
            self.w(mes2)
        
    def BEEP(self):
        '''BEEEEEP.
           -JF'''
        self.w('SYST:BEEP')
    
    def lock_interface(self):
        '''locks the front panel and buttons.
           -JF'''
        self.w('SYST:KLOC ON')
    
    def unlock_interface(self):
        '''unlocks the front panel and buttons.
           -JF'''
        self.w('SYST:KLOC OFF')
        
    def reset_trigger(self):
        '''resets the trigger system and places all trigger 
           sequences in idle state.
           -JF'''
        self.w('ABOR')
    
    def force_trigger(self):
        '''generates an immediate trigger.
           -JF'''
        self.w('*TRG')
    
    def trigger_level(self,level=1.4):
        '''sets trigger level in volts. Values range from 
           -5.0 to 5.0 V in steps of 0.1 V.
           -JF'''
        mes = 'TRIG:LEV %s' % level
        self.w(mes)
    
    def external_trigger_impedance(self, impedance=50):
        '''sets the impedance of external trigger input.
           impedance can be 50 Ohm or 1000 Ohm.
           -JF'''
        mes = 'TRIG:IMP %s' % impedance 
        self.w(mes)
    
    def enable_output(self,channel):
        '''sets the AWG DAC to run and enables the output on the specified channel.
           -JF'''
        mes = 'OUTP%s ON' % channel
        self.w('AWGC:RUN')
        self.w(mes)
        
    def disable_output(self,channel):
        '''turns the specified channel output off and stops the AWG DAC. -JF'''
        mes = 'OUTP%s OFF' % channel
        self.w(mes)
        self.w('AWGC:STOP')
