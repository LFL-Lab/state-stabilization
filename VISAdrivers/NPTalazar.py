''' A class for controlling the Alazar ATS9371 diitizer. When using this class,
call configureClock and configureTrigger functions to configure the board
before acquiring data. -JF'''

import atsapi as ats
from ctypes import c_uint16
from LFLpd.tools.datatools import unique_filename

class ADC():
    
    def __init__(self):
        '''Initialize the board and set the channel input parameters,
           doubtful that we'll even need to change the coupling. -JF'''
        self.board = ats.Board(systemId=1,boardId=1)
        self.board.inputControlEx(ats.CHANNEL_A,
                                  ats.DC_COUPLING,
                                  ats.INPUT_RANGE_PM_400_MV,
                                  ats.IMPEDANCE_50_OHM)
        self.board.inputControlEx(ats.CHANNEL_B,
                                  ats.DC_COUPLING,
                                  ats.INPUT_RANGE_PM_400_MV,
                                  ats.IMPEDANCE_50_OHM)
    ##########################
    # Configuration functions must be called before data acquisition
    ##########################
    def configureClock(self,MS_s=1000,INT_EXT="EXT"):
        ''' Configures the sampling rate and clock source.
           MS_s is sampling rate in mega-samples per second. Integer 300-1000.
           INT_EXT sets source. String "EXT" or "INT". default="EXT"
           -JF'''
        self.samplesPerSec = MS_s * 1e6
        if INT_EXT == "EXT":
            self.board.setCaptureClock(ats.EXTERNAL_CLOCK_10MHz_REF,
                                       self.samplesPerSec,
                                       ats.CLOCK_EDGE_RISING,
                                       1)
        elif INT_EXT == "INT":
            self.board.setCaptureClock(ats.INTERNAL_CLOCK,
                                       self.samplesPerSec,
                                       ats.CLOCK_EDGE_RISING,
                                       1)
        else:
            return ValueError
    
    def configureTrigger(self,source="EXT",trig_level_mV=20,delay=0,timeout=0):
        ''' Configures the board for specified trigger source and level.
           source should be one of {"A","B","EXT"} where A and B are channels
           trig_level_mV sets threshold voltage (in mV) for channel trigger.
           delay,timeout in seconds. default: no delay, no timeout.
           -JF'''
        if source == "A":
            atsSource = ats.TRIG_CHAN_A
            level = 128 + int(127 * trig_level_mV/400 + 0.5)
        elif source == "B":
            atsSource = ats.TRIG_CHAN_B
            level = 128 + int(127 * trig_level_mV/400 + 0.5)
        else:
            atsSource = ats.TRIG_EXTERNAL
            level = 150
            self.board.setExternalTrigger(ats.DC_COUPLING,
                                          ats.ETR_TTL)
        
        if level > 0 and level < 255:
            self.board.setTriggerOperation(ats.TRIG_ENGINE_OP_J,
                                           ats.TRIG_ENGINE_J,
                                           atsSource,
                                           ats.TRIGGER_SLOPE_POSITIVE,
                                           level,
                                           ats.TRIG_ENGINE_K,
                                           ats.TRIG_DISABLE,
                                           ats.TRIGGER_SLOPE_POSITIVE,
                                           128)
        else:
            return ValueError
        
        triggerDelay_samples = int(delay * self.samplesPerSec + 0.5)
        self.board.setTriggerDelay(triggerDelay_samples)
        
        triggerTimeout_clocks = int(timeout / 10e-6 + 0.5)
        self.board.setTriggerTimeOut(triggerTimeout_clocks)
        
        #sets the aux IO output to TTL high when acquiring data
        self.board.configureAuxIO(ats.AUX_OUT_TRIGGER,0)
    
    def startTriggeredCapture(self,samples,channel="A",returnfname=False):
        #acquisitionLength_sec = duration
        #samplesPerBuffer = 128
        preTriggerSamples = 0
        postTriggerSamples = samples
        recordsPerBuffer = 10
        buffersPerAcquisition = 10
        acquisitionLength_sec = float(samples/self.samplesPerSec)
        
        if channel == "A":
            channels = ats.CHANNEL_A
        elif channel == "B":
            channels = ats.CHANNEL_B
        elif channel == "AB":
            channels = ats.CHANNEL_A | ats.CHANNEL_B
        else:
            return ValueError
        
        channelCount = 0
        for c in ats.channels:
            channelCount += (c & channels == c)
        
        directory="E:/rawData/"
        dataFileName = unique_filename(directory,prefix='rawData',ext='bin')
        dataFile = open(dataFileName,'wb')
        
        # Save timestamp and parameters 

        with open(dataFileName[0:-4] + ".txt",'w') as f:
            from time import strftime
            f.write(strftime("%c")+'\n')
            f.write("Channels: " + channel + '\n')
            f.write("Acquisition duration: " + str(acquisitionLength_sec) + " seconds." + '\n')
            f.write("Samples per second: " + str(self.samplesPerSec*1e-6) + '\n')
        
        # compute bytes per record and per buffer
        memorySize_samples, bitsPerSample = self.board.getChannelInfo()
        bytesPerSample = (bitsPerSample.value + 7) // 8
        #bytesPerBuffer = bytesPerSample * samplesPerBuffer * channelCount
        samplesPerRecord = preTriggerSamples + postTriggerSamples
        bytesPerRecord = bytesPerSample * samplesPerRecord
        bytesPerBuffer = bytesPerRecord * recordsPerBuffer * channelCount
        
        #Calculate the number of buffers in the acquisition
        #samplesPerAcquisition = int(self.samplesPerSec * acquisitionLength_sec + 0.5)
        #buffersPerAcquisition = ((samplesPerAcquisition + samplesPerBuffer - 1) //
        #                         samplesPerBuffer)
        
        # TODO: Select number of DMA buffers to allocate
        bufferCount = buffersPerAcquisition
        
        # Allocate DMA buffers
        sample_type = c_uint16
        buffers = []
        for i in range(bufferCount):
            buffers.append(ats.DMABuffer(self.board.handle, sample_type, bytesPerBuffer))
        
        # Set the record size
        self.board.setRecordSize(preTriggerSamples, postTriggerSamples)
        recordsPerAcquisition = recordsPerBuffer * buffersPerAcquisition
        
        self.board.beforeAsyncRead(channels,
                                   -preTriggerSamples,  
                                   samplesPerRecord,
                                   recordsPerBuffer,       
                                   recordsPerAcquisition,        # recordsPerAcquisition
                                   ats.ADMA_EXTERNAL_STARTCAPTURE | ats.ADMA_NPT | ats.ADMA_FIFO_ONLY_STREAMING)
                                   
        # Post DMA buffers to board
        for buffer in buffers:
            self.board.postAsyncBuffer(buffer.addr, buffer.size_bytes)
        
        
        try:
            self.board.startCapture() # Start the acquisition
            buffersCompleted = 0
            bytesTransferred = 0
            bitShift = 4
            codeRange = (1<<(bitsPerSample.value - 1)) - 0.5
            while (buffersCompleted < buffersPerAcquisition and not
                   ats.enter_pressed()):
                # Wait for the buffer at the head of the list of available
                # buffers to be filled by the board.
                buffer = buffers[buffersCompleted % len(buffers)]
                self.board.waitAsyncBufferComplete(buffer.addr, timeout_ms=5000)
                buffersCompleted += 1
                #bytesTransferred += buffer.size_bytes
                
                sampleValues = buffer.buffer >> bitShift
                sampleValues.tofile(dataFile)
                self.board.postAsyncBuffer(buffer.addr, buffer.size_bytes)
        finally:
            self.board.abortAsyncRead()
            dataFile.close()
        if returnfname:
            return dataFileName