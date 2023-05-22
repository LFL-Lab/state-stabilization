#!/usr/bin/env python
# coding: utf-8

# This package includes
# * General setting
# * UHFQA connection
# * UHFQA general setting
# * UHFQA QA setting
# * UHFQA AWG setting
# * Scope module
# * Plot

# # General setting

# ## Import modules

# In[1]:


device_id_uhf = 'DEV2528'

import matplotlib.pyplot as plt
import time
import textwrap
import numpy as np
import scipy as sp
import enum
import scipy as scy
import zhinst.utils as ziut
import zhinst
from zhinst import ziPython

# ## Plot setting

# In[2]:


pi=np.pi
SMALL_SIZE = 16
MEDIUM_SIZE = 16
BIGGER_SIZE = 16

# plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
# plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
# plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
# plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
# plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
# plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
# plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


# ## ResultLoggingSource

# In[41]:


class ResultLoggingSource(enum.IntEnum):
    '''
    Constants for selecting result logging source
    '''
    TRANS = 0
    THRES = 1
    ROT = 2
    TRANS_STAT = 3
    CORR_TRANS = 4
    CORR_THRES = 5
    CORR_STAT = 6
    integration = 7


# # UHFQA connection


def create_api_sessions_uhf(device_uhf_id, use_discovery = 1, ip = '10.42.0.226'):
    '''
    create API sessions for UHFQA

    device_id_uhf:    device ID of UHF
    '''
#     global daq, device_uhf
    apilevel_example = 6  # The API level supported by this example.

    # Call a zhinst utility function that returns:
    # - an API session `daq` in order to communicate with devices via the data server.
    # - the device ID string that specifies the device branch in the server's node hierarchy.
    # - the device's discovery properties.
    # required_devtype = 'UHFQA'
    # required_options = ['QA','AWG']

    if not use_discovery:
        daq = zhinst.ziPython.ziDAQServer(ip, 8004, apilevel_example)
        daq.connectDevice(device_uhf_id, '1gbe')
        device_uhf = device_uhf_id
        daq.setInt(f'/zi/config/open', 1)
    else:
        daq, device_uhf,props = ziut.create_api_session(device_uhf_id, apilevel_example)
    return daq, device_uhf

def init(daq, device, output_range_uhf = 1.5, input_range_uhf = 0.3):
    '''
    Initialize device for UHFQA examples

    daq:               device daq
    device:            device ID
    output_range_uhf:  ouput range of UHF
    input_range_uhf:   input range of UHFQA
    '''
    # General setup
    parameters = [
        # Input and output ranges
        ('sigins/*/range', input_range_uhf),
        # ('sigouts/*/range', output_range_uhf),
        # Set termination to 50 Ohm
        ('sigins/*/imp50', 1),
        ('sigouts/*/imp50', 1),
        # Turn on both outputs
        ('sigouts/*/on', 1),
        # AWG in direct mode
        ('awgs/*/outputs/*/mode', 0),
        # DIO:
        # - output AWG waveform as digital pattern on DIO connector
        # ('dios/0/mode', 2),
        # - drive DIO bits 15 .. 0
        # ('dios/0/drive', 2),
        # Delay:
        # ('qas/0/delay', 0),
        # Deskew:
        # - straight connection: sigin 1 -- channel 1, sigin 2 -- channel 2
        # ('qas/0/deskew/rows/0/cols/0', 1),
        # ('qas/0/deskew/rows/0/cols/1', 0),
        # ('qas/0/deskew/rows/1/cols/0', 0),
        # ('qas/0/deskew/rows/1/cols/1', 1),
        # Results:
        # ResultLoggingSource.TRANS is crosstalk?
        # ('qas/0/result/length', 1.0),
        # ('qas/0/result/averages', 0),
        # ('qas/0/result/source', ResultLoggingSource.TRANS),
        # Statistics:
        # ('qas/0/result/statistics/length', 1.0),
        # Monitor length:
        # ('qas/0/monitor/length', 1024)
    ]

    # Number of readout channels
    num_readout_channels = 10

    # Rotation
    for i in range(num_readout_channels):
        parameters.append(('qas/0/rotations/{:d}'.format(i), 1+0j))

    # Transformation
    # - no cross-coupling in the matrix multiplication (identity matrix)
    for i in range(num_readout_channels):
        for j in range(num_readout_channels):
            parameters.append(('qas/0/crosstalk/rows/{:d}/cols/{:d}'.format(i, j), int(i == j)))

    # Threshold
    for i in range(num_readout_channels):
        parameters.append(('qas/0/thresholds/{:d}/level'.format(i), 1.0))

    # Update device
    daq.set([('/{:s}/{:s}'.format(device, node), value) for node, value in parameters])
    daq.sync()


def setup_ssb(daq,rr_IF):

    daq.setDouble('/dev2528/oscs/0/freq', rr_IF)
    daq.setInt('/dev2528/awgs/0/outputs/0/mode', 1)
    daq.setInt('/dev2528/awgs/0/outputs/1/mode', 1)

def setup_mixer_calib(daq):

    assert daq.get('/dev2528/awgs/0/outputs/0/mode')['dev2528']['awgs']['0']['outputs']['0']['mode']['value'][0] and daq.get('/dev2528/awgs/0/outputs/0/mode')['dev2528']['awgs']['0']['outputs']['0']['mode']['value'][0] and daq.get('/dev2528/oscs/0/freq')['dev2528']['oscs']['0']['freq']['value'][0] != 0,  "You need to setup SSB first"

    awg_program=textwrap.dedent("""

        const N = 32000;
        const a = 1.0;
        const pos = N/2;
        const w = N/4;

        // Readout pulse
        wave readoutPulse = gauss(N,a,pos,w);

        while(true) {
            resetOscPhase();
            startQA();
            playWave(readoutPulse,readoutPulse);
        }

    """)
    create_and_compile_awg(daq,  awg_program, seqr_index = 0, timeout = 1)

def awg_seq_readout(daq, cav_resp_time = 4e-6,base_rate = 450e6, amplitude_uhf = 1,rr_IF=5e6,readout_length = 2.2e-6,nPoints=1000,timeout=1):
    '''
    Create AWG sequence for spectroscopy experimentm compile it and upload it

    daq:                device daq
    device:             device ID
    sample_exponent:    a factor (base rate / 2^(sample_exponent)) to reduce sampling rate. Default is 450 MHz
    base_rate:          sampling rate
    amplitude_uhf:      amplitude setting in SeqC
    readout_length:     readout length in second
    averages_exponent:  average times in 2^(averages_exponent)
    result_length:      number of interested result
    '''

    awg_program=textwrap.dedent("""
    const fs = _c00_;
    const f_c = 1.8e9;      // clock rate
    const f_seq = f_c/8;     // sequencer instruction rate
    const dt = 1/f_seq;
    const cav_resp_time = _c1_;

    // Readout pulse
    wave readoutPulse = _c3_*ones(_c2_);

    while(true) {
        waitDigTrigger(1,1);
        startQA();
        playWave(1,readoutPulse);
    }
    """)
    awg_program = awg_program.replace('_c00_', str(base_rate))
    awg_program = awg_program.replace('_c1_', str(int(cav_resp_time*base_rate)))
    awg_program = awg_program.replace('_c2_', str(round(readout_length*base_rate)))
    awg_program = awg_program.replace('_c3_', str(amplitude_uhf))

    daq.setInt('/dev2528/awgs/0/auxtriggers/0/channel', 2) # sets the source of digital trigger 1 to be the signal at trigger input 3 (back panel)
    daq.setDouble('/dev2528/triggers/in/2/level', 0.1)
    daq.setInt('/dev2528/awgs/0/auxtriggers/0/slope', 1)
    create_and_compile_awg(daq,awg_program, seqr_index = 0, timeout = timeout)

def config_qa(daq,integration_length=2.2e-6,delay=300e-9,nAverages=128,rr_IF=5e6,sequence='pulse',result_length=1,source=7):
    # print('-------------Configuring QA-------------\n')
    base_rate=1.8e9
    bypass_crosstalk=0
    # set modulation frequency of QA AWG to some IF and adjust input range for better resolution
    daq.setDouble('/dev2528/sigins/0/range',0.5)
    # daq.setInt('/dev2528/oscs/0/freq'.format(device),int(rr_IF))
    # daq.setInt('/dev2528/awgs/0/outputs/0/mode', 1) # AWG set to modulation mode

    # QA setup Settings
    daq.setInt('/dev2528/qas/0/integration/sources/0', 0)
    daq.setInt('/dev2528/qas/0/delay',round(delay*base_rate))
    # if sequence =='spec':
    #     daq.setInt('/dev2528/qas/0/integration/mode'.format(device), 0) # 0 for standard (4096 samples max), 1 for spectroscopy
    # elif sequence =='pulse':
    #     daq.setInt('/dev2528/qas/0/integration/mode'.format(device), 0)
    daq.setInt('/dev2528/qas/0/integration/mode', 0) # 0 for standard (4096 samples max), 1 for spectroscopy
    daq.setInt('/dev2528/qas/0/integration/length', int(base_rate*integration_length))
    daq.setInt('/dev2528/qas/0/bypass/crosstalk', bypass_crosstalk)   #No crosstalk matrix
    daq.setInt('/dev2528/qas/0/bypass/deskew', 1)   #No crosstalk matrix
    # daq.setInt('/dev2528/qas/0/bypass/rotation'.format(device), 1)   #No rotation
    # x = np.linspace(0,integration_length,round(integration_length*base_rate))
    # weights_I = np.sin(2*np.pi*rr_IF*x)
    # weights_Q = np.zeros(round(integration_length*base_rate))
    weights_I = weights_Q = np.ones(round(integration_length*base_rate))
    # weights_Q = np.zeros(round(integration_length*base_rate))
    daq.setVector('/dev2528/qas/0/integration/weights/0/real', weights_I)
    daq.setVector('/dev2528/qas/0/integration/weights/0/imag', weights_Q)
    daq.sync()
    daq.setInt('/dev2528/qas/0/integration/trigger/channel', 7); # 0 for trigger input ch 1, 7 for internal triggering (QA AWG -> QA)

    # QA Monitor Settings
    daq.setInt('/dev2528/qas/0/monitor/trigger/channel', 7)
    daq.setInt('/dev2528/qas/0/monitor/averages',nAverages)
    daq.setInt('/dev2528/qas/0/monitor/length', 4096)
    # configure triggering (0=trigger input 1 7 for internal trigger)

    # QA Result Settings
    daq.setInt('/dev2528/qas/0/result/length', result_length)
    daq.setInt('/dev2528/qas/0/result/averages', nAverages)
    daq.setInt('/dev2528/qas/0/result/source', source) # 2 -> source = rotation | 7 = integration
    daq.setInt('/dev2528/qas/0/result/reset', 1)
    daq.setInt('/dev2528/qas/0/result/enable', 1)
    daq.setInt('/dev2528/qas/0/result/mode',0) # cyclic averaging
    daq.sync()


def config_sigin(daq, device, input_range_uhf = 1, ac = [0,0]):
    '''
    config signal input

    daq,                  daq ID
    device:               device ID
    input_range_uhf: input range of UHF, 10 mV to 1.5 V
    ac:   (0/1, 0/1):     (not AC couple/AC coupled for input ch1,  not AC couple/AC coupled for input ch2)
    '''

    daq.setInt('/%s/sigins/0/imp50' % device, 1)
    daq.setInt('/%s/sigins/1/imp50' % device, 1)
    daq.setInt('/%s/sigins/0/ac' % device, ac[0])
    daq.setInt('/%s/sigins/1/ac' % device, ac[1])
    daq.setDouble('/%s/sigins/0/range' % device, input_range_uhf)

def toggle_outputs(daq, device, channel=None):
    """
    Toggles signal output of UHFQA. If no channel specified toggles both.
    Arguments:
        daq (zhinst.ziDAQServer) -- Data Server Object
        device (String) -- device ID, e.g. "dev2266"
    Keyword Arguments:
        channel (int) -- list of channels to be toggled (in 0, 1)
                            (default: None)
    """

    if channel is not None:
        assert channel in [0, 1]
        channel = [channel]
    else:
        channel = [0, 1]
    for ch in channel:
        path = f"/{device}/sigouts/{ch}/on"
        if daq.getInt(path) == 0:
            daq.setInt(path, 1)
        elif daq.getInt(path) == 1:
            daq.setInt(path, 0)

def enable_output(daq, device, channel=[0, 1], enable=[1, 1]):
    """
    Enable or disable output channels.
    Arguments:
        daq (zhinst.ziDAQServer) -- Data Server Object
        device (String) -- device ID, e.g. "dev2266"
    Keyword Arguments:
        channel (array) -- channel array
                            (default: None)
        enable  (array) -- enable or disable arrage according to the channels. Default is 1 -> enable, 0 -> disable

    """


    for i in range(len(channel)):
        daq.setInt(f"/{device}/sigouts/{channel[i]}/on", enable[i])

def osc_f(daq, device, frequency=150e6):
    '''
    oscillator frequency setting

    daq:                         daq ID
    device:                      device ID
    frequency:      frequency of osc in Hz
    '''
    daq.setDouble('/{:s}/oscs/0/freq'.format(device), frequency)

def triggers(daq, device, in_out=['out','out'], source=[74,64], dio_ch = 0, dio_drive = 1):
    '''
    trigger sources setting and trigger in or out enabled for trigger port 1 and 2

    daq:                                daq ID
    device:                             device ID
    in_out:   ['in'/'out', 'in'/'out']: in or out enabled for trigger port 1 and 2
    source:            64:              QA result 1
                       74:              QA result Trigger
    '''
    daq.setInt(f'/{device:s}/dios/{dio_ch:d}/drive', dio_drive);
    for i in range(len(in_out)):
        daq.setInt('/{:s}/triggers/%s/%d/source'.format(device) % (in_out[i], i), source[i])   #Trigger out 1 to "QA Result Trigger"
        daq.setInt('/{:s}/triggers/%s/%d/drive'.format(device)  % (in_out[i], i), 1)

def set_qa_monitor(daq, device, monitor_length, averages):
    """
    Applies settings to the QA Monitor tab.
    Arguments:
        daq (zhinst.ziDAQServer) -- Data Server Object
        device (String) -- device ID, e.g. "dev2266"
        monitor_length (int) -- number of samples recorded in monitor tab
        averages (int) -- number of averages for monitor tab
    """

    settings = [
        ("qas/0/monitor/enable", 0),
        ("qas/0/monitor/length", monitor_length),
        ("qas/0/monitor/averages", averages),
        ("qas/0/monitor/enable", 1),
        ("qas/0/monitor/reset", 1),
    ]
    daq.set([(f"/{device}/{node}", value) for node, value in settings])

def set_qa_results(daq, device, result_length, result_averages, source="integration"):
    """
    Applies settings to the QA Results tab.
    Arguments:
        daq (zhinst.ziDAQServer) -- Data Server Object
        device (String) -- device ID, e.g. "dev2266"
        result_length (int) --  number of samples to be recorded
        result_averages (int) -- number of averages for results
    Keyword Arguments:
        source (str) -- specifies data source of QA
                        "integratio", "rotation" or "threshold"
                        (default: "integration")
    """

    if source == "integration":
        source = 7
    elif source == "rotation":
        source = 2
    elif source == "threshold":
        source = 1

    settings = [
        ("qas/0/result/enable", 0),
        ("qas/0/result/length", result_length),
        ("qas/0/result/averages", result_averages),
        ("qas/0/result/source", source),
        ("qas/0/result/enable", 1),
        ("qas/0/result/reset", 1),
    ]
    time.sleep(0.01)
    daq.set([(f"/{device}/{node}", value) for node, value in settings])

def qa_result_reset(daq, device, reset = 1):
    '''
    result reset of QA

    device:            device ID
    reset:      0:     do not reset
    reset:      1:     reset
    '''

    daq.setInt('/{:s}/qas/0/result/reset'.format(device), reset)

def qa_result_enable(daq, device, enable = 1):
    '''
    enable QA result

    device:            device ID
    enable:      0:    disable QA result
                 1:    enable QA result
    '''

    daq.setInt('/{:s}/qas/0/result/enable'.format(device), enable)

def create_sweep_data_dict(daq, device):
    '''
    create sweep data dictionary for ch1 and ch2

    device:         device ID
    '''

    channels = [0, 1]
    # Subscribe to result waves
    paths = []
    for ch in channels:
        path = '/{:s}/qas/0/result/data/{:d}/wave'.format(device, ch)
        paths.append(path)
    daq.subscribe(paths)
    sweep_data = dict()

    for path in paths:
        sweep_data[path] = np.array([])
    return sweep_data, paths

def acquisition_poll(daq, paths, num_samples, timeout):
    """ Polls the UHFQA for data. Taken from zhinst.examples.uhfqa.common
    Args:
        paths (list): list of subscribed paths
        num_samples (int): expected number of samples
        timeout (float): time in seconds before timeout Error is raised.
    """
    poll_length = 0.001  # s
    poll_timeout = 0  # ms
    poll_flags = 0
    poll_return_flat_dict = True

    # Keep list of recorded chunks of data for each subscribed path
    chunks = {p: [] for p in paths}
    gotem = {p: False for p in paths}

    # Poll data
    t = 0

    while t < timeout and not all(gotem.values()):

        dataset = daq.poll(poll_length, poll_timeout, poll_flags, poll_return_flat_dict)

        for p in paths:
            if p not in dataset:
                continue
            for v in dataset[p]:
                chunks[p].append(v['vector'])
                num_obtained = sum([len(x) for x in chunks[p]])
                if num_obtained >= num_samples:
                    gotem[p] = True
        t += poll_length

    if not all(gotem.values()):
        for p in paths:
            num_obtained = sum([len(x) for x in chunks[p]])
            print('Path {}: Got {} of {} samples'.format(p, num_obtained, num_samples))
        raise Exception('Timeout Error: Did not get all results within {:.1f} s!'.format(timeout))

    # Return dict of flattened data
    return {p: np.concatenate(v) for p, v in chunks.items()}

def stop_result_unit(daq, device, paths):
    '''
    stop QA result unit,

    daq:             dag ID
    device:          device ID
    paths:           data paths
    '''
    daq.unsubscribe(paths)
    daq.setInt('/{:s}/qas/0/result/enable'.format(device), 0)

def config_awg(daq, device, base_rate = 1.8e9, intermediate_frequency_uhf = 150e6,amplitude_scale_uhf = 1, input_range_uhf = 1.5, sample_exponent = 5,averages_exponent = 5, awg_mode = 1):

    '''
    config AWG

    daq:                         daq ID
    device:                      device ID
    base_rate:                   sampling rate
    intermediate_frequency_uhf:  UHF osc frequency in Hz
    amplitude_scale_uhf:         amplitude of AWG outupt 1 and 2, range from 0 to 1
    input_range_uhf:             input range of UHF, range from 10 mV to 1.5 V
    sample_exponent:             a factor (base rate / 2^(sample_exponent)) to reduce sampling rate
    averages_exponent:           averate times in 2^(averages_exponent)
    awg_mode:       0:           Plain: AWG Output goes directly to Signal Output.
                    1:           Modulation: AWG Output 1 (2) is multiplied with oscillator
                                 signal of demodulator 4 (8).
                    2:           Advanced: Output of AWG channel 1 (2) modulates demodulators
                                 1-4 (5-8) with independent envelopes.
    '''

    # Configfigure AWG

    daq.setDouble('/{:s}/awgs/0/outputs/0/amplitude'.format(device), amplitude_scale_uhf)
    daq.setDouble('/{:s}/awgs/0/outputs/1/amplitude'.format(device), amplitude_scale_uhf)
#    sample_rate = base_rate/2**sample_exponent
#    waveform_length = sample_rate*readout_length
#    num_averages = 2**averages_exponent
#    awg_seq(c0=sample_exponent, c1=readout_length, c20=result_length, c21=averages_exponent, c3=amplitude_uhf, contn=continuous)

    # Configure AWG
    daq.setInt('/{:s}/awgs/0/time'.format(device), sample_exponent)
    daq.setInt('/{:s}/awgs/0/outputs/*/mode'.format(device), awg_mode) # modulation mode =1
    daq.setInt('/{:s}/awgs/0/auxtriggers/0/slope'.format(device), 1)
    # Adjust input range
    daq.setDouble('/{:s}/sigins/*/range'.format(device), input_range_uhf)

    daq.asyncSetInt('/{:s}/awgs/0/single'.format(device), 1)
    # Now we're ready for readout. Enable result unit and start acquisition.
#    daq.setInt('/{:s}/qas/0/result/reset'.format(device), 1)
#    daq.setInt('/{:s}/qas/0/result/enable'.format(device), 1)
    daq.setDouble('/{:s}/oscs/0/freq'.format(device), intermediate_frequency_uhf)
    daq.asyncSetInt('/{:s}/awgs/0/single'.format(device), 1)

def awg_enable(daq, device, enable = 1):
    '''
    enable AWG

    device:            device ID
    enable:    0:      disable AWG
               1:      enable AWG
    '''
    daq.syncSetInt('/{:s}/awgs/0/enable'.format(device), enable)

def enable_awg(daq, device, enable=1):
    '''
    enable/disable AWG

    daq:             dag ID
    device:          device ID
    enable:     0:   diable AWG
                1:   enable AWG
    '''
    daq.asyncSetInt('/{:s}/awgs/0/single'.format(device), 1)
    daq.syncSetInt('/{:s}/awgs/0/enable'.format(device), enable)



def awg_seq_rb_1q(daq, device, trigger_delay_sec= 50e-9, wave_dur_sec=2e-6, amplitude_uhf=0.1):
    awg_program = textwrap.dedent("""
    const f_s = 1.8e9;
    const f_seq = f_s/8;

    const trigger_delay_sec = _c0_; // adjust this number if necessary to align the control and readout pulse
    const wave_dur_sec  = _c1_; // width of the readout pulse
    const amplitude_uhf = _c2_;

    wave w_I_qb1 = amplitude_uhf*ones(wave_dur_sec*f_s);
    wave w_Q_qb1 = amplitude_uhf*ones(wave_dur_sec*f_s);

    setTrigger(AWG_INTEGRATION_ARM);
    while (true) {
        waitDigTrigger(1, 1);
        wait(trigger_delay_sec*f_seq);
        setTrigger(AWG_INTEGRATION_ARM + AWG_INTEGRATION_TRIGGER +  AWG_MONITOR_TRIGGER);
        playWave(w_I_qb1, w_Q_qb1);
        setTrigger(AWG_INTEGRATION_ARM);
    }
    """)
    awg_program = awg_program.replace('_c0_', str(trigger_delay_sec))
    awg_program = awg_program.replace('_c1_', str(wave_dur_sec))
    awg_program = awg_program.replace('_c2_', str(amplitude_uhf))

    create_and_compile_awg(daq, device,  awg_program, seqr_index = 0, timeout = 2)

def create_and_compile_awg(daq, awg_program, seqr_index= 0, timeout=1,verbose=0):
    device = 'dev2528'
    awgModule = daq.awgModule()
    awgModule.set('device', device)
    awgModule.set('index', seqr_index)
    awgModule.execute()
    """Compile and upload awg_program as .elf file"""
    if verbose==0:
        print("Starting compilation.")
        awgModule.set('compiler/sourcestring', awg_program)
        compilerStatus = -1
        while compilerStatus == -1:
            compilerStatus = awgModule.getInt('compiler/status')
            time.sleep(0.1)
        compilerStatusString = awgModule.getString('compiler/statusstring')
        # print(f"Compiler messages:\n--------------\n{compilerStatusString}\n--------------")
        if compilerStatus == 1: # compilation failed
            print(awg_program)
            raise Exception("Compilation failed.")
        # if compilerStatus == 0:
        #     print("Compilation successful with no warnings.")
        # if compilerStatus == 2:
        #     print("Compilation successful with warnings.")
        # print("Waiting for the upload to the instrument.")
        elfProgress = 0
        elfStatus = 0
        lastElfProgressPrc = None
        while (elfProgress < 1.0) and (elfStatus != 1):
            elfProgress = awgModule.getDouble('progress')
            elfStatus = awgModule.getInt('elf/status')
            elfProgressPrc = round(elfProgress * 100);
            if elfProgressPrc != lastElfProgressPrc:
                # print(f'Upload progress: {elfProgressPrc:2.0f}%')
                lastElfProgressPrc = elfProgressPrc
            time.sleep(0.1)
    else:
        print("Starting compilation.")
        awgModule.set('compiler/sourcestring', awg_program)
        compilerStatus = -1
        while compilerStatus == -1:
            compilerStatus = awgModule.getInt('compiler/status')
            time.sleep(0.1)
        compilerStatusString = awgModule.getString('compiler/statusstring')
        print(f"Compiler messages:\n--------------\n{compilerStatusString}\n--------------")
        if compilerStatus == 1: # compilation failed
            raise Exception("Compilation failed.")
        if compilerStatus == 0:
            print("Compilation successful with no warnings.")
        if compilerStatus == 2:
            print("Compilation successful with warnings.")
        print("Waiting for the upload to the instrument.")
        elfProgress = 0
        elfStatus = 0
        lastElfProgressPrc = None
        while (elfProgress < 1.0) and (elfStatus != 1):
            elfProgress = awgModule.getDouble('progress')
            elfStatus = awgModule.getInt('elf/status')
            elfProgressPrc = round(elfProgress * 100);
            if elfProgressPrc != lastElfProgressPrc:
                print(f'Upload progress: {elfProgressPrc:2.0f}%')
                lastElfProgressPrc = elfProgressPrc
            time.sleep(0.1)
    if elfStatus == 0:
        print("Upload to the instrument successful.")
    if elfStatus == 1:
        raise Exception("Upload to the instrument failed.")

def awg_seq_reset_wo_qubit(daq, device):
    awg_program = textwrap.dedent("""    wave w_g = 0.1*ones(1024);
    wave w_e = 0.5*ones(1024);

    while (1) {

      if (getUserReg(0) == 0){
          playWave(w_g, w_g);
          startQAResult(1, 1);
          waitWave();
          wait(100);
          }
      else {
          playWave(w_e, w_e);
          startQAResult(1, 1);
          waitWave();
          wait(100);
         }
    }
       """ )
    create_and_compile_awg(daq, device,  awg_program, seqr_index = 0, timeout = 2)

def awg_seq_latency(daq, device):

    '''
    AWG sequence for latency demo of UHFQA only and compile it

    daq:                device daq
    device:             device ID
    '''

    awg_program = textwrap.dedent("""
    const length = 4096;
    wave w_I = ones(length);
    wave w_Q = ones(length);

    setTrigger(AWG_INTEGRATION_ARM);

    while(true) {
      repeat(100) {
        playWave(w_I, w_Q);
        setTrigger(AWG_INTEGRATION_ARM + AWG_INTEGRATION_TRIGGER + AWG_MONITOR_TRIGGER + 0b10);
        setTrigger(AWG_INTEGRATION_ARM);
        waitWave();
        wait(1024);

        playWave(-w_I, -w_Q);
        setTrigger(AWG_INTEGRATION_ARM + AWG_INTEGRATION_TRIGGER + AWG_MONITOR_TRIGGER + 0b10);
        setTrigger(AWG_INTEGRATION_ARM);
        waitWave();
        wait(1024);
          }
        }
    """)

    create_and_compile_awg(daq, device,  awg_program, seqr_index = 0, timeout = 2)

def sequence_multiplexed_readout(
    channels,
    frequencies,
    n_averages,
    state=None,
):
    """
    Returns an AWG sequence program (String) that specifies
    the sequence for multiplexed readout. Amplitudes and phases
    are hardcoded in the function for up to 10 channels and for
    ground and excited qubit states (simulated response of a
    readout resonator for qubit in either ground or excited state).
    Arguments:
        channels (int) -- indices of channels to create readout pulses for
        frequencies (float) -- frequencies (in Hz) of readout pulses
        n_averages )int) -- number of repetitions
    Keyword Arguments:
        state (int) -- states of measured channels to be simulated, 0 or 1
    Returns:
        (String) -- awg sequence program as string
    """

    # hard coded parameters for pulse parameters here... for ground and excited state (+ delta)
    amplitudes = np.array([0.13, 0.15, 0.16, 0.15, 0.14, 0.13, 0.17, 0.23, 0.19, 0.11])/20
    phases = np.zeros(10)
    deltas_amplitude = np.array([0.02, 0.01, -0.01, 0.02, -0.012, 0.0, 0.02, -0.012, 0.06, 0.03])
    deltas_phase = np.array([0.23, 0.31, 0.26, -0.171, 0.28, -0.31, 0.19, -0.21, 0.091, 0.29]) * np.pi/4

    n_channels = len(channels)
    assert len(frequencies) >= max(channels), "Not enough readout frequencies specified!"

    if state is None:
        state = [0] * n_channels

    for i, ch  in enumerate(channels):
        if frequencies[ch] < 0:
            frequencies[ch] = abs(frequencies[ch])
            # decide what to do here
        if state[i]:
            amplitudes[ch] = amplitudes[ch] * (1 + deltas_amplitude[ch])
            phases[ch] += deltas_phase[ch]

    # text snippet for the initialization of awg sequence
    awg_program_init = textwrap.dedent(
        """\
        const samplingRate = 1.8e9;
        // parameters for envelope
        const riseTime = 30e-9;
        const fallTime = 30e-9;
        const flatTime = 200e-9;
        const rise = riseTime * samplingRate;
        const fall = fallTime * samplingRate;
        const length = flatTime * samplingRate;
        const totalLength =  rise + length + fall;
        // define waveforms
        wave w_gauss_rise = gauss(2*rise, rise, rise/4);
        wave w_gauss_fall = gauss(2*fall, fall, fall/4);
        wave w_rise = cut(w_gauss_rise, 0, rise);
        wave w_fall = cut(w_gauss_fall, fall, 2*fall-1);
        wave w_flat = rect(length, 1.0);
        wave w_pad = zeros((totalLength-1)%16);
        // combine to total envelope
        wave readoutPulse = 1.0*join(w_rise, w_flat, w_fall, w_pad) + 0.0* w_gauss_rise;
        // init empty final waveforms
        wave w_I = zeros(totalLength);
        wave w_Q = zeros(totalLength);
    """
    )

    # text snippet for single pulse
    awg_program_singlePulse = textwrap.dedent(
        """\
        // modulate envelope for readout pulse *N*
        const f*N*_readout = _Frequency*N*_ ;
        wave w*N*_I =  _Amplitude*N*_ * readoutPulse * cosine(totalLength, 1, _Phase*N*_, f*N*_readout*totalLength/samplingRate);
        wave w*N*_Q = _Amplitude*N*_ * readoutPulse * sine(totalLength, 1, _Phase*N*_, f*N*_readout*totalLength/samplingRate);
        w_I = add(w_I, w*N*_I);
        w_Q = add(w_Q, w*N*_Q);
    """
    )

    # text snippet for main loop of .seqC
    awg_program_playWave = textwrap.dedent(
        """\
        // play waveform
        setTrigger(AWG_INTEGRATION_ARM);
        var result_averages = _nAverages_ ;
        repeat (result_averages) {
            playWave(w_I, w_Q);
            setTrigger(AWG_INTEGRATION_ARM + AWG_INTEGRATION_TRIGGER + AWG_MONITOR_TRIGGER + 1);
            setTrigger(AWG_INTEGRATION_ARM);
            waitWave();
            wait(1024);
        }
        setTrigger(0);
    """
    )

    # add all the pulses for N = ... readout channels
    # add each channel and replace indices for readout frequencies ...
    awg_program_pulses = ""
    for ch in channels:
        awg_program_pulses = awg_program_pulses + awg_program_singlePulse.replace(
            "*N*", str(ch)
        )

    # replace parameters in sequence program
    awg_program_pulses = awg_program_pulses.replace("_nChannels_", str(n_channels))
    for ch in channels:
        awg_program_pulses = awg_program_pulses.replace(
            f"_Frequency{ch}_", str(frequencies[ch])
        )
        awg_program_pulses = awg_program_pulses.replace(
            f"_Amplitude{ch}_", str(amplitudes[ch])
        )
        awg_program_pulses = awg_program_pulses.replace(
            f"_Phase{ch}_", str(phases[ch])
        )
    awg_program_playWave = awg_program_playWave.replace("_nAverages_", str(n_averages))

    return awg_program_init + awg_program_pulses + awg_program_playWave

def sequence_multiplexed_readout_non_simu(
    channels,
    frequencies,
    n_averages,
    time_origin_sec = 160e-6,
    trigger_delay_sec = 57e-9
):
    """
    Returns an AWG sequence program (String) that specifies
    the sequence for multiplexed readout. Amplitudes and phases
    are hardcoded in the function for up to 10 channels and for
    ground and excited qubit states (simulated response of a
    readout resonator for qubit in either ground or excited state).
    Arguments:
        channels (int) -- indices of channels to create readout pulses for
        frequencies (float) -- frequencies (in Hz) of readout pulses
        n_averages )int) -- number of repetitions
    Keyword Arguments:
        state (int) -- states of measured channels to be simulated, 0 or 1
    Returns:
        (String) -- awg sequence program as string
    """

    # hard coded parameters for pulse parameters here... for ground and excited state (+ delta)
    amplitudes = np.array([0.13, 0.15, 0.16, 0.15, 0.14, 0.13, 0.17, 0.23, 0.19, 0.11])/20
    phases = np.zeros(10)
    deltas_amplitude = np.array([0.02, 0.01, -0.01, 0.02, -0.012, 0.0, 0.02, -0.012, 0.06, 0.03])
    deltas_phase = np.array([0.23, 0.31, 0.26, -0.171, 0.28, -0.31, 0.19, -0.21, 0.091, 0.29]) * np.pi/4

    n_channels = len(channels)
    assert len(frequencies) >= max(channels), "Not enough readout frequencies specified!"

    if state is None:
        state = [0] * n_channels

#     for i, ch  in enumerate(channels):
# #         if frequencies[ch] < 0:
# #             frequencies[ch] = abs(frequencies[ch])
#             # decide what to do here
#         if state[i]:
#             amplitudes[ch] = amplitudes[ch] * (1 + deltas_amplitude[ch])
#             phases[ch] += deltas_phase[ch]

    # text snippet for the initialization of awg sequence

    awg_program_init = textwrap.dedent(
        """\
        const samplingRate = 1.8e9;
        const time_origin_sec = _c0_;
        cosnt trigger_delay_sec = _c1_;
        // parameters for envelope
        const riseTime = 30e-9;
        const fallTime = 30e-9;
        const flatTime = 200e-9;
        const rise = riseTime * samplingRate;
        const fall = fallTime * samplingRate;
        const length = flatTime * samplingRate;
        const totalLength =  rise + length + fall;
        // define waveforms
        wave w_gauss_rise = gauss(2*rise, rise, rise/4);
        wave w_gauss_fall = gauss(2*fall, fall, fall/4);
        wave w_rise = cut(w_gauss_rise, 0, rise);
        wave w_fall = cut(w_gauss_fall, fall, 2*fall-1);
        wave w_flat = rect(length, 1.0);
        wave w_pad = zeros((totalLength-1)%16);
        // combine to total envelope
        wave readoutPulse = 1.0*join(w_rise, w_flat, w_fall, w_pad) + 0.0* w_gauss_rise;
        // init empty final waveforms
        wave w_I = zeros(totalLength);
        wave w_Q = zeros(totalLength);
    """
    )

    # text snippet for single pulse
    awg_program_singlePulse = textwrap.dedent(
        """\
        // modulate envelope for readout pulse *N*
        const f*N*_readout = _Frequency*N*_ ;
        wave w*N*_I =  _Amplitude*N*_ * readoutPulse * cosine(totalLength, 1, _Phase*N*_, f*N*_readout*totalLength/samplingRate);
        wave w*N*_Q = _Amplitude*N*_ * readoutPulse * sine(totalLength, 1, _Phase*N*_, f*N*_readout*totalLength/samplingRate);
        w_I = add(w_I, w*N*_I);
        w_Q = add(w_Q, w*N*_Q);
    """
    )

    # text snippet for main loop of .seqC
    awg_program_playWave = textwrap.dedent(
        """\
        // play waveform
        setTrigger(AWG_INTEGRATION_ARM);
        var result_averages = _nAverages_ ;
        repeat (result_averages) {
            waitDigTrigger(1, 1);
            wait((time_origin_sec - trigger_delay_sec)*samplingRate/8);
            playWave(w_I, w_Q);
            setTrigger(AWG_INTEGRATION_ARM + AWG_INTEGRATION_TRIGGER + AWG_MONITOR_TRIGGER + 1);
            setTrigger(AWG_INTEGRATION_ARM);
            waitWave();
            wait(1024);
        }
        setTrigger(0);
    """
    )

    # add all the pulses for N = ... readout channels
    # add each channel and replace indices for readout frequencies ...
    awg_program_pulses = ""
    for ch in channels:
        awg_program_pulses = awg_program_pulses + awg_program_singlePulse.replace(
            "*N*", str(ch)
        )

    # replace parameters in sequence program
    awg_program_pulses = awg_program_pulses.replace("_nChannels_", str(n_channels))
    for ch in channels:
        awg_program_pulses = awg_program_pulses.replace(
            f"_Frequency{ch}_", str(frequencies[ch])
        )
        awg_program_pulses = awg_program_pulses.replace(
            f"_Amplitude{ch}_", str(amplitudes[ch])
        )
        awg_program_pulses = awg_program_pulses.replace(
            f"_Phase{ch}_", str(phases[ch])
        )
    awg_program_playWave = awg_program_playWave.replace("_nAverages_", str(n_averages))
    awg_program_playWave = awg_program_playWave.replace("_c0_", str(time_origin_sec))
    awg_program_playWave = awg_program_playWave.replace("_c1_", str(trigger_delay_sec))

    return awg_program_init + awg_program_pulses + awg_program_playWave

def compile_sequence(awg_module, awg_program):
    """
    Starts compilation of AWG sequence program dn loads it to device.
    Arguments:
        awg_module (awgModule) -- awgModule Object of AWG
        awg_program (String) -- specifies the awg sequence in .seqC format
    """
    awg_module.set("compiler/sourcestring", awg_program)
    while awg_module.getInt("compiler/status") == -1:
        time.sleep(1)
    assert awg_module.getInt("compiler/status") != 1, awg_module.getString(
        "/compiler/statusstring"
    )
    if awg_module.get("compiler/status") == 0:
        print("Compilation successful!")

def run_awg(daq, device):
    """
    Runs AWG sequence. Sets AWG to single shots and enables AWG.
    Arguments:
        daq (zhinst.ziDAQServer) -- Data Server Object
        device (String) -- device ID, e.g. "dev2266"
    """
    daq.asyncSetInt(f"/{device}/awgs/0/single", 1)
    daq.syncSetInt(f"/{device}/awgs/0/enable", 1)

def awg_output_amplitude(daq, device, amplitudes):
    '''
    AWG output amplitude setting

    device:          device ID
    amplitudes:      AWG output amplitude for ch1 and ch2, range from 0 to 1
    '''
    daq.setDouble('/{:s}/awgs/0/outputs/0/amplitude'.format(device), amplitudes)
    daq.setDouble('/{:s}/awgs/0/outputs/1/amplitude'.format(device), amplitudes)

def toggle_outputs(daq, device, channel=None):
    """
    Toggles signal output of UHFQA. If no channel specified toggles both.
    Arguments:
        daq (zhinst.ziDAQServer) -- Data Server Object
        device (String) -- device ID, e.g. "dev2266"
    Keyword Arguments:
        channel (int) -- list of channels to be toggled (in 0, 1)
                            (default: None)
    """

    if channel is not None:
        assert channel in [0, 1]
        channel = [channel]
    else:
        channel = [0, 1]
    for ch in channel:
        path = f"/{device}/sigouts/{ch}/on"
        if daq.getInt(path) == 0:
            daq.setInt(path, 1)
        elif daq.getInt(path) == 1:
            daq.setInt(path, 0)

def config_scope(daq, device, scope_length = 16384, scope_avg_weight = 1, scope_mode = 0,  input_sigs = [0, 0]):
    '''
    config scope

    device:             device ID
    base_rate:          sampling rate
    scope_length:       sampling points
    scope_avg_weight:   specify the averaging behaviour:
                        1:  weight, Averaging disabled.
                        >1: weight, Exponentially average the incoming scope records, updating the last scope record
                            in the history with the averaged record, see Section 3.11.5 of LabOneProgramming manual
    scope mode:     0:  Pass-through: scope segments assembled and returned unprocessed, non-interleaved.
                    1:  Moving average: entire scope recording assembled, scaling applied, averager if enabled
                        (see averager/weight), data returned in float non-interleaved format.
                    2:  Reserved for future use (average n segments).
                    3:  FFT, same as mode 1, except an FFT is applied to every segment of the scope recording
                        before averaging. See the fft/* parameters for FFT parameters.
    '''


    scope = daq.scopeModule()
    daq.setInt(f'/{device:s}/scopes/0/channel', 3)# 1: ch1; 2(DIG):ch2; 3: ch1 and ch2(DIG)
    daq.setInt(f'/{device:s}/scopes/0/channels/0/inputselect', input_sigs[0]) # 0: sigin 1; 1: sigin 2
    daq.setDouble('/%s/scopes/0/length' % device, scope_length)
    scope.set('scopeModule/averager/weight', scope_avg_weight)
    scope.set('scopeModule/mode', scope_mode)
    daq.sync()
    return scope

def init_scope(daq,device):
    scope = daq.scopeModule()
    scope.set('scopeModule/averager/weight', scope_avg_weight) # concecutive shots are averaged with an exp weight
    scope.set('scopeModule/mode', scope_mode)
    daq.setInt(f'/{device:s}/scopes/0/time', 0)
    daq.setInt(f'/{device:s}/scopes/0/trigchannel', 0)
    daq.setDouble(f'/{device:s}/scopes/0/trigdelay', 0) # in unit of second
    daq.setInt(f'/{device:s}/scopes/0/trigenable', 0)
    daq.setInt(f'/{device:s}/scopes/0/trigfalling', 0) # default is rise, if it is 1 falling is enabled
    daq.setInt(f'/{device:s}/scopes/0/trigforce', 0)
    daq.setDouble(f'/{device:s}/scopes/0/trigholdoff', 50e-3) # In units of second. Defines the time before the trigger is rearmed after a recording event
    daq.setInt(f'/{device:s}/scopes/0/trigholdoffmode', 0)
    daq.setInt(f'/{device:s}/scopes/0/trigholdoffcount', 1)
    daq.setDouble(f'/{device:s}/scopes/0/triglevel', 1e-2)
    daq.setDouble(f'/{device:s}/scopes/0/trigreference', 0.5)
    daq.setInt(f'/{device:s}/scopes/0/trigslope', 0)
    daq.getInt(f'/{device:s}/scopes/0/trigstate')
    #daq.get(f'/{device:s}/scopes/0/wave')
    daq.setInt(f'/{device:s}/scopes/0/segments/count', 1)
    daq.setInt(f'/{device:s}/scopes/0/segments/enable', 0)
    daq.setInt(f'/{device:s}/scopes/0/stream/rate', 5)# 1.8/2^n (n = 6 - 16)
    # daq.get(f'/{device:s}/scopes/0/stream/sample')
    daq.setInt(f'/{device:s}/scopes/0/triggate/enable', 0)
    daq.setInt(f'/{device:s}/scopes/0/triggate/inputselect', 0)
    daq.setInt(f'/{device:s}/scopes/0/trigholdoffmode', 0)
    daq.setDouble(f'/{device:s}/scopes/0/trighysteresis/absolute', 25e-3)
    daq.setInt(f'/{device:s}/scopes/0/trighysteresis/mode', 0) # 0, absolute; 1, relative
    daq.setDouble(f'/{device:s}/scopes/0/trighysteresis/relative', 0.15)
    daq.setInt(f'/{device:s}/scopes/0/channels/0/bwlimit', 0)
    # daq.setDouble(f'/{device:s}/scopes/0/channels/0/fullscale', 1)
    daq.setInt(f'/{device:s}/scopes/0/channels/0/inputselect', 0)
    # daq.setDouble(f'/{device:s}/scopes/0/channels/0/limitlower', -1)
    # daq.setDouble(f'/{device:s}/scopes/0/channels/0/limitupper', 1)
    daq.setDouble(f'/{device:s}/scopes/0/channels/0/offset', 0.1)
    daq.setInt(f'/{device:s}/scopes/0/stream/enables/0', 0)

def restart_avg_scope(scope):
    '''
    Set to 1 to reset the averager. The module sets averager/restart back to 0 automatically

    scope:     scope ID
    '''
    scope.set('scopeModule/averager/restart', 1)

def enable_scope(daq, device, enable = 1):
    '''
    Enable/disable scope

    daq:             device daq
    device:          device ID
    enable:    0:    disable scope
               1:    enable scope
    '''
    daq.setInt('/%s/scopes/0/enable' % device, enable)

def subscrib_scope(scope, device):
    '''
    subscrib scope

    scope:           scope ID
    device:          device ID
    '''
    scope.subscribe('/%s/scopes/*/wave' % device)

def execute_scope(scope):
    '''
    execute scope

    scope:           scope ID
    '''
    scope.execute()

def scope_progress(scope):
    '''
    scope progress

    scope:           scope ID
    '''
    return scope.progress()

def read_scope(scope):
    '''
    read scope

    scope:           scope ID
    '''
    return scope.read()

def finish_scope(scope):
    '''
    finish scope

    scope:           scope ID
    '''
    scope.finish()

def restart_and_execute_scope(scope, daq, device):
    '''
    restart and execute scope

    scope:           scope ID
    device:          device ID
    '''
    scope.set('scopeModule/averager/restart', 1)
    daq.setInt(f'/{device}/scopes/0/enable', 1)
    scope.subscribe(f'/{device}/scopes/0/wave')
    scope.execute()
    while int(scope.progress()) != 1:
        time.sleep(0.05)
        result = scope.read()
    spec = result['%s' % device]['scopes']['0']['wave'][0][0]['wave'][0]
    scope.finish()

    return spec

def do_plot1d(x_values, sweep_data, paths, readout_length, save_fig = 0):
    '''
    plot 1D data
    x_values:           x data
    sweep_data:         data from device
    paths:              data paths
    readout_length:     length of readout in second
    save_fig        0:  do not save
                    1:  do save
    '''

    y1label = 'Amplitude (V)'
    y2label = 'Phase (deg)'
    unwrap = 0

    if np.floor(x_values[0]/1e9) >= 1:
        x_values =  x_values/1e9
        xlabel = 'Frequency (GHz)'
    else:
        x_values =  x_values/1e6
        xlabel = 'Frequency (MHz)'


    base_rate = 1.8e9
    data_real = sweep_data[paths[0]]
    data_img  = sweep_data[paths[1]]
    data_cmpl = data_real+1j*data_img
    # plot
    fig, (ax1,ax2) = plt.subplots(nrows=2, sharex = True)

    ax1.set_ylabel(y1label)
#    ax1.set_xlabel(xlabel)
    ax1.semilogy(x_values, np.abs(data_cmpl), '-o', markersize = 3, c='C0')
    ax1.plot(x_values,np.abs(data_cmpl))
    ax2.set_ylabel(y2label)
    ax2.set_xlabel(xlabel)
    if unwrap:
        phase = np.unwrap(np.angle(data_cmpl))
        phase = sp.signal.detrend(phase)
    else:
        phase = np.angle(data_cmpl)
    ax2.plot(x_values, phase*180/pi, '-o', markersize = 3, c = 'C1')
    fig.align_ylabels([ax1,ax2])
#    fig.set_tight_layout(True)
    plt.show()
#     fig.canvas.manager.window.move(0, 500)
    if save_fig == 1:
        plt.savefig("spectroscopy "+ time.ctime(time.time()).replace(":","-") )

def set_integration_weights(daq, device,weights,channel,quadrature="real",demod_frequency=None):
    """
    Sets the integration weights of the UHFQA. The input signals
    are multiplied with the integrtion weights for each channel.
    Arguments:
        daq (zhinst.ziDAQServer) -- Data Server Object
        device (String) -- device ID, e.g. "dev2266"
        weights (double) -- list of double describing the integration
                              weights to be set, max. length is 4096
        channel (int) -- index of channel to set weights of
    Keyword Arguments:
        quadrature (str) -- quadrature of weights to be set,
                            either 'imag' or 'real' (default: 'real')
        demod_frequency (double) -- frequency for demodulation
                                    (default: None)
    """

    monitor_length = daq.getInt(f"/{device}/qas/0/monitor/length")
    integration_length = len(weights)
    assert integration_length <= 4096
    assert channel in range(10)
    assert quadrature in ["real", "imag"]

    # if weight is only one point, set constant weight for total length
    if len(weights) == 1:
        weights = weights * np.ones(monitor_length)

    # set lengths to the same, smallest value
    if integration_length > monitor_length:
        weights = weights[:monitor_length]
        integration_length = monitor_length
    if integration_length < monitor_length:
        monitor_length = integration_length

    # generate weights for digital demodulation
    if demod_frequency is not None:
        demod_weights = generate_demod_weights(integration_length, demod_frequency)
    else:
        demod_weights = np.ones(integration_length)

    # generate weights
    integration_weights = weights * demod_weights

    # reset
    daq.setInt(f"/{device}/qas/0/integration/length", 4096)
    daq.setVector(
        f"/{device}/qas/0/integration/weights/{channel}/{quadrature}",
        np.zeros(4096),
    )
    # set
    daq.setInt(f"/{device}/qas/0/integration/length", integration_length)
    daq.setVector(
        f"/{device}/qas/0/integration/weights/{channel}/{quadrature}",
        integration_weights,
    )

def generate_demod_weights(length, frequency, samplingRate=1.8e9, plot=False, phase=0):
    assert length <= 4096
    assert frequency > 0
    x = np.arange(0, length)
    y = np.sin(2 * np.pi * frequency * x / samplingRate + phase)
    return y

def reset_integration_weights(daq, device, channels=range(10)):
    """
    Resets the integration weights of the UHFQA to all zeros.
    If no channel specified all are reset.
    Arguments:
        daq (zhinst.ziDAQServer) -- Data Server Object
        device (String) -- device ID, e.g. "dev2266"
    Keyword Arguments:
        channels (int) --  list of indeces of channels to be reset
                             (default: range(10))
    """

    daq.setInt(f"/{device}/qas/0/integration/length", 4096)
    for ch in channels:
        daq.setVector(
            f"/{device}/qas/0/integration/weights/{ch}/real",
            np.zeros(4096),
        )
        daq.setVector(
            f"/{device}/qas/0/integration/weights/{ch}/imag",
            np.zeros(4096),
        )

def optimal_integration_weights(
    daq,
    device,
    awg,
    channel,
    frequencies,
    plot=False,
    delay=None,
):
    """
    Sets the optimal integration weights for specified channel.
    Measures IQ traces for channel being in ground/excited state
    and takes difference as optimal weights.
    Arguments:
        daq (zhinst.ziDAQServer) -- Data Server Object
        device (String) -- device ID, e.g. "dev2266"
        awg (awgModule) -- awgModule() Object of AWG
        channel (int) -- index of channel to set weights for
        frequencies (float) -- list of readout frequencies for all channels
    Keyword Arguments:
        plot (bool) -- if set, detailed plots are shown (default: False)
        delay (int) -- number of samples at the beginning of weights array
                       that are set to 0 (default: None)
    """


    daq.flush()
    monitor_length = daq.getInt(f"/{device}/qas/0/monitor/length")
    monitor_averages = daq.getInt(f"/{device}/qas/0/monitor/averages")

    reset_integration_weights(daq, device, channels=[channel])
    monitor_paths = [
        f"/{device}/qas/0/monitor/inputs/0/wave",
        f"/{device}/qas/0/monitor/inputs/1/wave",
    ]

    monitor_list = []

    ground_state = [0]
    excited_state = [1]

    for state in [ground_state, excited_state]:
        print(f"Channel {channel} in state |{','.join(str(num) for num in state)}>", flush=True)
        # set up AWG sequence program
        # readout pulse for only single channel!
        awg_program = sequence_multiplexed_readout(
            [channel],
            frequencies,
            monitor_averages,
            state=state
        )
        compile_sequence(awg, awg_program)

        # ensure all settings are synced and subscribe to monitor paths
        daq.sync()
        time.sleep(0.1)
        daq.subscribe(monitor_paths)

        # run AWG sequence and start acquisition
        run_awg(daq, device)

        # poll data
        monitor_list.append(acquisition_poll(daq, monitor_paths, monitor_length))
        print("\t\t--> Data acquired")

        # unsubscribe immediately after aquisition!
        daq.unsubscribe(monitor_paths)

    waves = []
    for polldata in monitor_list:
        for path, data in polldata.items():
            waves.append(data)

    wave_I_0 = waves[0]
    wave_Q_0 = waves[1]
    wave_I_1 = waves[2]
    wave_Q_1 = waves[3]

    weights_I = wave_I_1 - wave_I_0
    weights_Q = wave_Q_1 - wave_Q_0
    weights_I = weights_I / np.max(np.abs(weights_I))
    weights_Q = weights_Q / np.max(np.abs(weights_Q))

    if delay is not None:
        weights_I[:delay] = 0
        weights_Q[:delay] = 0

    set_integration_weights(
        daq, device, weights_I, channel, quadrature="real"
    )
    set_integration_weights(
        daq, device, weights_Q, channel, quadrature="imag"
    )

    if plot:
        # set up plot
        fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=[10, 8])

        ax1.grid("on")
        ax1.plot(wave_I_0, label="Input I", color=plt.cm.tab20(0))
        ax1.plot(wave_Q_0, label="Input Q", color=plt.cm.tab20(1))
        ax1.legend(frameon=False, loc=3)
        ax1.set_title("Qubit in ground state", position=[0.15, 0.7])
        ax1.set_xlim([0, monitor_length])

        ax2.grid("on")
        ax2.plot(wave_I_1, label="Input I", color=plt.cm.tab20(0))
        ax2.plot(wave_Q_1, label="Input Q", color=plt.cm.tab20(1))
        ax2.legend(frameon=False, loc=3)
        ax2.set_title("Qubit in excited state", position=[0.15, 0.7])
        ax2.set_xlim([0, monitor_length])

        ax3.axvline(delay, c="k", linewidth=0.5)
        ax3.grid("on")
        ax3.plot(weights_I, label="Weight I", color=plt.cm.tab20(0))
        ax3.plot(weights_Q, label="Weight Q", color=plt.cm.tab20(1))
        ax3.legend(frameon=False, loc=3)
        ax3.set_title("Integration weights for I/Q", position=[0.15, 0.7])
        ax3.set_xlim([0, monitor_length])
        ax3.set_xlabel("Samples")

        plt.show()
        plt.tight_layout()

def write_crosstalk_matrix(daq, device, matrix):
    """
    Writes the given matrix to the QA Setup crosstalk matrix of the UHFQA.
    Arguments:
        daq (zhinst.ziDAQServer) -- Connection to the Data Server
        device (String) -- device ID, e.g. "dev2266"
        matrix (2D array) -- crosstalk matrix to be written to the QA Setup tab
    """
    rows, cols = matrix.shape
    for r in range(rows):
        for c in range(cols):
            node = f"/{device}/qas/0/crosstalk/rows/{r}/cols/{c}"
            daq.setDouble(node, matrix[r, c])
    return

    '''
    plot 1D pulse experiment data

    sequence:          'Rabi', 'T1' or 'T2'
    complex_amplitude:  complex amplitude from the measurement
    x_vector:           x data
    fitting:     0:     do not fit
                 1:     do fit
    save_fig:    0:     do not save
                 1:     do save
    '''

    abs_camplitude = np.abs(complex_amplitude)
    phase_camplitude = np.angle(complex_amplitude)

    fig, (ax1,ax2) = plt.subplots(nrows=2, sharex=True)
    ax1.set_ylabel('Amplitude (V)')
    def rabi(x, x0, amp, offset):
        y = amp*(1-np.cos(pi*(x/x0)))/2+offset
        return y
    if fitting == 1:
        if sequence == "Rabi":
            best_vals, covar = scy.optimize.curve_fit(rabi, x_vector, abs_camplitude)
            ax1.plot(x_vector,rabi(x_vector, best_vals[0], best_vals[1], best_vals[2]))
            print('amplitude for pi pulse is ', best_vals[0])

        if sequence == ("T1" or "T2"):
            amp0 = abs_camplitude[0] - abs_camplitude[-1]
            init_vals = [0,1,amp0,0]
            best_vals, covar = scy.optimize.curve_fit(T1,t,y1, p0 = init_vals)
            ax1.plot(x_vector, T1(x_vector,best_vals[0],best_vals[1], best_vals[2], best_vals[3]))
            fit_para_str = sequence + ' time is ' + best_vals[1]
            print(fit_para_str)


    if sequence == "Rabi":
        label = 'Relative pulse amplitude'
#         ax1.set_xlabel(label)
        ax2.set_xlabel(label)
        ax1.set_title('Rabi oscillation with dA')
        best_vals, covar = scy.optimize.curve_fit(rabi, x_vector, abs_camplitude)

    if sequence == 'Rabi_t':
        label = 'Pulse duration (ns)'
#         ax1.set_xlabel(label)
        ax2.set_xlabel(label)
        ax1.set_title('Rabi oscillation with dt')
        x_vector = x_vector*1e9

    if sequence == "T1":
        label = 'Pulse separation ($\mu$s)'
#         ax1.set_xlabel(label)
        ax2.set_xlabel(label)
        ax1.set_title('T1')
        x_vector=x_vector*1e6
    if sequence[:2] == "T2":
        label = 'Pulse separation ($\mu$s)'
#         ax1.set_xlabel(label)
        ax2.set_xlabel(label)
        if sequence[3] == 'R':
            ax1.set_title('T2, Ramsey')
        elif sequence[3] == 'H':
            ax1.set_title('T2, Hanh echo')
        elif sequence[3] == 'C':
            ax1.set_title('T2, CPMG')
        x_vector=x_vector*1e6

    if sequence == 'Rabi_2p_t':
        label = 'Pulse duration (ns)'
#         ax1.set_xlabel(label)
        ax2.set_xlabel(label)
        ax1.set_title('Rabi $\pi/2-\pi/2$ with dt')
        x_vector = x_vector*1e9

    ax1.plot(x_vector, abs_camplitude, '-o', markersize = 3, c='C0')
    ax2.set_ylabel('Phase (deg)')
    ax2.plot(x_vector, phase_camplitude*180/pi,'-o', markersize = 3,  c='C1')
    fig.set_tight_layout(True)
    fig.align_ylabels([ax1,ax2])
    if save_fig == 1:
        plt.savefig("pulsed " + sequence + " "+ time.ctime(time.time()).replace(":","-") )