#!/usr/bin/env python
# coding: utf-8

# ## Import some packages

# In[49]:


import textwrap
import zhinst.ziPython as zp
import zhinst.utils as zu
import numpy as np
import time
from HDAWG import *
from UHFQA import *
import matplotlib
import matplotlib.pyplot as plt 


# Connect instruments

device_id_awg = 'dev8233'
device_id_qa = 'dev2528'
ip = '127.0.0.1'
port = 8006
awg, device_awg = create_api_sessions_hd(device_id_awg,use_discovery = 1, ip=ip)
qa, device_qa = create_api_sessions_uhf(device_id_qa,use_discovery =1, ip=ip)
init(awg,device_id_awg)
init(qa,device_id_qa,1,1.5)
zhinst.utils.disable_everything(awg, device_awg)
zhinst.utils.disable_everything(qa, device_qa)


# ## Generate Telegraph Noise Instance

# In[51]:


def gen_tel_noise(numPoints,tau,dt):
    
    signal = np.ones(numPoints)*(-1)**np.random.randint(0,2)
    for i in range(1,numPoints-1):
        if np.random.rand() < 1/(2*tau/dt)*np.exp(-1/(2*tau/dt)):
            signal[i+1] = - signal[i]
        else:
            signal[i+1] = signal[i]
    return signal


# ## Configure Scope


def config_qa_scope(daq, device, scope_length = 16384, scope_avg_weight = 20, scope_mode = 3, run_mode=1, streaming = False, input_sigs = [0, 0],channels=3):
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
    channels:       1:  channel 1
                    2:  channel 2
                    3:  channels 1 & 2
    
    run_mode:       1: continuous
                    2: single
                    
    '''
    
    
    scope = daq.scopeModule()
    daq.setInt(f'/{device:s}/scopes/0/channel', channels)# 1: ch1; 2(DIG):ch2; 3: ch1 and ch2(DIG)
    daq.setInt(f'/{device:s}/scopes/0/channels/0/inputselect', input_sigs[0]) # 0: sigin 1; 1: sigin 2
#     daq.setInt(f'/{device:s}/scopes/0/channels/1/inputselect', input_sigs[1])
    if run_mode == 1:
        daq.setInt('/%s/scopes/0/single' % device, 0)
    else:
        daq.setInt('/%s/scopes/0/single' % device, 1)
    
    if streaming:
        daq.setInt('/%s/scopes/0/stream/enables/*', 1)
    
        
    daq.setDouble('/%s/scopes/0/length' % device, scope_length)
    scope.set('scopeModule/averager/weight', scope_avg_weight)
    scope.set('scopeModule/mode', scope_mode) 
    daq.setDouble('/%s/sigins/*/range'%device, 1.5)
    daq.setInt('/%s/sigins/*/imp50' %device, 1)
    daq.setInt('/%s/scopes/0/trigenable' %device, 1)
    daq.setInt('/%s/scopes/0/trigchannel' %device, 1)
    daq.setDouble('/%s/scopes/0/triglevel' %device, 0.1)
    daq.setInt('/%s/scopes/0/time'%device, 0) # sets the sampling rate of the scope
    daq.setInt('/%s/scopes/0/enable' % device, 1) # arms the scope
    daq.sync()
    return scope


# ## Example

# In[87]:


"""
Zurich Instruments LabOne Python API Example

Demonstrate how to connect to a Zurich Instruments HDAWG and upload and run an
AWG program.
"""

# Copyright 2018 Zurich Instruments AG

import os
import time
import textwrap
import numpy as np
import zhinst.utils


def run_example(device_id):
    """
    Run the example: Connect to a Zurich Instruments HDAWG upload and run a
    basic AWG sequence program. It also demonstrates how to upload (replace) a
    waveform without changing the sequencer program.

    Requirements:

       HDAWG Instrument.

    Arguments:

      device_id (str): The ID of the device to run the example with. For
        example, `dev8006` or `hdawg-dev8006`.

    Returns:

      No return value.

    Raises:

      Exception: If the device is not an HDAWG.

      RuntimeError: If the device is not "discoverable" from the API.

    See the "LabOne Programming Manual" for further help, available:
      - On Windows via the Start-Menu:
        Programs -> Zurich Instruments -> Documentation
      - On Linux in the LabOne .tar.gz archive in the "Documentation"
        sub-folder.
    """


    zhinst.utils.api_server_version_check(awg)
    zhinst.utils.api_server_version_check(qa)
    
    # Create a base configuration: Disable all available outputs, awgs, demods, scopes,...
    zhinst.utils.disable_everything(awg, device_awg)
    zhinst.utils.disable_everything(qa, device_qa)
    
    # configure AWG
    # 'system/awg/channelgrouping' : Configure how many independent sequencers
    #   should run on the AWG and how the outputs are grouped by sequencer.
    #   0 : 4x2 with HDAWG8; 2x2 with HDAWG4.
    #   1 : 2x4 with HDAWG8; 1x4 with HDAWG4.
    #   2 : 1x8 with HDAWG8.
    awg.setInt(f"/{device_awg}/system/awg/channelgrouping", 1)
    set_awg(awg,device_id,1,1,1)
    
    # Some basic device configuration to output the generated waveforms.
    
    for i in range(2):
        awg.setInt(f'/{device_id_awg}/sigouts/{i}/on', 1)
        awg.setDouble(f'/{device_id_awg}/sigouts/{i}/range',2.0)
        
    # Ensure that all settings have taken effect on the device before continuing.
    awg.sync()
    
    #configure QA
    scope = config_qa_scope(qa,device_qa,scope_length=100000,scope_avg_weight = 1.0,scope_mode = 1, run_mode = 2)  

    # Number of points between the 2 pi/2 pulses
    numPoints = 16384

    # Define an AWG program as a string stored in the variable awg_program, equivalent to what would
    # be entered in the Sequence Editor window in the graphical UI.

    awg_program = textwrap.dedent(
        """
        const rate = 1.2e9;
        const length_ACprepulse = 1024;
        const length_free_evol_max = c1;
        const length_pi2 = 128;
        const length_ACpostpulse = 256;
        const dt = 512;
        const L2 = length_ACpostpulse-length_pi2;
        const amp_ac_stark = 0.4;
        const amp_telegraph = 0.2;

        wave telegraph_noise = ones(length_free_evol_max);
        wave AC_prepulse = amp_ac_stark*ones(length_ACprepulse);
        wave AC_noise = randomGauss(length_free_evol_max,1,amp_ac_stark, 0.05);
        wave AC_postpulse = amp_ac_stark*ones(length_ACpostpulse);
        
        
        wave pi2 = ones(length_pi2);

        wave w_stark_pre = join(AC_prepulse,amp_ac_stark*ones(length_pi2));
        wave w_qubit_pre = join(zeros(length_ACprepulse),pi2);

        wave w_stark_post = AC_postpulse;
        wave w_qubit_post = join(pi2,zeros(L2));
        setRate(rate);

        cvar i;
        for (i=31; i<length_free_evol_max;i=i+dt) {

        playWave(1,w_qubit_pre,2,w_stark_pre);
        playWaveIndexed(1,amp_telegraph*telegraph_noise,2,AC_noise,0,i);
        playWave(1,w_qubit_post,2,w_stark_post);
        //setTrigger(1);
        //setTrigger(0);
        wait(50);
        }
        """
    )

    awg_program = awg_program.replace("c1",str(numPoints))
    
    tic = time.perf_counter()
    # Transfer the AWG sequence program. Compilation starts automatically.
    create_and_compile_awg(awg, device_id, awg_program, 0, 5)
    toc = time.perf_counter()
    print(f"Compilation took {toc-tic:0.6f} seconds")

    # generate telegraph noise instance, convert it to HDAWG-friendly format and load it
    tic = time.perf_counter()
    whiteNoise =  np.random.normal(0.4, 0.05, numPoints)
    telegraph_noise = 0.2*gen_tel_noise(numPoints,1,1/(1.2*10**3))    
    waveform_native = zhinst.utils.convert_awg_waveform(telegraph_noise,whiteNoise)
    index = 1  # there are 3 waveforms in total (see waveform sub-tab of AWG)
    path = f"/{device_id}/awgs/0/waveform/waves/{index}"
    awg.setVector(path, waveform_native)
    toc = time.perf_counter()
    print(f"Loading new waveforms took {toc-tic:0.6f} seconds")
    
#     fig,ax = plt.subplots()
#     t = np.linspace(0,13.653,numPoints)
#     ax.plot(t,telegraph_noise)
#     ax.set(xlabel='time',ylabel='B(t)')
#     plt.show()      
    print(f"Enabling the AWG: Set /{device_awg}/awgs/0/userregs/0 to 1 to trigger waveform playback.")
    # This is the preferred method of using the AWG: Run in single mode continuous waveform playback is best achieved by
    # using an infinite loop (e.g., while (true)) in the sequencer program.
    awg.setInt(f"/{device_awg}/awgs/0/enable", 1)
          
run_example(device_id_awg)


awg.setInt(f"/{device_awg}/awgs/0/enable", 1)





