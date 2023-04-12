"""
Created on Thu Apr 14 11:35:33 2022

@author: Evangelos Vlachos <evlachos@usc.edu>
"""

import time
import textwrap
import numpy as np
import zhinst.utils as ziut
import zhinst.ziPython
from scipy.signal import gausspulse
import zhinst.toolkit as zt

device_hd_id='DEV8233'
use_discovery = 1
ip = '127.0.0.1'

def create_api_sessions_hd(device_hd_id, use_discovery = 1, ip = '10.42.0.225'):
    '''
    connect to HDAWG

    device_hd_id:               device HDAWG ID
    use_discovery              0: remote way using ip_address
                               1: used locally with USB or local Ethernet
    '''
#     global daq_hd, device_hd

    apilevel_example = 6  # The API level supported by this example.

    # Call a zhinst utility function that returns:
    # - an API session `daq` in order to communicate with devices via the data server.
    # - the device ID string that specifies the device branch in the server's node hierarchy.
    # - the device's discovery properties.

    # required_devtype = 'HDAWG'
    if not use_discovery:
        daq_hd = zhinst.ziPython.ziDAQServer(ip, 8004, apilevel_example)
        daq_hd.connectDevice(device_hd_id, '1gbe')
        device_hd = device_hd_id
        daq_hd.setInt(f'/zi/config/open', 1)
    else:
        daq_hd, device_hd, _ = ziut.create_api_session(device_hd_id, apilevel_example)

    return daq_hd, device_hd

def init(daq, device, range_hd = 1):
    '''
    Initialize device for UHFQA examples

    daq:             daq ID
    device:          device ID
    range_hd:        output range of HD
    '''
    # General setup

    # Configure the HDAWG to use one sequencer with the same waveform on all output channels.
    daq.setInt('/{}/system/awg/channelgrouping'.format(device), 1)

    # Some basic device configuration to output the generated wave.

    exp_setting = [
        ['/%s/sigouts/*/on'                % device, 1],
        ['/%s/awgs/0/outputs/0/modulation/mode'       % device, 0],
        ['/%s/awgs/0/time'                 % device, 1], # set AWG sampling rate to 1.2GHz
        ['/%s/awgs/0/userregs/0'           % device, 0],
        ['/%s/awgs/0/outputs/0/modulation/mode' % device, 0], # Output 1 modulated with Sine 1
        ['/%s/awgs/0/outputs/1/modulation/mode' % device, 0], # Output 2 modulated with Sine 2
        ['/%s/sines/0/phaseshift'   % device, 0],   # Sine 1 phase
        ['/%s/sines/1/phaseshift'  % device, 90],   # Sine 2 phase
        ['/%s/sines/0/oscselect'    % device, 0],   # Select oscillator 1 for Sine 1
        ['/%s/sines/1/oscselect'    % device, 0],   # Select oscillator 1 for Sine 2
        ['/%s/sines/0/amplitudes/*' % device, 0],   # Turn off CW signals
        ['/%s/sines/1/amplitudes/*' % device, 0],   # Turn off CW signals
        ['/%s/sines/0/enables/*'    % device, 0],   # Turn off CW signals
        ['/%s/sines/1/enables/*'    % device, 0],    # Turn off CW signals
        ['/%s/oscs/*/freq'  % device, 0] #set modulation frequency to 0
    ]
    daq.set(exp_setting)
    time.sleep(0.1)
    daq.setInt('/%s/system/clocks/referenceclock/source' % device, 0)
    daq.setDouble('/%s/system/clocks/sampleclock/freq'   % device, 2.4e+09)


def init_wfms(awg,device_awg,nPoints,nWfs):
    awg_program = ""
    for i in range(nWfs):
        awg_program = awg_program + "wave w%d = 0.01*%d*ones(%d);\nplayWave(w%d);\n" %(i,i,2*nPoints,i)
    print('Initializing waveforms')
    create_and_compile_awg(awg, device_awg, awg_program, seqr_index = 0, timeout = 1)
    enable_awg(awg,device_awg,enable=1)


def awg_seq(awg, fs=1.2e9, amp_q = 1,nSteps=100, nPoints=1024,pi2Width=100,nPointsPre=900,nPointsPost=120,
            measPeriod=400e-6,sequence='rabi',qubit_drive_dur=20e-6,mu=0,sigma=0,B0=0,
            Tmax=2e-6,nAverages=128,active_reset=False,axis='X'):
    """

    Creates and uploads the sequence file to be used in the measurement.

    Args:
        awg (class): The API class instance of the AWG.
        fs (float, optional): Sampling rate of the AWG. Defaults to 1.2e9.
        amp_q (float, optional): amplitude of pi2 pulses (ramsey) or qubit drive (rabi). Defaults to 1.
        nSteps (int, optional): Number of points in the sequence. Defaults to 100.
        pi2Width (int, optional): Pi/2 pulse length in units of dt=1/fs. Defaults to 100.
        nPointsPre (int, optional): Number of points in the AC pre-pulse in units of dt. Defaults to 900.
        nPointsPost (int, optional): Number of points in the AC post-pulse in units of dt. Defaults to 120.
        measPeriod (float, optional): Delay between experiments in units of seconds. Defaults to 400e-6.
        sequence (string, optional): The type of experiment. Defaults to 'rabi'.
        qubit_drive_dur (int, optional): duration of ON pulse in Rabi and spectroscopy measurements. Defaults to 30e-6.
        mu (float, optional): amplitude of AC stark tone. Defaults to 0.
        sigma (float, optional): Amplitude of white noise waveform. Defaults to 0.
        B0 (float, optional): Amplitude of Generalized Markovian Noise waveform. Defaults to 0.
        Tmax (float, optional): Maximum experiment duration. Defaults to 2e-6.
        nAverages (int, optional): DESCRIPTION. Defaults to 128.
        active_reset (boolean, optional): DESCRIPTION. Defaults to False.

    Returns:
        None.

    """

    active_reset_program = ('''
        waitDigTrigger(1);
        wait(1);
        wait(2000);
        //playZero(256,AWG_RATE_37P5MHZ);
        if (getDigTrigger(2) == 0) {
            //playZero(32);
            wait(10);
        } else {
            _apply_reset_
            }
        //playZero(302,AWG_RATE_9P4MHZ);
        wait(10000);
      ''')

    if sequence == 'qubit spec':
        awg_program = textwrap.dedent("""
        const f_s = _c0_;
        const f_c = 2.4e9;      // clock rate
        const f_seq = f_c/8;     // sequencer instruction rate
        const measInt_fs = 1.17e6; // sampling rate during passive reset period
        const dt = 1/f_seq;
        const trig_interval = _c1_; // one cycle
        const period_wait_sample = floor(_c1_*measInt_fs);
        const wave_dur_sample  = _c2_;
        wave w = _c3_*ones(wave_dur_sample);

        wave w_marker = 2*marker(256,1);

       _add_AC_stark_

        // Beginning of the core sequencer program executed on the HDAWG at run time


        repeat(_c4_) {
            // OFF Measurement
            playZero(wave_dur_sample,AWG_RATE_600MHZ);
            playWave(1,w_marker);
            playZero(period_wait_sample,AWG_RATE_1P2MHZ);
            // ON measurement
            playWave(1,w);
            playWave(1,w_marker);
            playZero(period_wait_sample,AWG_RATE_1P2MHZ);
                    }
            """)

        if mu != 0:
            qubit_drive_tone = amp_q*np.ones(qubit_drive_dur)
            qubit_drive_tone = qubit_drive_tone[...,None]
            AC_stark_tone = mu*np.ones(qubit_drive_dur)
            AC_stark_tone = AC_stark_tone[...,None] # turns row vector into column vector
            wfm_2D_arr = np.hstack((qubit_drive_tone,AC_stark_tone))
            fileName = "qubit_spec_wfm"
            np.savetxt("C:/Users/LFL/Documents/Zurich Instruments/LabOne/WebServer/awg/waves/"+fileName+".csv", wfm_2D_arr, delimiter = ",")
            txt = "wave wfm = \"%s\";\n"%(fileName)
            txt += "wave ACprepulse = %f*ones(%d);\n"%(mu,nPointsPre)
            txt += "wave ACpostpulse = %f*ones(%d);\n"%(mu,nPointsPost)
            awg_program = awg_program.replace('_add_white_noise_',txt)
            txt_loop = "playWave(2,ACprepulse);\n"
            txt_loop += "playWave(wfm);\n"
            txt_loop += "playWave(2,ACpostpulse);\n"

            awg_program = awg_program.replace('_add_AC_stark_',txt)
            awg_program = awg_program.replace('playZero(wave_dur_sample,AWG_RATE_600MHZ)','playWave(2,%f*ones(%d))'%(mu,nPointsPre+qubit_drive_dur+nPointsPost))
            awg_program = awg_program.replace('playWave(1,w);',txt_loop)

        else:
           awg_program = awg_program.replace('_add_AC_stark_','')

        awg_program = awg_program.replace('_c0_', str(fs))
        awg_program = awg_program.replace('_c1_', str(measPeriod))
        awg_program = awg_program.replace('_c2_', str(qubit_drive_dur))
        awg_program = awg_program.replace('_c3_', str(amp_q))
        awg_program = awg_program.replace('_c4_', str(nAverages))


    elif sequence =='rabi':
        awg_program = textwrap.dedent(
        """
        const f_s = _c0_;
        const f_c = 2.4e9;      // clock rate
        const measInt_fs = 1.17e6; // sampling rate during passive reset period
        const f_seq = f_c/8;     // sequencer instruction rate
        const dt = 1/f_seq;
        const trigger_interval= _c1_; // one meas cycle in sec
        const period_wait_sample = floor(_c1_*measInt_fs);
        var i=0;

        wave w_marker = 2*marker(256,1);

        _add_white_noise_
        // Beginning of the core sequencer program executed on the HDAWG at run time
        repeat(_c2_){
            for (i=0; i<_c3_; i++) {
                    _add_AC_pre_pulse_
                    executeTableEntry(i);
                    _add_AC_post_pulse_
                    playWave(1,w_marker);
                    playZero(period_wait_sample,AWG_RATE_1P2MHZ);
          }
        }
        """)


        if mu != 0 or B0 != 0:
            qubit_drive_tone = amp_q*np.ones(qubit_drive_dur)
            qubit_drive_tone = qubit_drive_tone[...,None]
            AC_stark_tone = mu*np.ones(qubit_drive_dur)
            AC_stark_tone = AC_stark_tone[...,None] # turns row vector into column vector
            wfm_2D_arr = np.hstack((qubit_drive_tone,AC_stark_tone))
            fileName = "rabi_wfm"
            np.savetxt("C:/Users/LFL/Documents/Zurich Instruments/LabOne/WebServer/awg/waves/"+fileName+".csv", wfm_2D_arr, delimiter = ",")
            txt = "//Make waveforms\nwave wfms = \"%s\";\n"%(fileName)+"assignWaveIndex(wfms,0);\n"
            txt += "wave ACprepulse = %f*ones(%d);\nwave qubit_channel_pre_pulse = zeros(%d);\nassignWaveIndex(qubit_channel_pre_pulse,ACprepulse,1);\n\n"%(mu,nPointsPre,nPointsPre)
            awg_program = awg_program.replace('_add_white_noise_',txt)
            txt2 = 'executeTableEntry(_c5_);'
            txt3 = txt2
            awg_program = awg_program.replace('_add_AC_pre_pulse_',txt2)
            awg_program = awg_program.replace('_add_AC_post_pulse_',txt3)
        else:
            if axis == 'X':
                awg_program = awg_program.replace('_add_white_noise_',"wave drive_pulse=_c4_*ones(_c5_);\n"+"assignWaveIndex(drive_pulse,0);\n")
            elif axis == 'Y':
                awg_program = awg_program.replace('_add_white_noise_',"wave drive_pulse=_c4_*ones(_c5_);\n"+"assignWaveIndex(zeros(N),drive_pulse,0);\n")
            awg_program = awg_program.replace('_add_AC_pre_pulse_','')
            awg_program = awg_program.replace('_add_AC_post_pulse_','')

        awg_program = awg_program.replace('_c0_', str(fs))
        awg_program = awg_program.replace('_c1_', str(measPeriod))
        awg_program = awg_program.replace('_c2_',str(nAverages))
        awg_program = awg_program.replace('_c3_',str(nSteps))
        awg_program = awg_program.replace('_c4_',str(amp_q))
        awg_program = awg_program.replace('_c5_',str(nPoints))
        awg.setInt('/dev8233/triggers/out/0/source',4)

    elif sequence =='ramsey':

        awg_program = textwrap.dedent("""
        // Define experimental variables
        const f_s = _c0_;
        const f_c = 2.4e9;      // clock rate
        const measInt_fs = 1.17e6; // sampling rate during passive reset period
        const f_seq = f_c/8;     // sequencer instruction rate
        const dt = 1/f_seq;         // one clock cycle in sec
        const trigger_interval= _c1_; // one meas cycle in sec
        const free_evol_dur  = _c2_;    // max free evolution time in sec
        const period_wait_sample = floor(_c1_*measInt_fs);
        const N  = floor(_c2_*f_s);
        var i;

        wave pi2pulse = _c3_*ones(_c4_);
        wave w_marker = marker(256,1);

        _add_white_noise_
        _active_reset_pulses_
        // Beginning of the core sequencer program executed on the HDAWG at run time
          repeat(_c5_){
            for (i=0; i<_c6_; i++) {
                     _add_AC_pre_pulse_
                     _add_AC_post_pulse_
                    executeTableEntry(i);
                    playWave(1,w_marker);
                    playZero(period_wait_sample,AWG_RATE_1P2MHZ);
                    _active_reset_
        }
        }
        """)

        if mu != 0 or sigma != 0 or B0 != 0:
            fileName = "ramsey_wfm"
            txt = "//Load custom waveform\nwave wfms = \"%s\";\n"%(fileName)+"assignWaveIndex(wfms,0);\n"

            if mu != 0:
                txt += "//Make pre-pulse\nwave AC_stark_tone_pre = %f*ones(%d);\n"%(mu,nPointsPre)+"wave pi2pulse_pre_zeros = zeros(%d);\n"%(nPointsPre-pi2Width)+"wave pi2pulse_pre=join(pi2pulse_pre_zeros,pi2pulse);\nassignWaveIndex(pi2pulse_pre,AC_stark_tone_pre,1);\n\n"
                txt += "//Make post-pulse\nwave AC_stark_tone_post = %f*ones(%d);\n"%(mu,nPointsPost)+"wave pi2pulse_post_zeros = zeros(%d);\n"%(nPointsPost-pi2Width)+ "wave pi2pulse_post=join(pi2pulse,pi2pulse_post_zeros);\nassignWaveIndex(pi2pulse_post,AC_stark_tone_post,2);\n\n"
                txt2 = 'executeTableEntry(_c6_);'
                txt3 = 'executeTableEntry(_c6_+1);'
                awg_program = awg_program.replace('_add_AC_pre_pulse_',txt2)
                awg_program = awg_program.replace('_add_AC_post_pulse_',txt3)

                if active_reset == True:
                    awg_program = awg_program.replace('_active_reset_pulses_','//Make reset pulses\nwave pipulse = join(pi2pulse_pre,pi2pulse_post);\nwave ac_pulse = join(AC_stark_tone_pre,AC_stark_tone_post);\nassignWaveIndex(pipulse,ac_pulse,3);\n')
                    active_reset_program = active_reset_program.replace('_apply_reset_','executeTableEntry(_c6_+2);')
                    active_reset_program = active_reset_program.replace('_c6_',str(nSteps))
                    active_reset_program = active_reset_program.replace('_do_nothing_','')
                    awg_program = awg_program.replace('_active_reset_',active_reset_program)
                    awg_program = awg_program.replace('playZero(period_wait_sample,AWG_RATE_1P2MHZ);','')
                elif active_reset == False:
                    awg_program = awg_program.replace('_active_reset_pulses_','')
                    awg_program = awg_program.replace('_active_reset_','')

            elif mu == 0:
                awg_program = awg_program.replace('_add_AC_pre_pulse_','executeTableEntry(_c6_);')
                awg_program = awg_program.replace('_add_AC_post_pulse_','executeTableEntry(_c6_);')
                if active_reset == True:
                    awg_program = awg_program.replace('_active_reset_pulses_','//Make reset pulses\nwave pipulse = join(pi2pulse,pi2pulse);\nassignWaveIndex(pi2pulse,1);\nassignWaveIndex(pipulse,2);\n')
                    active_reset_program = active_reset_program.replace('_apply_reset_','executeTableEntry(_c6_+1);')
                    active_reset_program = active_reset_program.replace('_c6_',str(nSteps))
                    active_reset_program = active_reset_program.replace('_do_nothing_','')
                    awg_program = awg_program.replace('_active_reset_',active_reset_program)
                    awg_program = awg_program.replace('playZero(period_wait_sample,AWG_RATE_1P2MHZ);','')
                elif active_reset == False:
                    awg_program = awg_program.replace('_active_reset_pulses_','assignWaveIndex(pi2pulse,1);\n')
                    awg_program = awg_program.replace('_active_reset_','')
            awg_program = awg_program.replace('_add_white_noise_',txt)

        elif B0 == 0 and mu == 0:
            awg_program = awg_program.replace('_add_white_noise_',"")
            awg_program = awg_program.replace('executeTableEntry(i);',"playZero(i*_c7_,AWG_RATE_%dMHZ);"%(int(fs/1e6)))
            awg_program = awg_program.replace('_add_AC_pre_pulse_','playWave(pi2pulse);')
            awg_program = awg_program.replace('_add_AC_post_pulse_','playWave(pi2pulse);')

            if active_reset == True:
                awg_program = awg_program.replace('_active_reset_pulses_','wave pipulse = join(pi2pulse,pi2pulse);\n')
                active_reset_program = active_reset_program.replace('_apply_reset_','playWave(1,pipulse);\n')
                active_reset_program = active_reset_program.replace('_do_nothing_','')
                awg_program = awg_program.replace('_active_reset_',active_reset_program)
                awg_program = awg_program.replace('playZero(period_wait_sample,AWG_RATE_1P2MHZ);','')
            else:
                awg_program = awg_program.replace('_active_reset_pulses_','')
                awg_program = awg_program.replace('_active_reset_','')


        awg_program = awg_program.replace('_c0_', str(fs))
        awg_program = awg_program.replace('_c1_', str(measPeriod))
        awg_program = awg_program.replace('_c2_', str(Tmax))
        awg_program = awg_program.replace('_c3_', str(amp_q))
        awg_program = awg_program.replace('_c4_',str(pi2Width))
        awg_program = awg_program.replace('_c5_',str(nAverages))
        awg_program = awg_program.replace('_c6_',str(nSteps))
        awg_program = awg_program.replace('_c7_',str(round(Tmax/nSteps*fs)))


        awg.setInt('/dev8233/triggers/out/0/source',4)

    elif sequence == 'echo':

        awg_program = textwrap.dedent("""
        // Define experimental variables
        const f_s = _c0_;
        const f_c = 2.4e9;      // clock rate
        const f_seq = f_c/8;     // sequencer instruction rate
        const measInt_fs = 1.17e6; // sampling rate during passive reset period
        const dt = 1/f_seq;        // one clock cycle in sec
        const trigger_interval= _c1_; // one meas cycle in sec
        const free_evol_dur  = _c2_;    // max free evolution time in sec
        const period_wait_sample = floor(_c1_*measInt_fs);
        const N  = floor(_c2_*f_s);
        var i;

        wave w_marker = marker(256,1);
        wave pi2pulse = _c3_*ones(_c4_);
        wave pipulse = 0.245*ones(2*_c4_);

        _add_white_noise_

        // Beginning of the core sequencer program executed on the HDAWG at run time
        repeat(_c5_){
            for (i=0; i<_c6_; i++) {
                _add_AC_pre_pulse_
                executeTableEntry(i);_add_mid_pulse_executeTableEntry(i);
                _add_AC_post_pulse_
                playWave(1,w_marker);
                playZero(period_wait_sample,AWG_RATE_1P2MHZ);
                _active_reset_
                }

         }
         """)

        fileName = "echo_wfm"
        if mu != 0 and B0 == 0:
            txt = "//Make pre-pulse\nwave AC_stark_tone_pre = %f*ones(%d);\n"%(mu,nPointsPre)+"wave pi2pulse_pre_zeros = zeros(%d);\n"%(nPointsPre-pi2Width)+"wave pi2pulse_pre=join(pi2pulse_pre_zeros,pi2pulse);\nassignWaveIndex(pi2pulse_pre,AC_stark_tone_pre,1);\n\n"
            txt += "//Make post-pulse\nwave AC_stark_tone_post = %f*ones(%d);\n"%(mu,nPointsPost)+"wave pi2pulse_post_zeros = zeros(%d);\n"%(nPointsPost-pi2Width)+ "wave pi2pulse_post=join(pi2pulse,pi2pulse_post_zeros);\nassignWaveIndex(pi2pulse_post,AC_stark_tone_post,2);\n\n"
            txt += "//Make mid-pulse\nwave AC_stark_mid_pulse = %f*ones(2*_c4_);\n"%(mu)+"assignWaveIndex(pipulse,AC_stark_mid_pulse,3);\n"
            txt += "\nwave wfms = \"%s\";\n"%(fileName)+"assignWaveIndex(wfms,0);\n"
            awg_program = awg_program.replace('_add_white_noise_',txt)
            txt1 = "playWave(join(pi2pulse_pre,pipulse,pi2pulse_post),join(AC_stark_tone_pre,AC_stark_mid_pulse,AC_stark_tone_post));"
            txt2 = 'executeTableEntry(_c6_);'
            txt3 = 'executeTableEntry(_c6_+1);'
            # txt4 = "playWave(pipulse,%f*ones(2*_c4_));"%(mu)
            txt4 = 'executeTableEntry(_c6_+2);'
            awg_program = awg_program.replace('_add_first_point_',txt1)
            awg_program = awg_program.replace('_add_AC_pre_pulse_',txt2)
            awg_program = awg_program.replace('_add_AC_post_pulse_',txt3)
            awg_program = awg_program.replace('_add_mid_pulse_',txt4)
        elif B0 == 0 and mu == 0:
            awg_program = awg_program.replace('_add_white_noise_',"")
            awg_program = awg_program.replace('executeTableEntry(i);',"playZero(i*_c7_,AWG_RATE_%dMHZ);"%(int(fs/1e6)))
            awg_program = awg_program.replace('_add_AC_pre_pulse_','playWave(pi2pulse);')
            awg_program = awg_program.replace('_add_mid_pulse_','playWave(pipulse);')
            awg_program = awg_program.replace('_add_AC_post_pulse_','playWave(pi2pulse);')
        elif B0 != 0:
            txt = "wave wfms = \"%s\";\n"%(fileName)+"assignWaveIndex(wfms,0);\n"
            fileName = "echo_wfm"
            awg_program = awg_program.replace('_add_white_noise_',txt)
            awg_program = awg_program.replace('_add_AC_pre_pulse_',"playWave(1,pi2pulse);")
            awg_program = awg_program.replace('_add_AC_post_pulse_',"playWave(1,pi2pulse);")
            awg_program = awg_program.replace('_add_mid_pulse_','playWave(2,pipulse);')

        if active_reset == True:
            active_reset_program = active_reset_program.replace('_apply_reset_','playWave(1,pipulse);\n')
            active_reset_program = active_reset_program.replace('_do_nothing_','')
            awg_program = awg_program.replace('_active_reset_',active_reset_program)
            awg_program = awg_program.replace('playZero(period_wait_sample,AWG_RATE_1P2MHZ);','')
        else:
            awg_program = awg_program.replace('_active_reset_','')

        awg_program = awg_program.replace('_c0_', str(fs))
        awg_program = awg_program.replace('_c1_', str(measPeriod))
        awg_program = awg_program.replace('_c2_', str(Tmax))
        awg_program = awg_program.replace('_c3_', str(amp_q))
        awg_program = awg_program.replace('_c4_',str(pi2Width))
        awg_program = awg_program.replace('_c5_',str(nAverages))
        awg_program = awg_program.replace('_c6_',str(nSteps))
        awg_program = awg_program.replace('_c7_',str(round(Tmax/nSteps*fs)))
        awg.setInt('/dev8233/triggers/out/0/source',4)

    elif sequence == 'echo_v2':

        awg_program = textwrap.dedent("""
        // Define experimental variables
        const f_c = 2.4e9;      // clock rate
        const f_seq = f_c/8;     // sequencer instruction rate
        const dt = 1/f_seq;        // one clock cycle in sec
        const measInt_fs = 1.17e6; // sampling rate during passive reset period
        const trigger_interval= _c1_; // one meas cycle in sec
        const period_wait_sample = floor(_c1_/dt);
        var i;

        wave w_marker = marker(256,1);


        // Beginning of the core sequencer program executed on the HDAWG at run time
        repeat(_c5_){
            _main_body_
         }
         """)

        txt = ''
        for i in range(nSteps):
            txt += 'playWave("echo_wfm_%03d");\nplayWave(1,w_marker);\nwait(period_wait_sample);\n'%i

        awg_program = awg_program.replace('_main_body_',txt)
        awg_program = awg_program.replace('_c1_', str(measPeriod))
        awg_program = awg_program.replace('_c5_',str(nAverages))

    elif sequence =='T1':
        awg_program = textwrap.dedent("""
        const f_s = _c0_;
        const f_c = 2.4e9;      // clock rate
        const f_seq = f_c/8;     // sequencer instruction rate
        const measInt_fs = 1.17e6; // sampling rate during passive reset period
        const dt = 1/f_seq;
        const trigger_interval= _c1_; // one meas cycle in sec
        const tmax  = _c2_;    // max waiting time
        const period_wait_sample = floor(_c1_*measInt_fs);
        const N  = floor(_c2_*f_s);
        var i;

        wave w_marker = 2*marker(256,1);
        wave pipulse = _c3_*ones(_c4_);
        _add_white_noise_
        // Beginning of the core sequencer program executed on the HDAWG at run time

        repeat(_c5_){
            for (i=0; i<_c6_; i++) {
                    _add_AC_pre_pulse_
                    executeTableEntry(i);
                    _add_AC_post_pulse_
                    playWave(1,w_marker);
                    playZero(period_wait_sample,AWG_RATE_1P2MHZ);
                        }

          }

        """)


        if mu != 0:
            fileName = "T1_wfm"
            txt = "wave pipulse_pre_zeros = zeros(%d);\n"%(nPointsPre-2*pi2Width)+"wave pipulse_pre=join(pipulse_pre_zeros,pipulse);"
            txt += "wave wfms = \"%s\";\n"%(fileName)+"assignWaveIndex(wfms,0);\n"
            txt += "wave AC_stark_tone_pre = %f*ones(%d);\n\nassignWaveIndex(pipulse_pre,AC_stark_tone_pre,1);\n\n"%(mu,nPointsPre)
            txt += "wave AC_stark_tone_post = %f*ones(%d);\n\n\nassignWaveIndex(zeros(%d),AC_stark_tone_post,2);\n\n"%(mu,nPointsPost,nPointsPost)
            awg_program = awg_program.replace('_add_white_noise_',txt)
            txt1 = "executeTableEntry(_c6_);"
            txt2 = "executeTableEntry(_c6_+1);"
            awg_program = awg_program.replace('_add_AC_pre_pulse_',txt1)
            awg_program = awg_program.replace('_add_AC_post_pulse_',txt2)
        else:
            awg_program = awg_program.replace('_add_white_noise_',"")
            awg_program = awg_program.replace('executeTableEntry(i);',"playZero(i*_c7_,AWG_RATE_1200MHZ);")
            awg_program = awg_program.replace('_add_AC_pre_pulse_','playWave(pipulse);')
            awg_program = awg_program.replace('_add_AC_post_pulse_','')

        awg_program = awg_program.replace('_c0_', str(fs))
        awg_program = awg_program.replace('_c1_', str(measPeriod))
        awg_program = awg_program.replace('_c2_', str(Tmax))
        awg_program = awg_program.replace('_c3_', str(amp_q))
        awg_program = awg_program.replace('_c4_',str(2*pi2Width))
        awg_program = awg_program.replace('_c5_',str(nAverages))
        awg_program = awg_program.replace('_c6_',str(nSteps))
        awg_program = awg_program.replace('_c7_',str(round(Tmax/nSteps*fs)))
        awg.setInt('/dev8233/triggers/out/0/source',4)

    elif sequence == 'single_shot':

        awg_program = textwrap.dedent("""
        const f_s = _c0_;
        const f_c = 2.4e9;      // clock rate
        const f_seq = f_c/8;     // sequencer instruction rate
        const measInt_fs = 1.17e6; // sampling rate during passive reset period
        const dt = 1/f_seq;
        const trigger_interval= _c1_; // one meas cycle in sec
        const tmax  = _c2_;    // max waiting time
        const period_wait_sample = floor(_c1_*measInt_fs);
        const N  = floor(_c2_*f_s);
        var i=0;

        wave w_marker = 2*marker(512,1);
        wave pi2pulse = _c3_*ones(_c4_);
        wave pipulse = _c3_*ones(2*_c4_);
        _add_white_noise_
        // Beginning of the core sequencer program executed on the HDAWG at run time

        repeat(_c5_) {
            // OFF Measurement
            playZero(N,AWG_RATE_600MHZ);
            playWave(1,w_marker);
            playZero(period_wait_sample,AWG_RATE_1P2MHZ);
            //waitDigTrigger(1);
            //wait(1);
            //playZero(48,AWG_RATE_37P5MHZ);
            //if (getDigTrigger(2) == 0) {
              //      playZero(32);
            //} else {
              //  playWave(1,pi_pulse,2,AC_tone);
                //}
            //playZero(32,AWG_RATE_2P34MHZ);
            // ON measurement
            playWave(1,pipulse);
            playWave(1,w_marker);
            playZero(period_wait_sample,AWG_RATE_1P2MHZ);
            //waitDigTrigger(1);
            //wait(1);
            //playZero(48,AWG_RATE_37P5MHZ);
            //if (getDigTrigger(2) == 0) {
              //      playZero(32);
            //} else {
              //  playWave(1,pi_pulse,2,AC_tone);
                //}
            //playZero(32,AWG_RATE_2P34MHZ);
}

        """)



        if mu != 0 or sigma != 0:
            txt = "wave AC = _c6_*ones(_c4_);"
            txt = "//Make pre-pulse\nwave AC_stark_tone_pre = %f*ones(%d);\n"%(mu,nPointsPre)+"wave pi2pulse_pre_zeros = zeros(%d);\n"%(nPointsPre-pi2Width)+"wave pi2pulse_pre=join(pi2pulse_pre_zeros,pi2pulse);\n"
            txt += "//Make post-pulse\nwave AC_stark_tone_post = %f*ones(%d);\n"%(mu,nPointsPost)+"wave pi2pulse_post_zeros = zeros(%d);\n"%(nPointsPost-pi2Width)+ "wave pi2pulse_post=join(pi2pulse,pi2pulse_post_zeros);\nwave pi_pulse = join(pi2pulse_pre,pi2pulse_post);\nwave AC_tone = join(AC_stark_tone_pre,AC_stark_tone_pre);\n"
            txt1 = "playWave(2,_c6_*ones(N));"
            txt2 = "playWave(1,pi_pulse,2,AC_tone);"
            awg_program = awg_program.replace('_add_white_noise_',txt)
            awg_program = awg_program.replace('playZero(N,AWG_RATE_600MHZ);',txt1)
            awg_program = awg_program.replace('playWave(1,pipulse);',txt2)
        else:
            awg_program = awg_program.replace('_add_white_noise_',"")


        awg_program = awg_program.replace('_c0_', str(fs))
        awg_program = awg_program.replace('_c1_', str(measPeriod))
        awg_program = awg_program.replace('_c2_', str(1e-6))
        awg_program = awg_program.replace('_c3_', str(amp_q))
        awg_program = awg_program.replace('_c4_',str(int(pi2Width)))
        awg_program = awg_program.replace('_c5_',str(nAverages))
        awg_program = awg_program.replace('_c6_',str(mu))
        awg.setInt('/dev8233/triggers/out/0/source',4)

    elif sequence == 'mixer-calib':
        awg_program = textwrap.dedent("""

            const N = 1024;
            wave w_const = _c0_*ones(N);
            wave w_zeros = zeros(N);

            while (true) {
                playWave(1,2,w_const,1,2,w_zeros);
                waitWave();
                }

                                      """)

        awg_program = awg_program.replace('_c0_', str(amp_q))

    create_and_compile_awg(awg, awg_program, seqr_index = 0, timeout = 10)

def create_and_compile_awg(awg, awg_program, seqr_index= 0, timeout=1,verbose=0):
    """
    Compiles and uploads the sequence file into the AWG

    Args:
        awg (class): AWG class instance.
        awg_program (string): Sequence file.
        seqr_index (int, optional): Which AWG to upload the sequence file to. Defaults to 0.
        timeout (float, optional): How long to wait for sequence file upload before time out. Defaults to 1.
        verbose (TYPE, optional): DESCRIPTION. Defaults to 0.

    Raises:
        Exception: DESCRIPTION.

    Returns:
        None.

    """

    awgModule = awg.awgModule()
    awgModule.set('device', 'dev8233')
    awgModule.set('index', seqr_index)
    awgModule.execute()
    """Compile and upload awg_program as .elf file"""
    if verbose==0:
        # print("Starting compilation.")
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
    if elfStatus == 0 and verbose == 1:
        print("Upload to the instrument successful.")
    if elfStatus == 1:
        raise Exception("Upload to the instrument failed.")


def set_triggers_out(daq, device, trigger_ch = 0, source = 0, trigger_delay = -268.2e-15):
    '''
    set trigger outputs

    daq:             daq ID
    device:          device ID
    trigger_ch:      physical channel from 0 to 7/4
    source:          trigger source
                     Allowed Values:
                     0 Trigger output is assigned to AWG Trigger 1, controlled by AWG sequencer commands.
                     1 Trigger output is assigned to AWG Trigger 2, controlled by AWG sequencer commands.
                     2 Trigger output is assigned to AWG Trigger 3, controlled by AWG sequencer commands.
                     3 Trigger output is assigned to AWG Trigger 4, controlled by AWG sequencer commands.
                     4 Output is assigned to Output 1 Marker 1.
                     5 Output is assigned to Output 1 Marker 2.
                     6 Output is assigned to Output 2 Marker 1.
                     7 Output is assigned to Output 2 Marker 2.
                     8 Output is assigned to Trigger Input 1.
                     9 Output is assigned to Trigger Input 2.
                    10 Output is assigned to Trigger Input 3.
                    11 Output is assigned to Trigger Input 4.
                    12 Output is assigned to Trigger Input 5.
                    13 Output is assigned to Trigger Input 6.
                    14 Output is assigned to Trigger Input 7.
                    15 Output is assigned to Trigger Input 8.
                    17 Output is set to high.
                    18 Output is set to low.
    delay:           Trigger delay, controls the fine delay of the trigger output. The resolution is 78 ps.

    '''

    daq.setInt(f'/{device:s}/triggers/out/{trigger_ch:d}/source', source)
    daq.setDouble(f'/{device:s}/triggers/out/{trigger_ch:d}/delay', trigger_delay)

def set_sine_amp(daq, device, sine_ch = [0, 1], sine_amp=[0.1, 0.1]):
    for i in sine_ch:
        daq.setDouble(f'/{device}/sines/{i}/amplitudes/{i}', sine_amp[i])

def set_output_amp_scale(daq, device, output_ch=[0, 1], output_amp = [1, 1]):
    for i in output_ch:
        daq.setDouble(f'/{device}/awgs/0/outputs/{i}/gains/{i}', output_amp[i])

def set_phase(daq, device, diff_phase = 0, phase_I = 0, channel_I = 0, channel_Q = 1):

    '''
    set sine phase of I and Q channel (phase of sine generator ch1 and ch2)

    daq:                daq ID
    device:             device ID
    diff_phase:         phase_Q - phase_I, range from 0 to 360 deg
    phase_I:            phase of I signal, range from 0 to 360 deg
    channel_I:          output channel for I signal, int from 0 to 7/3
    channel_Q:          output channel for Q signal, int from 0 to 7/3
    '''

    set_phase_I_str='/{:s}/sines/%d/phaseshift' % int(channel_I)
    set_phase_Q_str='/{:s}/sines/%d/phaseshift' % int(channel_Q)
    daq.setDouble(set_phase_I_str.format(device), float(phase_I))

    daq.setDouble(set_phase_Q_str.format(device), float(phase_I + diff_phase))

def set_off_I(daq, device, off_I = 0, channel_I = 0):
    '''
    set output offset for I signal

    daq:                daq ID
    device:             device ID
    off_I:              offset of waveform ouTput for I signal
    channel_I:          output channel for I signal
    '''
    set_offset_str='/{:s}/sigouts/%d/offset' % int(channel_I)
    daq.setDouble(set_offset_str.format(device), off_I)

def set_off_Q(daq, device, off_Q = 0,  channel_Q = 1):
    '''
    set output offset for Q signal

    daq:                daq ID
    device:             device ID
    off_Q:              offset of waveform ouTput for Q signal
    channel_Q:          output channel for Q signal
    '''
    set_offset_str='/{:s}/sigouts/%d/offset' % int(channel_Q)
    daq.setDouble(set_offset_str.format(device), off_Q)

def set_dA(daq, device, dA_IQ_set = 0.0, amp = 100e-3, channel_I = 0, channel_Q = 1):
    '''
    set relative amplitude difference of sine wave amplitude for I and Q signals

    daq:                daq ID
    device:             device ID
    dA_IQ_set:          relative amplitude of sine wave amplitude for I and Q signal in %
    amp:                amplitude of sine wave in Volt
    channel_I:          output channel for I signal
    channel_Q:          output channel for Q signal
    '''

    dA_IQ = (dA_IQ_set)/100.0
    ampl1 = (1 + dA_IQ)*amp
    ampl2 = (1 - dA_IQ)*amp
    set_dA_I_str = '/{:s}/sines/%d/amplitudes/0' % int(channel_I)
    set_dA_Q_str = '/{:s}/sines/%d/amplitudes/1' % int(channel_Q)
    daq.setDouble(set_dA_I_str.format(device), float(ampl1))

    daq.setDouble(set_dA_Q_str.format(device), float(ampl2))

def config_output(daq, device, osc_ch = 0, osc_f = 0, sines_ch = [0, 1], sines_phaseshift = [0, 90], sines_amplitude = [0, 0],\
                  sines_enable = [0, 0], awg_group = 1, sr_exp_reduction = 0, outp_range = [2, 2], outp_offset = [0, 0], \
                  modula=[0,0]):


    '''
config output of HDAWG

daq:                daq ID
device:             device ID
osc_ch:             oscillator channel
osc_f:              oscillator frequency
sines_ch:           sine signal channels for I and Q
sines_phaseshift:   phase shift of sine signals
sines_amplitude:    amplitudes of sine signals
sines_enable:       disable/enable sine signals
awg_group:     0:   4 x 2 or (2 x 2) channel grouping, 4 sequencers
               1:   2 x 4 or (1 x 4) channel grouping, 2 sequencers
               2:   1 x 8 only for HD8, 1 sequencer
sr_exp_reduction:   n (int n from 0 to 13): sampling rate setting base_rate/2^n
outp_range:         range of HDAWG waveform output in V
outp_offset:        offset of HDAWG waveform output in V
'''
    daq.setInt(f'/{device:s}/system/awg/channelgrouping', awg_group)

    if awg_group == 0:
        daq.setInt('/{:s}/awgs/0/time'.format(device), sr_exp_reduction)

    for i in sines_ch:
        sines_osc_str        = '/{:s}/sines/%d/oscselect'     % i
        sines_phaseshift_str = '/{:s}/sines/%d/phaseshift'    % i
        sines_amplitude_str  = '/{:s}/sines/%d/amplitudes/%d' % (i, i)
        outp_offset_str      = '/{:s}/sigouts/%d/offset'      % i
        outp_range_str       = '/{:s}/sigouts/%d/range'       % i
        outp_on_str          = '/{:s}/sigouts/%d/on'          % i
        sines_enable_str     = '/{:s}/sines/%d/enables/%d'    % (i, i)
        if int(i % 2) == 0:
            osc_f_str        = '/{:s}/oscs/%d/freq'           % int(np.floor(i/2))


        daq.setInt(sines_osc_str.format(device), osc_ch)
        daq.setDouble(osc_f_str.format(device),   osc_f)
        daq.setDouble(sines_phaseshift_str.format(device), sines_phaseshift[i])
        daq.setDouble(sines_amplitude_str.format(device),   sines_amplitude[i])
        daq.setDouble(outp_offset_str.format(device),  outp_offset[i])
        daq.setDouble(outp_range_str.format(device),    outp_range[i])
        daq.setInt(outp_on_str.format(device),                      1)
        daq.setInt(sines_enable_str.format(device),   sines_enable[i])
        daq.setInt(f'/{device}/awgs/0/outputs/{i}/modulation/mode', modula[i])


def enable_awg(daq, device, enable = 1, single = 1, awgs = [0]):
    '''
    enable/disable AWG

    daq:               daq ID
    device:            device ID
    enable:     0/1:   disable/enable AWG
    awg_return:   0:   AWG sequence run in return mode
            1:   AWG sequence run in single mode

    '''
    [daq.setInt(f'/{device:s}/awgs/{core:d}/single', single) for core in awgs]
    [daq.syncSetInt(f'/{device:s}/awgs/{core:d}/enable', enable) for core in awgs] # Run the HDAWG AWG


# # HDAWG-PC option

# ## enable_precompensation

# In[ ]:



def enable_precompensation(daq, device, tc_exponential = [[20e-9]], amp_exponential = [[-30e-3]], tc_highpass = [919e-9], clearing_highpass = [1], delay_bounce = [0], amp_bounce = [0], parameters_fir = [[510e-3, 0, 250e-3, 0, 250e-3]], exponentials_num = [1], channels = [0], enable_precompensation = 1, enable_exponentials = 1, enable_highpass = 1, enable_bounce = 0, enable_fir = 1):

    '''
precompensation setting

daq:                   daq ID
device:                device ID
tc_exponential:        time constant in second of exponential filter
amp_exponential:       amplitude of exponential filter
tc_highpass:           time constant in second of highpass filter
clearing_highpass: 0:  level
                   1:  rising
                   2:  falling
                   3:  both rising and falling
delay_bounce:          delay in second of bounce filter
amp_bounce:            amplitude of bounce filter
parameters_fir:        FIR parameters
exponentials_num:      number of exponential filter
channels:              which channels using precompensation
enable_precompensation:
                   0:  disable
                   1:  enable
enable_exponentials:   0/1: disable/enable
enable_highpass:       0/1: disable/enable
enable_bounce:         0/1: disable/enable
enable_fir:            0/1: disable/enable
'''

    for i in range(len(channels)):
        daq.setInt('/{:s}/sigouts/%d/precompensation/status/reset'.format(device) % channels[i], 1)
        if enable_exponentials == 1:
            for k in range(exponentials_num[i]):

                daq.setDouble('/%s/sigouts/%d/precompensation/exponentials/%d/timeconstant'                            % (device, channels[i], k), tc_exponential[i][k])

                daq.setDouble('/%s/sigouts/%d/precompensation/exponentials/%d/amplitude'                            % (device, channels[i], k), amp_exponential[i][k])
                daq.setInt('/%s/sigouts/%d/precompensation/exponentials/%d/enable'                            % (device, channels[i], k), enable_exponentials)
        if enable_highpass == 1:

            daq.setDouble('/%s/sigouts/%d/precompensation/highpass/0/timeconstant'                            % (device, channels[i]), tc_highpass[i])
            daq.setInt('/%s/sigouts/%d/precompensation/highpass/0/clearing/slope'                            % (device, channels[i]), clearing_highpass[i])
            daq.setInt('/%s/sigouts/%d/precompensation/highpass/0/enable'                            % (device, channels[i]), enable_highpass)
        if enable_bounce == 1:
            daq.setInt('/%s/sigouts/%d/precompensation/bounces/0/enable'                            % (device, channels[i]), enable_bounce)
            daq.setDouble('/%s/sigouts/%d/precompensation/bounces/0/amplitude'                            % (device, channels[i]), amp_bounce[i])
            daq.setDouble('/%s/sigouts/%d/precompensation/bounces/0/delay'                            % (device, channels[i]), delay_bounce[i])

        if enable_fir == 1:

#            coef_fir[0:len(parameters_fir[i]), channels[i]]=np. transpose(parameters_fir[i])

#            for j in range(len(parameters_fir[i])):
            daq.setVector('/%s/sigouts/%d/precompensation/fir/coefficients'                        % (device, channels[i]), np.asarray(parameters_fir[i]))
            daq.setInt('/%s/sigouts/%d/precompensation/fir/enable'                            % (device, channels[i] ), enable_fir)
        daq.setInt('/%s/sigouts/%d/precompensation/enable'  % (device, channels[i]), enable_precompensation)


# # HDAWG-CNT option

# ## config_pusle_counters

# In[ ]:


def config_pusle_counters(daq, device, cnt = 0, cnt_mode = 0, cnt_period = 100e-6, cnt_input = 32, cnt_gate = 0,                           trigger_edge = [1, 0], cnt_operation = 0, cnt_integrate = 0, enable = 1):
    '''
    config pulse counters

    daq:                   daq ID
    device:                device ID
    cnt:                   pulse counter ID, range from 0 to 7
    cnt_mode:              run mode
                     1:    Free runing
                     2:    gated free running
                     3:    gated
                     4:    time tagging
    cnt_period:            period used for for the free running and gated running modes also sets the hold-off
                           time for the time tagging mode
    cnt_input:             counter signal source
    cnt_gate:              signal source used for enabling the counter in the Gated Free Running and Gated modes
                           Allowed Values:
                           0 Trigger/Ref Input 1 (front panel).
                           1 Trigger/Ref Input 2 (front panel).
                           2 Trigger Input 3 (rear panel).
                           3 Trigger Input 4 (rear panel).
                           4 AWG Trigger 1.
                           5 AWG Trigger 2.
                           6 AWG Trigger 3.
                           7 AWG Trigger 4.
    trigger_edge:          trigger edge
                 [0/1,0/1]:[disable/enable rise, disable/enable fall]
    cnt_operation:         Select the arithmetic operation (addition, subtraction) applied to the counter unit outputs.\
                           'Other counter' refers to the grouping of the counter units: 1 with 2, and 3 with 4.
                           Allowed Values:
                           0 None
                           1 Add Other Counter
                           2 Subtract Other Counter
    cnt_integrate:         Sum up counter values over time
    enable:                0/1  disable/enable the pulse counter
    '''
    daq.setInt(f'/{device:s}/cnts/{cnt:d}/mode', cnt_mode)
    daq.setDouble(f'/{device:s}/cnts/{cnt:d}/period', cnt_period)
    daq.setInt(f'/{device:s}/cnts/{cnt:d}/inputselect', cnt_input)
    daq.setInt(f'/{device:s}/cnts/{cnt:d}/gateselect', cnt_gate)
    daq.setInt(f'/{device:s}/cnts/{cnt:d}/trigrising', trigger_edge[0])
    daq.setInt(f'/{device:s}/cnts/{cnt:d}/trigfalling', trigger_edge[1])
    daq.setInt(f'/{device:s}/cnts/{cnt:d}/operation', cnt_operation)
#     daq.setInt(f'/{device:s}/cnts/{cnt:d}/sample', sample)
    daq.setInt(f'/{device:s}/cnts/{cnt:d}/integrate', cnt_integrate)
    daq.setInt(f'/{device:s}/cnts/{cnt:d}/enable', enable)


# ## data_acquisation_cnts

# In[ ]:


def data_acquisation_cnts(daq, device, cnt = 0, time_recording = 1, timeout = 100):
    '''
    data acquisation for pulse counters

    daq:                   daq ID
    device:                device ID
    cnt:                   pulse counter ID, range from 0 to 7
    time_recording:        recording time in second
    timeout:               timeout in second
    '''

    #path = f'/{device:s}/cnts/{cnt:d}/value'
    path = f'/{device:s}/cnts/{cnt:d}/sample'

    daq.subscribe(path)
#     daq.sync() clear buffer and the
    data = daq.poll(time_recording, timeout, flat=True)
    return data
# def daq_module()
# daq_module = ziDAQ('dataAcquisitionModule');
# ziDAQ('set', daq_module, 'triggernode', '/dev8027/triggers/streams/0/sample.TrigIn1');
# ziDAQ('set', daq_module, 'preview', 1);
# ziDAQ('set', daq_module, 'device', 'dev8027');
# ziDAQ('set', daq_module, 'historylength', 100);
# ziDAQ('set', daq_module, 'bandwidth', 0);
# ziDAQ('set', daq_module, 'hysteresis', 0.01);
# ziDAQ('set', daq_module, 'level', 0.1);
# ziDAQ('set', daq_module, 'type', 6);
# ziDAQ('set', daq_module, 'save/directory', 'C:\Users\chunyans\Documents\Zurich Instruments\LabOne\WebServer');
# ziDAQ('set', daq_module, 'clearhistory', 1);
# ziDAQ('set', daq_module, 'clearhistory', 1);
# ziDAQ('set', daq_module, 'bandwidth', 0);
# ziDAQ('set', daq_module, 'clearhistory', 1);
# ziDAQ('set', daq_module, 'bandwidth', 0);
# daq_hd,device_hd = create_api_sessions_hd('DEV8027')
# config_pusle_counters(daq_hd, device_hd, cnt = 0, cnt_mode = 0, cnt_period = 100e-6, cnt_input = 32, cnt_gate = 0, \
#                           trigger_edge = [1, 0], cnt_operation = 0, cnt_integrate = 0, enable = 1)
# def multi_awg_seq(daq, device, index = [0], experiment = 'Rabi'):

# def set_MFM(daq, device, index)

# def set_sweeper()



# In[ ]:


# import numpy as np
# import matplotlib.pyplot as plt
# cnt = 0
# daq, device = create_api_sessions_hd('dev8027')
# awg_seq(daq, device, result_length = 1, amp_q = 0.1, \
#             time_increment = 50e-6, sample_exponent = 5, time_origin_sec = 160e-6, \
#             period_wait_sec_pulse = 400e-6, period_wait_sec_spec = 10e-6, \
#             exp = 'T2')

# # please find the setting detail in p217-219
# trigger_edge = [0, 1]
# daq.setInt(f'/{device:s}/cnts/{cnt:d}/mode', 3)
# daq.setInt(f'/{device:s}/cnts/{cnt:d}/gateselect', 1)
# daq.setInt(f'/{device:s}/cnts/{cnt:d}/inputselect', 32)
# daq.setDouble(f'/{device:s}/cnts/{cnt:d}/period', 1e-3)
# daq.setInt(f'/{device:s}/cnts/{cnt:d}/trigrising', trigger_edge[0])
# daq.setInt(f'/{device:s}/cnts/{cnt:d}/trigfalling', trigger_edge[1])
# daq.setInt(f'/{device:s}/cnts/{cnt:d}/integrate', 10)
# daq.setInt(f'/{device:s}/cnts/{cnt:d}/operation', 0)
# daq.setInt(f'/{device:s}/cnts/{cnt:d}/enable', 1)
# N = 1000
# a = np.zeros(N)
# b=daq_hd.record(f'/{device_hd:s}/cnts/0/value', )
# # for i in range(N):
# #     a[i] = daq_hd.getInt(f'/{device_hd:s}/cnts/0/value')

# plt.plot(a, '.')
# help(daq_hd)


# In[ ]:


# daq_hd, device_hd = create_api_sessions_hd('dev8027')
# path = f'/{device_hd:s}/cnts/0/value'
# daq_hd.subscribe(path)
# while True:
#     ret = daq_hd.poll(1.0, 100, flat=True)
#     if path in ret:
#         print(ret[path])




# In[ ]:



# exp_str = ['qubit spec', 'Rabi', 'Rabi_t', 'Rabi_n2pi', 'Rabi_2p_t', 'T1', 'T2 Ramsey', 'T2 Hanh echo', 'T2 CPMG', \
#            'reset qubit','precompensation_demo']

# sequence = exp_str[5]
# daq_hd, device_hd=create_api_sessions_hd('dev8138')
# awg_seq(daq_hd, device_hd, result_length = 50, amp_q = 0.5, \
#             time_increment = 50e-9, sample_exponent = 0, time_origin_sec = 10e-6, \
#             period_wait_sec_pulse = 10e-6, period_wait_sec_spec = 10e-6, \
#             exp = sequence,trigger_type = 0, wave_dur_sec = 100e-9, n_pi_CPMG=10, seqr_index = 0, timeout = 5)


# In[ ]:


# daq_hd, device_hd=create_api_sessions_hd('dev8138')
# sweep_str = ['sweep N', 'sweep amp', 'sweep t', 'sweep phase']
# n = 3
# sequence = sweep_str[n]
# if n == 0:
#     nSteps = [5, 5]
#     start   = [0, 0]
#     step    = [1, 1]
# elif n == 1:
#     nSteps = [5, 5]
#     start   = [0, 0]
#     step    = [0.1, 0.1]
# elif n == 2:
#     nSteps = [5, 5]
#     start   = [32, 48]
#     step    = [16, 16]
# elif n == 3:
#     nSteps = [5, 5]
#     start   = [0, 0]
#     step    = [30, 30]

# awg_seq_cz(daq_hd, device_hd, nSteps= nSteps, period_wait_sec_pulse = 1e-6, time_origin_sec = 1e-6, \
#                     ampl = 0.5, n_dc =2, ampl_dc1 = 0.1, ampl_dc2 = 0.2,  t1 = 32, t2 = 48, phase = 0,\
#                     start = start, step = step, wave_dur_sec = 100e-9, exp = sequence, timeout = 10)


# In[ ]:


# import zhinst.utils as ziut
# import matplotlib.pyplot as plt
# n_combination=2000
# max_samples = 12000
# period_wait_sec = 100e-9

# wave_long, marker_long=false_rb_seq_generation1(amplx=0.5, amply=0.3, samples_x= 200, samples_y=200, max_samples=max_samples, n_combination=n_combination,\
#                              n_avg_seq = 8, max_n_seq=16)
# print(len(wave_long))
# print('false rb seq generated')


# In[ ]:


# plt.close()
# plt.figure()
# plt.plot(wave_long[0:12000*20])
# plt.plot(marker_long[0:12000*20])


# In[ ]:


# timeout=50
# daq_hd, device_hd=create_api_sessions_hd('dev8138')
# awg_seq_rb_1q(daq_hd, device_hd, n_combination=n_combination, max_samples=max_samples, period_wait_sec=period_wait_sec, timeout=timeout)
# new_wave = wave_long
# marker = marker_long
# new_wave1 = ziut.convert_awg_waveform(wave1=new_wave, markers=marker)
# daq_hd.setVector(f'/{device_hd}/awgs/0/waveform/waves/0', new_wave1)


# In[ ]:


# import time as time
# daq_hd, device_hd=create_api_sessions_hd('dev8138')
# max_length=12000
# n_combination=20

# awg_seq_dynm_demo(daq_hd, device_hd, n_repeat=1, max_samples=max_samples,\
#                   period_wait_sec=100e-9, timeout=2)
# enable_awg(daq_hd, device_hd, enable=1, single=1)
# print('set awg seq')
# for i in range(n_combination):
#     wave_sub   = wave_long[i*max_length:(i+1)*max_length]
#     marker_sub = marker_long[i*max_length:(i+1)*max_length]
#     new_wave1 = ziut.convert_awg_waveform(wave1=wave_sub, markers=marker_sub)
#     a=time.clock()
#     daq_hd.setVector(f'/{device_hd}/awgs/0/waveform/waves/0', new_wave1)
#     enable_awg(daq_hd, device_hd, enable=0, single=1)
#     b=time.clock()

# time_spend = (b-a)*n_combination
# print(time_spend)


# In[33]:


# daq_hd, device_hd = create_api_sessions_hd('dev8027')
# new_wave = ziut.convert_awg_waveform(0.2*np.ones(128))
# new_wave_1 = ziut.convert_awg_waveform(0.5*np.ones(128))

# daq_hd.setVector(f'/{device_hd}/awgs/0/waveform/waves/0', new_wave)
# daq_hd.setVector(f'/{device_hd}/awgs/0/waveform/waves/1', new_wave_1)
# enable_awg(daq_hd, device_hd, enable=1, single=0)


# In[ ]:


# samples = 100000000
# seqc=textwrap.dedent('''playWave(zeros(100000000));''')
# create_and_compile_awg(daq_hd,device_hd,seqc,50)
# new_wave = np.ones(samples)*0.5
# new_wave1 = zhinst.utils.convert_awg_waveform(new_wave)
# daq_hd.setVector(f'/{device_hd}/awgs/0/waveform/waves/0', new_wave1)

