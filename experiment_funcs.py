"""
Created on Thu Apr 14 11:35:33 2022

@author: Evangelos Vlachos <evlachos@usc.edu>
"""

# # Import Modules

import time
import sys
import matplotlib.pyplot as plt
import numpy as np
import UHFQA as qa
from tqdm import tqdm
import HDAWG as hd
#import USB6501 as usb
import zhinst.utils as ziut
import comTablefuncs as ctfuncs
import json
import csv
# import keyboard as kb
import os
# from plot_functions_nb import Volt2dBm,Watt2dBm
import itertools
from VISAdrivers.sa_api import *
import collections
from scipy import optimize,fft
# import statsmodels.api as sm
# from PyTektronixScope import PyTektronixScope
import seaborn as sns; sns.set() # styling

sns.set_style('ticks')


def setup_active_reset(awg,daq,threshold=5):
    """
    Sets up active reset. The discrimination threshold has to be previously established via a Rabi Measurement.

    Args:
        awg (class): API class instance of HDAWG.
        daq (class): API class instance of UHFQA.
        threshold (int, optional): The value that is going to be used to discriminate the measurement results. Defaults to 5.

    Returns:
        None.

    """
    # Configure AWG settings
    # select trigger sources
    daq.setInt('/dev8233/awgs/0/auxtriggers/0/channel', 0)
    daq.setInt('/dev8233/awgs/0/auxtriggers/1/channel', 1)
    # Select trigger slope. First trigger is QA Result Trigger (rise), second is QA Result (level)
    daq.setInt('/dev8233/awgs/0/auxtriggers/0/slope', 1)
    daq.setInt('/dev8233/awgs/0/auxtriggers/1/slope', 0)
    # sets trigger level
    daq.setDouble('/dev8233/triggers/in/0/level', 0.3)
    daq.setDouble('/dev8233/triggers/in/1/level', 0.3)

    #Configure QA settings
    # select trigger sources
    daq.setInt('/dev2528/triggers/out/0/source', 74) # QA result trigger. The QA sends a trigger to HDAWG Ch. 1 when measurement is done.
    daq.setInt('/dev2528/triggers/out/1/source', 64) # QA result. Sends trigger to HDAWG based on the measurement result.
    # set trigger mode to output ("drive")
    daq.setInt('/dev2528/triggers/out/0/drive', 1)
    daq.setInt('/dev2528/triggers/out/1/drive', 1)
    # set trigger levels to 3 V
    # daq.setDouble('/dev2528/triggers/in/0/level', 3)
    # daq.setDouble('/dev2528/triggers/in/1/level', 3)
    # sets QA result threshold
    daq.setDouble('/dev2528/qas/0/thresholds/0/level', threshold)
    daq.sync()

def readoutSetup(daq,sequence='time',readout_pulse_length=1.2e-6,rr_IF=5e6,cav_resp_time=0.5e-6):
    """
    Configures the AWG of the UHFQA for readout.

    Args:
        qa (class): API class instance of QA.
        sequence (string, optional): Whether we are executing spectroscopy (spec) or time measurement. Defaults to 'time'.
        readout_pulse_length (float, optional): Readout pulse length in seconds. Defaults to 1.2e-6.
        rr_IF (float, optional): Intermediate mod/demod frequency in Hz. Defaults to 5e6.
        cav_resp_time (float, optional): The response time of the cavity in seconds. Defaults to 0.5e-6.

    Returns:
        None.

    """

    fs = 450e6
    readout_amp = 0.7
    # prepares QA for readout | loads readout Sequence into QA AWG and prepares QA
    print('-------------Setting up Readout Sequence-------------')
    if sequence =='spec':
        qa.awg_seq_readout(daq, readout_length=readout_pulse_length,rr_IF=rr_IF,nPoints=roundToBase(readout_pulse_length*fs),base_rate=fs,cav_resp_time=cav_resp_time,amplitude_uhf=readout_amp)
        daq.setInt('/dev2528/awgs/0/time',2)
    elif sequence =='time':
        qa.awg_seq_readout(daq, readout_length=readout_pulse_length,rr_IF=rr_IF,nPoints=roundToBase(readout_pulse_length*fs),base_rate=fs,cav_resp_time=cav_resp_time,amplitude_uhf=readout_amp)
        daq.setInt('/dev2528/awgs/0/time',2)

def pulsed_spec_setup(daq,awg,nAverages,qubit_drive_amp,mu=0,sigma=0,qubit_drive_dur=20e-6,result_length=1,integration_length=2e-6,nPointsPre=0,nPointsPost=0,delay=500,measPeriod=100e-6):
    '''Setup HDAWG and UHFQA for pulsed spectroscopy'''
    hd.awg_seq(awg,sequence='qubit spec',fs=0.6e9,nAverages=nAverages,qubit_drive_dur=qubit_drive_dur,measPeriod=measPeriod,mu=mu,sigma=sigma,amp_q= qubit_drive_amp,nPointsPre=nPointsPre,nPointsPost=nPointsPost)
    awg.setInt('/dev8233/awgs/0/time',2) # sets AWG sampling rate to 600 MHz
    qa.config_qa(daq,sequence='spec',integration_length=integration_length,nAverages=1,result_length=result_length,delay=delay)

def single_shot_setup(daq,awg,nAverages=1024,qubit_drive_amp=0.1,mu=0,sigma=0,fs=0.6e9,result_length=1,integration_length=2e-6,pi2Width=100e-9,measPeriod=400e-6):
    '''Setup HDAWG for single shot experiment'''
    pi2Width = int(pi2Width*fs)
    print('-------------Setting HDAWG sequence-------------')
    hd.awg_seq(awg,sequence='single_shot',fs=fs,nAverages=nAverages,mu=mu,sigma=sigma,amp_q=qubit_drive_amp,measPeriod=measPeriod,pi2Width=pi2Width)
    awg.setInt('/dev8233/awgs/0/time',2) # sets AWG sampling rate to 600 MHz

def ssb_setup(awg,qb_IF=50e6):

    awg.setDouble('/dev8233/oscs/0/freq', qb_IF)
    awg.setDouble('/dev8233/sines/1/phaseshift', 90)
    awg.setInt('/dev8233/awgs/0/outputs/0/modulation/mode', 3)
    awg.setInt('/dev8233/awgs/0/outputs/1/modulation/mode', 4)
    awg.setDouble('/dev8233/sigouts/0/range', 0.6)
    awg.setDouble('/dev8233/sigouts/1/range', 0.6)

def IQ_imbalance(awg,g,phi):
    """
    Calculates the correction matrix (C**-1) and applies it to the AWG

    Args:
        awg: awg instance to apply correction
        g (TYPE): amplitude imbalance.
        phi (TYPE): phase imbalance.

    Returns:
        list: correction matrix.

    """

    c = np.cos(phi)
    s = np.sin(phi)
    N = 1 / ((1-g**2)*(2*c**2-1))
    corr = np.array([float(N * x) for x in [(1-g)*c, (1+g)*s, (1-g)*s, (1+g)*c]]).reshape(2,2)

    for i in range(2):
        for j in range(2):
            awg.setDouble(f'/dev8233/awgs/0/outputs/{j}/gains/{i}', corr[i,j])

def seq_setup(awg,sequence='rabi',nAverages=128,nPoints=1024,pulse_length_start=32,
              fs=2.4e9,nSteps=100,pulse_length_increment=16,Tmax=0.3e-6,amp_q=1,sigma=0,pi2Width=0,piWidth_Y=0,
              pipulse_position=20e-9,measPeriod=200e-6,instance=0,B0=0,active_reset=False,axis='X',noise_rate=1,qb_IF=50e6):
    """
    Function Description
    --------------------

    Sets up the right sequence to the AWG along with the corresponding command table. It also sets up Single Sideband Modulation.

    Parameters
    ----------
    awg : class
        awg to write the sequence to
    'dev8233' : string
        serial number of AWG. Default is dev8233
    sequence : string
        Type of experiment. Available options are 'rabi','ramsey','echo', and 'T1'. The default is 'rabi'.
    nAverages : integer
    nPoints : integer
        Number of points in the waveform. The default is 1024.
    fs : float
        Sampling rate of the AWG. The default is 2.4e9.
    nSteps : integer
        Number of points in the sequence. The default is 100.
    pulse_length_increment : integer
        dt in AWG samples (1/fs). The default is 16.
    Tmax : float
    active_reset: Boolean
    amp_q : TYPE, optional
        qubit drive and B0 (if applicable) amplitude. The default is 1.
    pi2Width : TYPE, optional
        pi/2 duration The default is 50e-9.
    measPeriod : TYPE, optional
        Waiting interval between measurements. Must be at least 2*T_1. The default is 200e-6.
    qb_IF:    float, optional
        Qubit intermediate modulation frequency in Hz. Default is 50e6.

    """
    fs_base = 2.4e9

    awg.setDouble('/dev8233/oscs/0/freq', qb_IF) # set the IF frequency of the modulator
    # Generate and compile program
    print('-------------Setting HDAWG sequence-------------')
    bt = time.time()
    awg.setInt('/dev8233/awgs/0/time',(int(fs_base/fs-1))) # set sampling rate of AWG
    hd.awg_seq(awg,sigma=sigma,B0=B0,axis=axis,fs=fs,nSteps=nSteps,nPoints=nPoints,
               pi2Width=pi2Width, Tmax=Tmax,amp_q=amp_q,nAverages=nAverages,
               sequence=sequence,measPeriod=measPeriod,active_reset=active_reset)
    et = time.time()
    print('HDAWG compilation duration: %.1f s'%(et-bt))

    if sequence == 'echo_v2':
        ct = {}
    else:
        # create and upload command table
        ct=ctfuncs.ct_pulse_length(n_wave=nSteps, pulse_length_start=pulse_length_start, pulse_length_increment=pulse_length_increment,
                                   pipulse=2*pi2Width,active_reset=active_reset,sequence=sequence,noise_rate=noise_rate)
        awg.setVector("/dev8233/awgs/0/commandtable/data", json.dumps(ct))

    awg.sync()
    return nSteps,ct

def single_shot(daq,awg,cav_resp_time=1e-6,measPeriod=400e-6,integration_length=2.3e-6,mu=0,sigma=0,rr_IF=30e6,pi2Width=100e-9,qubit_drive_amp=1,readout_drive_amp=0.1,setup=0,nAverages=128):
    '''
    DESCRIPTION: Executes single shot experiment.

    '''
    result_length =  2*nAverages
    fsAWG = 600e6
    base_rate = 1.8e9

    readout_pulse_length = integration_length + cav_resp_time + 1e-6

    if not setup:
        single_shot_setup(daq,awg,pi2Width=pi2Width,result_length=result_length,fs=fsAWG,mu=mu,sigma=sigma,integration_length=integration_length,nAverages=nAverages,qubit_drive_amp=qubit_drive_amp,measPeriod=measPeriod)
        readoutSetup(daq,sequence='spec',readout_pulse_length=readout_pulse_length,cav_resp_time=cav_resp_time)
        time.sleep(0.1)

    sweep_data, paths = qa.create_sweep_data_dict(daq, 'dev2528')
    data_pi = []
    data_OFF = []

    bt = time.time()
    qa.qa_result_reset(daq, 'dev2528')
    hd.enable_awg(awg, 'dev8233',enable=0,awgs=[0])
    qa.config_qa(daq,sequence='single shot',nAverages=1,integration_length=integration_length,result_length=result_length,delay=cav_resp_time)
    daq.sync()

    qa.qa_result_enable(daq, 'dev2528')
    qa.enable_awg(daq, 'dev2528') # start the readout sequence
    hd.enable_awg(awg,'dev8233',enable=1,awgs=[0])

    print('Start measurement')
    data = qa.acquisition_poll(daq, paths, result_length, timeout = 3*nAverages*measPeriod) # transfers data from the QA result to the API for this frequency point
    # seperate OFF/ON data and average
    data_OFF = np.append(data_OFF, [data[paths[0]][k] for k in even(len(data[paths[0]]))])/(integration_length*base_rate)
    data_pi =  np.append(data_pi, [data[paths[0]][k] for k in odd(len(data[paths[0]]))])/(integration_length*base_rate)


    hd.enable_awg(awg, 'dev8233',enable=0,awgs=[0])
    qa.stop_result_unit(daq, 'dev2528', paths)
    qa.enable_awg(daq, 'dev2528', enable = 0)

# ----------------------------------------------------------------------------------
    et = time.time()
    duration = et-bt
    print(f'Measurement time: %.1f s'%duration)
#-----------------------------------

    return data_OFF,data_pi

def spectroscopy(daq,awg,qubitLO=0,cav_resp_time=1e-6,integration_length=2e-6,AC_pars=[0.0,0],qubit_drive_amp=1,
                 readout_drive_amp=0.1,setup=True,nAverages=128,frequencies=np.linspace(3.7,3.95,1001)):
    '''
    DESCRIPTION: Executes qubit spectroscopy.

    '''
    result_length =  2*nAverages
    fsAWG = 600e6
    base_rate = 1.8e9
    mu = AC_pars[0]
    sigma = AC_pars[1]

    readout_pulse_length = integration_length + cav_resp_time + 2e-6
    # daq.setDouble('/dev2528/sigouts/0/amplitudes/0', readout_drive_amp)
    nPointsPre = nPointsPost = roundToBase(500e-9*fsAWG,base=16)
    qubit_drive_dur = roundToBase(60e-6*fsAWG,base=16)
    if setup:
        pulsed_spec_setup(daq, awg, result_length=result_length,mu=mu,sigma=sigma,qubit_drive_dur=qubit_drive_dur,integration_length=integration_length,nAverages=nAverages,qubit_drive_amp=qubit_drive_amp,nPointsPre=nPointsPre,nPointsPost=nPointsPost,delay=cav_resp_time)
        readoutSetup(daq, sequence='spec',readout_pulse_length=readout_pulse_length,cav_resp_time=cav_resp_time)
        time.sleep(0.1)

    # initialize signal generators and set power

    print('Start measurement')
    sweep_data, paths = qa.create_sweep_data_dict(daq, 'dev2528')
    data_ON = []
    data_OFF = []

    qa.enable_awg(daq, 'dev2528') # start the readout sequence
    bt = time.time()
    j = 0
    for f in frequencies:
        qubitLO.RF_OFF()
        qubitLO.set_freq(f)
        qubitLO.RF_ON()
        qa.qa_result_reset(daq, 'dev2528')
        qa.qa_result_enable(daq, 'dev2528')
        hd.enable_awg(awg,'dev8233',awgs=[0]) #runs the drive sequence
        data = qa.acquisition_poll(daq, paths, result_length, timeout = 60) # transfers data from the QA result to the API for this frequency point
        # seperate OFF/ON data and average
        data_OFF = np.append(data_OFF, np.mean([data[paths[0]][k] for k in even(len(data[paths[0]]))]))
        data_ON =  np.append(data_ON, np.mean([data[paths[0]][k] for k in odd(len(data[paths[0]]))]))

        sys.stdout.write('\r')
        sys.stdout.write(f'progress:{int((j+1)/len(frequencies)*100)}%')
        sys.stdout.flush()
        j = j + 1


    data = (data_ON-data_OFF)/(integration_length*base_rate)
    I_data= data.real
    Q_data = data.imag

    power_data = np.abs(I_data*I_data.conjugate()+Q_data*Q_data.conjugate())

    hd.enable_awg(awg, 'dev8233',enable=0,awgs=[0])
    qa.stop_result_unit(daq, 'dev2528', paths)
    qa.enable_awg(daq, 'dev2528', enable = 0)

# ----------------------------------------------------------------------------------
    et = time.time()
    duration = et-bt
    print(f'Measurement time: {duration} s')
#-----------------------------------

    return power_data,I_data,Q_data

def pulse(daq,awg,setup=[0,0,0],Tmax=0.3e-6,nSteps=61,prePulseLength=1500e-9,postPulseLength=1500e-9,nAverages=128,amp_q=1,
          sequence='rabi',mu=0,sigma=0,stepSize=2e-9, B0=0,tauk=0,nu=0,cav_resp_time=0.5e-6,piWidth_Y=0,AC_freq=5e-9, source=2,verbose=0,
          pipulse_position=20e-9,integration_length=2.3e-6,qubitDriveFreq=3.8135e9,pi2Width=0,rr_IF = 30e6,sampling_rate=1.2e9,
          measPeriod=300e-6,sweep=0,instance=0,active_reset=False,threshold=500e-3,white_noise_instance=[],
          noise_instance=[],phi=0,axis='X',wk=0,noise_rate=1):

    '''
    DESCRIPTION:            Runs a single pulsed experiment (Rabi,Ramsey,T1,Echo)
    -----------------------------------------------------------------------------------------------------------------
    setup[0]:                   If setup[0]=0, the right seqc programs are loaded into the HDAWG. If setup[0] = 2, the noise waveforms are substituted and there is no compilation (used for sweeps or statistics)
    setup[1]:                   If setup[1]=0, the right seqc programs are loaded into the QA AWG for readout
    setup[2]:                   If setup[2]=0, QA is configured
    Tmax:                       max length of drive (rabi) or pi2 pulse separation
    amp_q:               amplitude of qubit drive channel
    sequence:                   Which experiment to perform (see description)
    pi2Width:                   Length of pi2 pulse in seconds
    instance:                   Which instance of telegraph noise to use. Used for sweeps
    rr_IF:                      The IF of the readout mixer
    'pipulse_position':         Where to insert the pipulse (only applicable for echo with telegraph noise). A higher number means the pipulse is applied sooner
    cav_resp_time:              The time it takes for the cavity to ring up/down. Since we are using square integration pulse, we don't want to include the edges of the pulse
    integration_length:         How long the QA integrates for. The readout pulse is 2 microseconds longer than integration+cavity_response
    phi:                        Whether to include a random phase in the cos term of the noise signal (only applicable for noiseType II)
    '''
    fs = sampling_rate
    base_rate = 1.8e9       # sampling rate of QA (cannot be changed in standard mode)
    readout_pulse_length = integration_length + cav_resp_time + 0.5e-6
    if sweep == 0:
        verbose = 1
    cores = [0]

    # stops AWGs and reset the QA
    hd.enable_awg(awg,'dev8233',enable=0,awgs=cores)
    qa.enable_awg(daq, 'dev2528',enable=0)
    qa.qa_result_reset(daq, 'dev2528')

    if setup[0] == 0:

        nPoints,nSteps,pulse_length_increment,pulse_length_start = calc_nSteps(sequence=sequence,fsAWG=fs,
                                                                               stepSize=stepSize,Tmax=Tmax,verbose=verbose)

        if sequence == 'ramsey' or sequence == 'echo' or sequence == 'T1':
            create_wfm_file(awg,mu,sigma,B0,nu,tauk, noise_rate=noise_rate,sequence=sequence, nPoints=nPoints, white_noise_instance=white_noise_instance,
                            noise_instance=noise_instance,Tmax=Tmax,phi=phi,sweep=sweep, wk=wk)

        if sequence == 'echo_v2':
            create_echo_wfms(awg,fs,mu,sigma,B0,nPoints,Tmax,amp_q,pi2Width,nSteps,pulse_length_increment)

        nSteps,ct = seq_setup(awg,noise_rate=noise_rate,sequence=sequence,axis=axis,piWidth_Y=piWidth_Y,pipulse_position=pipulse_position,nSteps=nSteps,
                              nPoints=nPoints,fs=fs,pulse_length_start=pulse_length_start,pulse_length_increment=pulse_length_increment,
                              instance=instance,amp_q=amp_q,nAverages=nAverages,pi2Width=pi2Width,
                              Tmax=Tmax,sigma=sigma,B0=B0,measPeriod=measPeriod,active_reset=active_reset)

    elif setup[0] == 2:
        bt = time.time()
        # replace waveforms, don't recompile program
        nPoints,nSteps,pulse_length_increment,pulse_length_start = calc_nSteps(sequence=sequence,fsAWG=fs,piWidth_Y=piWidth_Y,
                                                                               stepSize=stepSize,Tmax=Tmax,verbose=verbose)
        if B0 == 0:
            noise_instance = np.zeros(nPoints)
        if mu != 0 or sigma != 0:
            white_noise = np.random.normal(mu, sigma, nPoints)
            waveforms_native = ziut.convert_awg_waveform(wave1=noise_instance,wave2=white_noise)
        else:
            waveforms_native = ziut.convert_awg_waveform(wave1=noise_instance)
        path = '/dev8233/awgs/0/waveform/waves/0'
        awg.setVector(path,waveforms_native)
        awg.sync()
        et = time.time()
        print('Replacing waveforms took: %.1f ms'%(1e3*(et-bt)))

    if setup[1] == 0:
        readoutSetup(daq,readout_pulse_length=readout_pulse_length,sequence='pulse',rr_IF=rr_IF,cav_resp_time=cav_resp_time) # setup QA AWG for readout
    if setup[2] == 0:
        qa.config_qa(daq,sequence='pulse',nAverages=nAverages,rr_IF=rr_IF,integration_length=integration_length,result_length=nSteps,delay=cav_resp_time,source=source) # configure UHFQA result unit
        daq.sync()
    if active_reset == True:
        setup_active_reset(awg, daq,threshold=threshold)

    # Determine whether command table is used for error checking later on
    if mu != 0 and sequence != 'echo_v2':
        use_ct = 1
    elif mu == 0 and (sequence == 'rabi' or (sequence == 'ramsey' and sigma != 0)):
        use_ct = 1
    else:
        use_ct = 0

    measTime = calc_timeout(nAverages, measPeriod, stepSize, nSteps)
    if sweep == 0:
        print('Estimated Measurement Time (with/without active reset): %.1f/%.1f sec'%(int(1/8*measTime),measTime))
    else:
        pass

    # Checks whether the right command table was uploaded to the HDAWG
    ct_awg = json.loads(daq.get("/dev8233/awgs/0/commandtable/data",flat=True)["/dev8233/awgs/0/commandtable/data"][0]['vector'])
    if setup[0] == 0 and use_ct == 1:
        if ct_awg != ct:
            print('Error! Invalid Command Table used for Measurement\nCommand Table Sent to AWG\n\n%s\n\nCommand Table in AWG\n\n%s'%(ct,ct_awg))
            sys.exit()

    if sequence == 'ramsey':
        result_length = nSteps-1
    else:
        result_length = nSteps
    if active_reset == False:
        timeout = 1.2*measTime
    elif active_reset == True:
        timeout = 0.2*1.2*measTime

    sweep_data, paths = qa.create_sweep_data_dict(daq, 'dev2528')

    qa.enable_awg(daq, 'dev2528',enable=1) # start the readout sequence

    qa.qa_result_enable(daq, 'dev2528') # arm the qa

    str_meas = time.time()
    hd.enable_awg(awg,'dev8233',enable=1,awgs=cores) #runs the drive sequence
    data = qa.acquisition_poll(daq, paths, num_samples = result_length, timeout = timeout) # retrieve data from UHFQA

    for path, samples in data.items():
        sweep_data[path] = np.append(sweep_data[path], samples)

    # reset QA result unit and stop AWGs
    qa.stop_result_unit(daq, 'dev2528', paths)
    hd.enable_awg(awg, 'dev8233', enable = 0,awgs=cores)
    qa.enable_awg(daq, 'dev2528', enable = 0)

    end_meas = time.time()
    if sweep == 0:
        print('\nmeasurement duration: %.1f s' %(end_meas-str_meas))
    else:
        pass


    data = sweep_data[paths[0]][0:result_length]/(integration_length*base_rate)
    if source == 2 or source == 1:
        results = data
    elif source == 7:
        I = data.real
        Q = data.imag
        results = [[I],[Q]]

    #Generate time array points using the command table (if applicable)
    t = np.zeros(result_length)
    if use_ct == 1:
        for i in range(result_length):
            t[i] = ct_awg['table'][i]['waveform']['length']/fs
    else:
        t = np.linspace(pulse_length_start/fs,Tmax,len(data))

    if sequence=='echo' or sequence == 'echo_v2':
        t = 2*t

    return t,results,nSteps

def calc_nSteps(sequence='ramsey',fsAWG=1.2e9,stepSize=10e-9,Tmax=5e-6,verbose=1):
    """
    Calculates the number of steps in the sequence and number of points in the waveform used in the experiment. The final stepsize of the sequence
    might be different than the one specified due to the fact that the AWG has a granularity of the waveform is 16 samples.

    Args:
        sequence (string, optional): Type of experiment (rabi,ramsey,T1, echo). Defaults to 'ramsey'.
        fsAWG (float, optional): HDAWG sampling rate. Defaults to 1.2e9.
        stepSize (float, optional): dt for sequence. Defaults to 10e-9.
        Tmax (float, optional): maximum experiment array length. For example if Tmax = 5e-6 for ramsey experiment then maximum pulse separation is 5us. Defaults to 5e-6.
        verbose (boolean, optional):  Defaults to 1.

    Returns:
        nPoints (int): number of waveform points.
        nSteps (int): number of sequence steps.
        pulse_length_increment (int): sequence step size in units of samples 1/fsAWG.
        pulse_length_start (int): starting sequence point in units of samples.

    """

    if sequence == 'rabi':
        base = 4 # this increases the sampling rate for Rabi beyond the 16-sample limit imposed by the AWG. This is only possible for Rabi as we don't care about zero padding
    elif sequence == 'echo_v2':
        base = 1
    else:
        base = 16
    pulse_length_start = int(roundToBase(stepSize*fsAWG,base=base))
    # if pulse_length_start < 32 and sequence != 'rabi':
    #     print('Smallest Waveform Length is 32 samples. The first point in this sequence has %d samples'%(pulse_length_start))
        # sys.exit()
    # if sequence == 'echo' and B0 != 0:
    #     pulse_length_start = roundToBase(64+int(piWidth_Y*fsAWG))
    #     # if pulse_length_start % 16 != 0:
    #     #     print('Pulse length start point (%d) is not multiple of 16'%pulse_length_start)
    #     #     sys.exit()
    # elif sequence == 'rabi':
    #     pulse_length_start = 32
    # else:
    #     pulse_length_start = 32
    pulse_length_increment = roundToBase(fsAWG*stepSize,base=base)
    nPoints = roundToBase(Tmax*fsAWG,base=pulse_length_increment) # this ensures there is an integer number of time points
    nSteps = int((nPoints-pulse_length_start)/pulse_length_increment) + 1 # 1 is added to include the first point
    # if sequence == 'echo':
    #   pulse_length_increment = pulse_length_increment / 2
    if verbose == 1:
        print("dt is %.1f ns (%d pts) ==> f_s = %.1f MHz \nNpoints = %d | n_steps is %d | Pulse length start = %.1f ns (%d pts)" %(pulse_length_increment/fsAWG*1e9,pulse_length_increment,1e-6*fsAWG/pulse_length_increment,nPoints,nSteps,pulse_length_start*1e9/fsAWG,pulse_length_start))
    else:
        pass

    if nSteps > 1024:
        raise Exception('Error: The maximum number of steps is 1024')

    return nPoints,nSteps,pulse_length_increment,pulse_length_start

def pull_wfm(sweep_name,nu,tauk,sequence='ramsey'):
    """
    Retrieves waveform instances from csv file. Only applicable for parameter sweeps

    Args:
        sweep_name (str): sweep identification number.
        nu (float): memory kernel modulation frequency in Hz.
        tauk (float): mean switching time in seconds.
        sequence (str, optional): type of experiment. Options are 'ramsey' and 'echo'. Defaults to 'ramsey'.

    Returns:
        noise_realizations (array): 2D waveform array.

    """

    path = "E:\\generalized-markovian-noise\\%s\\sweep_data\\%s\\%s\\noise_instances"%('CandleQubit_6',sequence,sweep_name)
    filename = 'nu_%d_Hz_tau_%d_ns.csv' %(round(nu*1e3),round(tauk*1e3))
    # filename = 'RTN_tau_%d_ns.csv' %(round(tau*1e3))
    print(os.path.join(path,filename))
    noise_realizations = np.loadtxt(os.path.join(path,filename),dtype=float,delimiter=',')

    return noise_realizations

def create_wfm_file(awg,mu,sigma,B0,nu,tauk,sequence="ramsey",nPoints=1000,white_noise_instance=[],sweep=0,noise_instance=[],Tmax=5e-6,phi=0,meas_device='CandleQubit_6',wk=0,noise_rate=1):
    '''Generates waveform text file to be used by the HDAWG sequence file. The waveforms are either generated by the function or passed as an argument (used for sweeps)'''
    # create RTN noise or pull instance from file (only for parameter sweeps)
    if B0 != 0 and len(noise_instance) == 0:
        print('Generating New Waveform')
        t = np.linspace(0,Tmax,nPoints)
        if wk == 0:
            qubit_free_evol = B0 * np.cos(2*np.pi*nu*1e3*t+phi*2*np.pi*np.random.rand()) * gen_tel_noise(nPoints, tauk, dt=Tmax/nPoints)
        else:
            qubit_free_evol = B0 * gen_WK_sig(fs=2*nPoints/Tmax, nu=nu*1e3, tauk=tauk*1e-6, Tmax=Tmax)
    elif B0 != 0 and len(noise_instance) > 0:
        print('Loading Waveform from csv File')
        qubit_free_evol = noise_instance
    else:
        qubit_free_evol = np.zeros(nPoints)

    qubit_free_evol = qubit_free_evol[...,None] # transpose waveforms such that they are columns (necessary such that they are readable by AWG seqc files)

    # create white noise instance
    if mu != 0 or sigma != 0:
        if sweep == 0:
            white_noise = np.random.normal(loc=mu, scale=sigma, size=nPoints)
        elif sweep == 1:
            white_noise = white_noise_instance

        white_noise = white_noise[...,None]
        wfm_arr = np.hstack((qubit_free_evol,white_noise)) # stack waveform columns horizontally
    else:
        wfm_arr = qubit_free_evol
    # save file
    np.savetxt("C:/Users/LFL/Documents/Zurich Instruments/LabOne/WebServer/awg/waves/"+sequence+"_wfm"+".csv", wfm_arr, delimiter = ",") # save file where it can be called by the AWG sequence program

def gen_WK_sig(fs,nu,tauk,Tmax):

    '''generate waveform using Wiener-Kinchin method'''

    N = int(fs*Tmax)+1
    dt = 1/fs
    df = 1/Tmax
    t = np.linspace(-Tmax/2,Tmax/2,N)

    # define autocorrelation and compute PSD
    autocorr = np.cos(2*np.pi*nu*t)*np.exp(-np.abs(t)/tauk)
    psd = 2/N*fft.fft(autocorr)
    freqs = fft.fftfreq(N,dt*1e6)[:round(N/2)]
    power = np.sqrt(np.abs(psd*df))
    # freq = np.fft.fftfreq(N,dt)[:round(2*nu/df)]
    # plt.plot(freq,psd[:round(2*nu/df)])
    # plt.show()

    # generate array of random phases
    phi_arr = np.random.rand(int(len(power)/2)) * 2*np.pi
    P_psd = np.zeros(len(phi_arr),dtype=complex)

    for i in range (1,len(phi_arr)):
        P_psd[i] = power[i]*np.exp(1j*phi_arr[i])

    # construct the 2 sided, symmetric fourier spectrum
    rand_psd = np.concatenate(([power[0]],np.conjugate(P_psd),np.flip(P_psd)))
    # inverse fourier transform to get timeseries
    signal = np.fft.ifft(rand_psd)
    signal = signal[int(len(signal)/2)+1:]

    if np.random.rand() < 0.5:
        signal = - signal
    else:
        pass

    return np.real(signal)/max(np.real(signal)),np.abs(psd[:round(N/2)]),freqs,autocorr[round(N/2)+1:]

def calc_autocorr(sig):
    '''Calculates the autocorrelation of the given signal'''
    return sm.tsa.acf(sig,nlags=len(sig))

def gen_noise_realizations(par1_arr=np.linspace(0,10,100),par2_arr=[0],numRealizations=3,nPoints=1000,T_max=5e-6,sweep_count=1,
                           meas_device='CandleQubit_6',sequence='ramsey',wk=False,plot=False):
    """
    Generates noise realizations and saves them to a csv file for parameter sweep

    Args:
        par1_arr (array, optional): array of tauk values. Defaults to np.linspace(0,10,100).
        par2_arr (array, optional): arary of nu values. Defaults to [0].
        numRealizations (TYPE, optional): DESCRIPTION. Defaults to 3.
        nPoints (int, optional): DESCRIPTION. Defaults to 1000.
        T_max (float, optional): DESCRIPTION. Defaults to 5e-6.
        sweep_count (int, optional): DESCRIPTION. Defaults to 1.
        meas_device (TYPE, optional): DESCRIPTION. Defaults to 'CandleQubit_6'.
        sequence (str, optional): DESCRIPTION. Defaults to 'ramsey'.
        wk (boolean, optional): whether to use Wiener-Kinchin method to generate waveforms. Defaults to 0.
        plot (boolean, optional): whether to plot noise waveforms. Defaults to 0.

    Returns:
        None.

    """

    if len(par2_arr) > 1 or par2_arr[0] != 0:
        phi = 1
    else:
        phi = 0
    numPoints_par1 = len(par1_arr)
    numPoints_par2 = len(par2_arr)
    t = np.linspace(0,T_max,nPoints)
    parent_dir = 'E:\\generalized-markovian-noise\\%s\\sweep_data\\%s\\'%(meas_device,sequence)
    directory = 'sweep_%03d\\noise_instances'%(sweep_count)
    path = os.path.join(parent_dir,directory)
    noise_arr = np.zeros((numRealizations,nPoints))
    for i in range(numPoints_par2):
        for k in range(numPoints_par1):
            filename = "nu_%d_Hz_tau_%d_ns.csv" % (round(par2_arr[i]*1e3),round(par1_arr[k]*1e3))
            with open(os.path.join(path,filename),"w",newline="") as datafile:
                writer = csv.writer(datafile)
                for j in range(numRealizations):
                    if len(par2_arr) > 1 or par2_arr != 0:
                        if wk:
                            noise_arr[j,:] = np.cos(2*np.pi*par2_arr[i]*1e3*t + phi*2*np.pi*np.random.rand()) * gen_tel_noise(nPoints, par1_arr[k], dt = T_max/nPoints)
                        elif wk:
                            noise,psd,freqs,autocorr = gen_WK_sig(fs=2*nPoints/T_max, nu=par2_arr[i]*1e3, tauk=par1_arr[k]*1e-6, Tmax=1e-3)
                            noise_arr[j,:] = noise[:nPoints]/max(noise[:nPoints])
                    elif len(par2_arr) <= 1 and par2_arr[0] == 0:
                        noise_arr[j,:] = gen_tel_noise(nPoints, par1_arr[k], dt = T_max/nPoints)

                writer.writerows(noise_arr)

    if plot:
        fig = plt.figure(figsize=(4,8),dpi=300)
        ax1 = fig.add_subplot(3,1,1) # noise realization plot
        ax2 = fig.add_subplot(3,1,2) # mean autocorrelation plot
        ax3 = fig.add_subplot(3,1,3) # PSD plot (real & imag)

        ac = np.zeros(nPoints)
        # compute autocorrelations and average over noise realizations
        for i in range(numRealizations):
            ac += calc_autocorr(noise_arr[i,:])
        ac = ac/numRealizations

        ax1.plot(t[:100]*1e6,noise_arr[1,:100])
        ax1.set_ylabel('$\sigma_x(t)$')
        ax1.set_title('Noise Realization - $\\tau_k$ = %.1f $\mu$s | $\\nu$ = %d kHz'%(par1_arr[0],par2_arr[0]))

        ax2.plot(t[:100]*1e6,autocorr[:100],'r',label='Analytic')
        ax2.plot(t[:100]*1e6,ac[:100],'b',label='Numeric')
        ax2.set_title('Autocorrelation')
        ax2.set_xlabel('Time ($\mu$s)')
        ax2.legend()

        ax3.plot(freqs[0:10000],np.real(psd)[0:10000],'ro',label='Re',linestyle='None')
        ax3.plot(freqs[0:10000],np.imag(psd)[0:10000],'bx',label='Im',linestyle='None')
        ax3.set_xlabel('Frequency(MHz)')
        ax3.legend()
        ax3.set_title('PSD')
        fig.tight_layout()
        plt.show()

def gen_tel_noise(numPoints,tau,dt):
    '''Generates a single instance of telegraph noise'''
    signal = np.ones(numPoints)*(-1)**np.random.randint(0,2)
    for i in range(1,numPoints-1):
        if np.random.rand() < 1/(2*tau*1e-6/dt)*np.exp(-1/(2*tau*1e-6/dt)):
            signal[i+1] = - signal[i]
        else:
            signal[i+1] = signal[i]
    return signal

def roundToBase(nPoints,base=16):
    '''Make the AWG happy by uploading a wfm whose points are multiple of 16'''
    y = base*round(nPoints/base)
    if y==0:
        y = base*round(nPoints/base+1)
    return y

def convertB0_kHz_to_mV(B0_kHz):
    return (B0_kHz-6.75)/62.1

def convertB0_mV_to_kHz(B0_mV):
    return 62.1*B0_mV+6.75

def odd(n):
    return range(1,n,2)

def even(n):
    return range(0,n,2)

def qa_monitor_avg(daq,length,averages):

    settings = [
        ("qas/0/monitor/enable", 0),
        ("qas/0/monitor/length", length),
        ("qas/0/monitor/averages", averages),
        ("qas/0/monitor/enable", 1),
        ("qas/0/monitor/reset", 1),
        ('qas/0/monitor/trigger/channel', 7)
    ]
    daq.set([(f"/{'dev2528'}/{node}", value) for node, value in settings])
    # Signals to measure

    paths = []
    for channel in range(2):
        path = f"/{'dev2528':s}/qas/0/monitor/inputs/{channel:d}/wave"
        paths.append(path)
    daq.setInt('/dev2528/qas/0/monitor/reset', 1)
    daq.setInt('/dev2528/qas/0/monitor/enable', 1)

    daq.sync()
    time.sleep(0.2)
    # daq.setInt('dev2528/qas/0/monitor/trigger/channel',1)
    daq.subscribe(paths)

    # Perform acquisition
    print("Acquiring data...")
    data = qa.acquisition_poll(daq, paths, length,timeout=60)
    daq.setInt('/dev2528/qas/0/monitor/enable', 0)
    fig = plt.figure(figsize=(12,6))
    plt.plot(data[paths[0]])
    plt.plot(data[paths[1]])
    print(len(data[paths[0]]))
    print(len(data[paths[1]]))
    plt.title(f'Input signals after {averages:d} averages')
    plt.xlabel('nPoints')
    plt.ylabel('Amp (V)')
    plt.grid()
    plt.show()

    return data

def scope_meas(awg,daq,length=8192,nAverages=128,samp_rate=1.8e9,trigLevel=0.1):

    '''Executes a measurement with the UHFQA'''
    #setup and initialize scope
    scope = qa.config_scope(daq,'dev2528',scope_length=length,scope_avg_weight=1,scope_mode=0)
    scope = daq.scopeModule()
    daq.setInt('/dev2528/scopes/0/channel', 3)# 1: ch1; 2(DIG):ch2; 3: ch1 and ch2(DIG)
    daq.setInt('/dev2528/scopes/0/channels/0/inputselect', 0) # 0: sigin 1; 1: sigin 2
    daq.setDouble('/dev2528/scopes/0/length', length)
    scope.set('scopeModule/averager/weight', 1)
    scope.set('scopeModule/mode', 2)
    daq.setDouble('dev2528/scopes/0/length', length*nAverages)
    daq.setInt('/dev2528/scopes/0/single', 1)
    daq.setInt('/dev2528/scopes/0/trigchannel', 0)
    daq.setInt('/dev2528/scopes/0/trigenable', 1)
    daq.setDouble('/dev2528/scopes/0/trigholdoff', 50e-6) #In units of second. Defines the time before the trigger is rearmed after a recording event
    daq.setDouble('/dev2528/scopes/0/triglevel', trigLevel) #in Volts
    # daq.setInt('/dev2528/scopes/0/segments/enable', 1)
    daq.setInt('/dev2528/scopes/0/time', int(np.log2(1.8e9/samp_rate))) #set sampling rate
    daq.setInt('/dev2528/scopes/0/channel', 3) # enables both signal inputs
    # daq.setInt('/dev2528/scopes/0/segments/count',nAverages)
    daq.setDouble('/dev2528/sigins/0/range',0.4)
    daq.sync()

    qa.restart_avg_scope(scope)
    qa.enable_scope(daq,'dev2528',enable=1)
    qa.subscrib_scope(scope,'dev2528')
    qa.execute_scope(scope)

    hd.enable_awg(awg,'dev8233',enable=1,awgs=[0])

    while int(scope.progress()) != 1:
        # time.sleep(0.05)
        result = scope.read()

    ch1Data = result['%s' % 'dev2528']['scopes']['0']['wave'][0][0]['wave'][0]/2**15
    ch2Data = result['%s' % 'dev2528']['scopes']['0']['wave'][0][0]['wave'][1]/2**15
    avgCh1Data = np.zeros(length)
    avgCh2Data = np.zeros(length)

    # for k in range(nAverages):
    #     avgCh1Data = avgCh1Data + ch1Data[k:k+length]
    #     avgCh2Data = avgCh2Data + ch2Data[k:k+length]

    # lines = plt.plot(avgCh1Data)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(ch1Data)
    ax.plot(ch2Data)

    qa.enable_scope(daq,'dev2528',enable=0)
#    scope.set('scopeModule/clearhistory', 0)
    qa.finish_scope(scope)

    # return scope,avgCh1Data,avgCh2Data
    return scope,ch1Data,ch2Data

def calc_timeout(nAverages,measPeriod,dt,nSteps):
    '''Calculates timeout for experiment. If measurement takes more than timeout seconds, then the program stops and gives an error'''
    t = 0
    for i in range(nSteps):
        t += (dt*i+measPeriod)*nAverages
    return t

def init_arrays(numRealizations=128,interval=2,nPointsBackground=200,nPoints=200):
    '''Initializes arrays to be used for storing exp data during parameter sweep'''
    bData = np.zeros((int(numRealizations/interval),nPointsBackground),dtype=float)
    data = np.zeros((numRealizations,nPoints),dtype=float)
    return bData,data

# def set_AWG_output_amplitude(range)

def calc_sweep_time(par1,par2,measTimeBackground=1,measTime=25,nMeasBackground=100,nMeas=100):
    return (measTimeBackground*nMeasBackground+measTime*nMeas)*len(par1)*len(par2)


def create_echo_wfms(awg,fs=1.2e9,mu=0,sigma=0,B0=0,nPoints=1024,Tmax=5e-6,amp=0,pi2Width=50e-9,nSteps=101,pulse_length_increment=32):

    '''
    DESCRIPTION: Generates a series of waveforms to be uploaded into the AWG. The output is a series of csv files.
    '''

    start = time.time()
    ACpre = mu*np.ones(roundToBase(1500e-9*fs))
    pi2 = amp*np.ones(int(fs*pi2Width))
    pi2pre = 0 * ACpre
    ac_noise = np.random.normal(mu, sigma, 2*nPoints)
    tel_noise = np.zeros(2*nPoints)
    for i in range(nSteps):
        ch1_wfm = np.concatenate((pi2pre,pi2,tel_noise[0:i*pulse_length_increment],pi2,pi2,tel_noise[i*pulse_length_increment:2*i*pulse_length_increment],pi2,pi2pre))
        ch2_wfm = np.concatenate((ACpre,mu*pi2/amp,ac_noise[0:i*pulse_length_increment],mu*pi2/amp,mu*pi2/amp,ac_noise[i*pulse_length_increment:2*i*pulse_length_increment],mu*pi2/amp,ACpre))
        ch1_wfm = ch1_wfm[...,None]
        ch2_wfm = ch2_wfm[...,None]
        wfm_2D_arr = np.hstack((ch1_wfm,ch2_wfm))
        np.savetxt("C:/Users/LFL/Documents/Zurich Instruments/LabOne/WebServer/awg/waves/"+"echo_"+"wfm_%03d"%(i)+".csv", wfm_2D_arr, delimiter = ",")

    end = time.time()
    print('Generating echo Waveforms took %.1f' %(end-start))

def rabi_ramsey(daq,awg,qubitLO,qubitDriveFreq=3.8e9,AC_pars=[0,0],plot=1):
    '''---------------------------------------Do Rabi, calibrate pi-pulse, then do Ramsey----------------------------------------------'''

    optionsRabi= {
        'nAverages':        128,
        'Tmax':             0.2e-6,
        'amp_q':     1.0,
        'sequence':         'rabi',
        'channel':          0,
        'measPeriod':       200e-6,
        'qubitDriveFreq':   qubitDriveFreq,
        'AC_pars':          AC_pars
        }

    t,ch1Data,ch2Data,nPoints = pulse(daq,awg,qubitLO,setup=[0,0,0],**optionsRabi)
    pi_pulse,error = pf.pulse_plot1d(sequence='rabi',dt=optionsRabi['Tmax']*1e6/nPoints,qubitDriveFreq=optionsRabi['qubitDriveFreq'],amp_q=optionsRabi['amp_q'],x_vector=t, y_vector=ch1Data,fitting=1,AC_pars=optionsRabi['AC_pars'],plot=0)

    optionsRamsey = {
        'nAverages':        128,
        'Tmax':             10e-6,
        'stepSize':         100e-9,
        'pi2Width':         1/2*pi_pulse*1e-9,
        'amp_q':     optionsRabi['amp_q'],
        'sequence':         'ramsey',
        'channel':          0,
        'measPeriod':       200e-6,
        'qubitDriveFreq':   optionsRabi['qubitDriveFreq'],
        'AC_pars':          optionsRabi['AC_pars'],
        'RT_pars':          [0,0]

        }


    t,ch1Data,ch2Data,nPoints = pulse(daq,awg,qubitLO,setup=[0,0,0],**optionsRamsey)
    # plot data
    detuning,T_phi,error = pf.pulse_plot1d(sequence='ramsey',dt=optionsRamsey['Tmax']*1e6/nPoints,qubitDriveFreq=optionsRamsey['qubitDriveFreq'],amp_q=optionsRamsey['amp_q'],x_vector=t, y_vector=ch1Data,fitting=1,AC_pars=optionsRamsey['AC_pars'],pi2Width=optionsRamsey['pi2Width'],RT_pars=optionsRamsey['RT_pars'],plot=plot)

    return detuning,T_phi,error


def snr(sa,fc,thres):
    """
    Calculates the SNR

    Args:
        sa (class): The Spectrum Analyzer Object
        fc (double): The frequency of the signal in Hz
        thres (int): The reference level of the spectrum analyzer

    Returns:
        snr (double) The SNR.

    """

    # configure SA
    sa_config_acquisition(device = sa, detector = SA_AVERAGE, scale = SA_LOG_SCALE)
    sa_config_center_span(sa, fc, 0.5e6)
    sa_config_level(sa, thres)
    sa_config_gain_atten(sa, SA_AUTO_ATTEN, SA_AUTO_GAIN, True)
    sa_config_sweep_coupling(device = sa, rbw = 1e3, vbw = 1e3, reject=0)

    # Initialize SA
    sa_initiate(sa, SA_SWEEPING, 0)
    query = sa_query_sweep_info(sa)
    sweep_length = query["sweep_length"]
    start_freq = query["start_freq"]
    bin_size = query["bin_size"]

    freqs = np.array([start_freq + i * bin_size for i in range(sweep_length)],dtype=float)

    signal = sa_get_sweep_64f(sa)['max']
    plt.plot(1e-9*freqs,signal)
    plt.xticks(np.linspace(min(1e-9*freqs), max(1e-9*freqs),5))
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Power (dBm)')
    plt.show()

    max_ind = np.argmax(signal)
    max_val = np.max(signal)
    mask = np.logical_or (freqs < freqs[max_ind]-10e3, freqs > freqs[max_ind]+10e3)
    noisetemp = signal[mask]
    avg_noise = np.mean(noisetemp)
    snr = max_val-avg_noise


    print("SNR: %.1f\nNoise Floor: %.1f dBm"%(snr,avg_noise))

    return snr

def condition(x): return x > 5

def line(x,a,b):
    return a*x+b

def save_data(device_name,project,exp='rabi',iteration=1,par_dict={},data=[]):

    dir_path = f'D:\\{project}\\{device_name}\\{exp}-data'
    filename = f'\\{exp}'+f'_data_{iteration}.csv'

    if os.path.exists(dir_path):
        pass
    else:
        print(f'Directory not found; making new directory at {dir_path}')
        os.makedirs(dir_path)

    with open(dir_path+filename,"w",newline="") as datafile:
        writer = csv.writer(datafile)
        writer.writerow(par_dict.keys())
        writer.writerow(par_dict.values())
        writer.writerow(data[0,:])
        writer.writerow(data[1,:])