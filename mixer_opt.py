# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 14:55:08 2022

@author: Evangelos Vlachos evlachos@usc.edu

Mixer optimization functions
"""


# # Import Modules

import time
import tqdm
import plot_functions as pf
import numpy as np
#import USB6501 as usb
import zhinst.utils as ziut
import comTablefuncs as ctfuncs
# import keyboard as kb
from VISAdrivers.sa_api import *
import collections
from scipy import optimize
from HDAWG import awg_seq
from UHFQA import setup_mixer_calib

def get_power(sa,inst,fc=4e9,threshold=-50,config=False,plot=False,output=False):
    """
    Measures the power at a specific frequency using the spectrum analyzer. Can calculate the ON/OFF ratio if desired.

    Args:
        sa (class): The API class instance corresponding to the spectrum analyzer.
        inst (class): The API class instance corresponding to the HDAWG or UHFQA.
        fc (double): The frequency at which we want to measure the power.
        threshold (int, optional): The reference level for the SA. Defaults to -50.
        plot (boolean, optional): To plot or not the data. Defaults to False.

    Returns:
        OFF_power (double): The power (leakage) at the frequency specified by fc.

    """

    # configure SA
    if config:
        sa_config_acquisition(device = sa, detector = SA_AVERAGE, scale = SA_LOG_SCALE)
        sa_config_center_span(sa, fc, 0.5e6)
        sa_config_gain_atten(sa, SA_AUTO_ATTEN, SA_AUTO_GAIN, True)
        sa_config_sweep_coupling(device = sa, rbw = 1e3, vbw = 1e3, reject=0)
        sa_config_level(sa, threshold)
        sa_initiate(sa, SA_SWEEPING, 0)
        query = sa_query_sweep_info(sa)
        sweep_length = query["sweep_length"]
        start_freq = query["start_freq"]
        bin_size = query["bin_size"]
        freqs = np.array([start_freq + i * bin_size for i in range(sweep_length)],dtype=float)

    #get OFF power (leakage)
    signal = sa_get_sweep_64f(sa)['max']
    power = np.max(signal)

    if plot:
        pf.power_plot(freqs, signal, power, fc=fc)
        if output:
            print(f'{power} dBm at {fc/1e9} GHz')
    return power

def config_sa(sa,fc,span=5e6,reference=-30):
        """
        Prepares spectrum analyzer for measurement
        Parameters
        ----------
        sa :
            Handle for spectrum analyzer
        freq : float
            Center frequency of span.
        span : float, optional
            DESCRIPTION. The default is 5e6.
        reference : float, optional
            Upper power threshold of SA in dBm. The default is -30.
        Returns
        -------
        None.
        """

        sa_config_level(sa, reference) # sets sensitivity
        sa_config_center_span(sa, freq, span) # sets center frequency
        sa_initiate(sa, SA_SWEEPING, 0)
        query = sa_query_sweep_info(sa)
        sweep_length = query["sweep_length"]
        start_freq = query["start_freq"]
        bin_size = query["bin_size"]
        freqs = np.array([start_freq + i * bin_size for i in range(sweep_length)],dtype=float)

def min_leak(sa,inst,device='dev8233',mode='fine',mixer='qubit',threshold=-50,f_LO=3.875e9,amp=0.2,channels=[0,1],measON=False,plot=False):
    """

    DESCRIPTION:
        Optimizes mixer at given frequency

    INPUTS:
        sa (class): API class instance of spectrum analyzer.
        inst (class): The class instance of the instrument (HDAWG or UHFQA) that controls the I and Q channels of the mixer we want to optimize.
        mode(string): Coarse or fine tuning. In "coarse" tuning mode the parameter steps are large.
        mixer (string): The mixer we want to optimize. Options are "qubit","resonator",and "stark". Defaults to 'qubit'.
        f_LO (float): The local oscillator (LO) frequency of the mixer. Defaults to 3.875e9.
        f_IF (float): The intermediate (IF) frequency of the mixer. Defaults to 50e6.
        amp (float): Amplitude of ON Pulse.
        channels (list): The AWG channel used for I/Q in the experimental setup.
        measON (boolean): Whether or not to measure the ON power of the mixer.
        plot (boolean): Whether or not to plot the leakage as a function the parameters.
    """

    start = time.time()
    if mode == 'coarse':
        span=20e-3
        dV=1e-3
    elif mode == 'fine':
        span=2e-3
        dV=0.1e-3

    # generate arrays for optimization parameters
    vStart = np.zeros(2)
    for i in range(len(vStart)):
        vStart[i] = inst.get(f'/{device}/sigouts/{channels[i]}/offset')[f'{device}']['sigouts'][f'{channels[i]}']['offset']['value']
        inst.sync()
    VoltRange1 = np.arange(vStart[0]-span/2,vStart[0]+span/2,dV)
    VoltRange2 = np.arange(vStart[1]-span/2,vStart[1]+span/2,dV)

    vStart[i] = inst.get(f'/{device}/sigouts/{channels[i]}/offset')[f'{device}']['sigouts'][f'{channels[i]}']['offset']['value']
    inst.sync()

    VoltRange1 = np.arange(vStart[0]-span/2,vStart[0]+span/2,dV)
    VoltRange2 = np.arange(vStart[1]-span/2,vStart[1]+span/2,dV)

    L1 = len(VoltRange1)
    L2 = len(VoltRange2)
    power_data = np.zeros((L1,L2))

    config_sa(sa,fc=f_LO,reference=threshold)

    # Sweep individual channel voltages and find leakage
    with tqdm(total = L1*L2) as progress_bar:
        for i,V1 in enumerate((VoltRange1)):
            for j,V2 in enumerate((VoltRange2)):
                inst.setDouble(f'/{device}/sigouts/{channels[0]}/offset',V1)
                inst.setDouble(f'/{device}/sigouts/{channels[1]}/offset',V2)
                inst.sync()
                power_data[i,j] = get_power(sa, inst, f_LO,threshold=threshold,plot=False,config=False)

    # find index of voltage corresponding to minimum LO leakage
    argmin = np.unravel_index(np.argmin(power_data), power_data.shape)

    opt_I = VoltRange1[argmin[0]]
    opt_Q = VoltRange2[argmin[1]]
    # set voltages to optimal values
    inst.set('/%s/sigouts/%d/offset'%(device,channels[0]),opt_I)
    inst.sync()
    inst.set('/%s/sigouts/%d/offset'%(device,channels[1]),opt_Q)
    inst.sync()
    print(f'optimal I_offset = {round(opt_I*1e3,1)} mV, optimal Q_offset = {round(1e3*opt_Q,1)} mV')

    end = time.time()
    print('%s mixer Optimization took %.1f seconds'%(mixer,(end-start)))

    # get LO leakage for optimal DC values
    OFF_power = get_power(sa, inst, f_LO,threshold=threshold,plot=False)

    if measON:
        offset = inst.get(f'/{device}/sigouts/{channels[0]}/offset')[f'{device}']['sigouts'][f'{channels[0]}']['offset']['value'][0]
        #get ON power
        sa_config_level(sa, 0)
        sa_initiate(sa, SA_SWEEPING, 0)
        inst.set(f'/{device}/sigouts/{channels[0]}/offset', amp)
        inst.sync()
        signal_ON = sa_get_sweep_64f(sa)['max']
        ON_power = np.max(signal_ON)
        inst.set(f'/{device}/sigouts/{channels[0]}/offset', offset)
        inst.sync()
    else:
        pass

    if device == 'dev8233':
         element = 'qubit'
    elif device == 'dev2528':
         element = 'readout'

    if plot:
        pf.plot_mixer_opt(VoltRange1, VoltRange2, power_data,cal='LO',element=element,fc=fc)

def suppr_image(sa,inst,mode='fine',mixer='qubit',threshold=-50,f_LO=3.875e9,f_IF=50e6,channels=[0,1],sb='lsb',gen=0,plot=False):
        """

        DESCRIPTION:
            Minimizes power at sideband at given frequency. The default sideband we want to minimize is the lower sideband for now.
            In later versions we can add the ability to choose which sideband to optimize.

        INPUTS:
            sa (class): API class instance of spectrum analyzer.
            inst (class): The class instance of the instrument (HDAWG or UHFQA) that controls the I and Q channels of the mixer we want to optimize.
            mode(string): Coarse or fine tuning. In "coarse" tuning mode the parameter steps are large.
            mixer (string): The mixer we want to optimize. Options are "qubit","resonator",and "stark". Defaults to 'qubit'.
            f_LO (float): The local oscillator (LO) frequency of the mixer. Defaults to 3.875e9.
            f_IF (float): The intermediate (IF) frequency of the mixer. Defaults to 50e6.
            channels (list): The AWG channel connected to the I port of the mixer you want to optimize.
            sb (str): Which sideband is the image. Default is the lower sideband ('lsb')
            gen (int): The sine generator used for modulation.
            plot (boolean): Whether or not to plot the image power as a function the parameters.
        """

        if sb == 'lsb':
            f_im = f_LO - f_IF
        elif sb == 'usb':
            f_im = f_LO + f_IF

        start = time.time()
        if mode == 'coarse':
            span=20e-3
            dp = 1
            da = 1e-3
        elif mode == 'fine':
            span=2e-3
            dp = 0.1
            da = 0.1e-3

        if str(inst) == 'awg':
            device = 'dev8233'
        elif str(inst) == 'daq':
            device = 'dev2528'

        # get current values of phase and amplitude
        p0 = inst.get(f'/{device}/sines/1/phaseshift')[f'{device}']['sines'][f'{gen}']['phaseshift']['value']
        a0 = inst.get(f'/{device}/awgs/0/outputs/0/gains/0')[f'{device}']['awgs']['0']['outputs'][f'{channels[0]}']['gains'][f'{channels[1]}']['value'][0]
        # generate arrays for optimization parameters based on current values of phi and a used
        phiArr = np.arange(p0-span/2,p0+span/2,dp)
        ampArr = np.arange(a0-span/2,a0+span/2,da)

        imPower1 = np.zeros(len(pArr))
        imPower2 = np.zeros(len(ampArr))

        inst.enable_awg(inst,device,enable=0) # stop AWG
        # upload and run AWG sequence program
        if device == 'dev8233':
            awg_seq(inst,sequence='mixer_calib')

        elif device == 'dev2528':
            setup_mixer_calib(inst)

        inst.enable_awg(inst,device,enable=1)

        L1 = len(phiArr)
        L2 = len(ampArr)
        power_data = np.zeros((L1,L2))

        config_sa(sa,fc=f_LO,reference=threshold)

        # Sweep individual channel voltages and find leakage
        with tqdm(total = L1*L2) as progress_bar:
            for i,amp in enumerate((ampArr)):
                for j,phi in enumerate((phiArr)):
                    inst.setDouble(f'/{device}/sines/1/phaseshift',phi)
                    inst.setDouble(f'/{device}/sigouts/{channels[1]}/offset',amp)
                    inst.sync()
                    power_data[i,j] = get_power(sa, inst, f_LO,threshold=threshold,plot=False,config=False)

        # find index of voltage corresponding to minimum LO leakage
        argmin = np.unravel_index(np.argmin(power_data), power_data.shape)

        opt_phi = phiArr[argmin[0]]
        opt_amp = ampArr[argmin[1]]
        # set voltages to optimal values
        inst.setDouble(f'/{device}/sines/1/phaseshift',opt_phi)
        inst.setDouble(f'/{device}/sigouts/{channels[1]}/offset',opt_amp)
        inst.sync()
        print(f'optimal phi = {round(opt_phi,3)}, optimal amp = {round(1e3*opt_amp,1)}')


        # find index of voltage corresponding to minimum LO leakage
        min_ind1 = np.argmin(imPower1)
        min_ind2 = np.argmin(imPower2)

        end = time.time()
        print('%s mixer Optimization took %.1f seconds'%(mixer,(end-start)))

        if plot:
            pf.plot_mixer_opt(phiArr, ampArr, power_data,cal='sb',element='qubit',fc=freq)
        # print('Ch1 Voltage (mV):%.2f\nCh2 Voltage (mV):%.2f\nOFF Power (dBm): %.1f\nON Power (dBm): %.1f'%(VoltRange1[min_ind1]*1e3,VoltRange2[min_ind2]*1e3,OFF_power,ON_power))


