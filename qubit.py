# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 15:43:05 2023

@author: lfl
"""

from scipy.signal.windows import gaussian
import csv
import glob
import os
from tqdm import tqdm
import time
import json
import numpy as np
from VISAdrivers.sa_api import *
import instrument_funcs as instfuncs
from zhinst.utils import create_api_session,convert_awg_waveform
from zhinst.toolkit import Session,CommandTable,Sequence,Waveforms
from zhinst.toolkit.waveform import Wave, OutputType
from scipy.interpolate import interp1d
import scipy as scy
import scipy.fftpack
import matplotlib.pyplot as plt
import seaborn as sns; sns.set() # styling
from scipy.signal import butter,lfilter,freqz,find_peaks,peak_widths
import re
import pandas as pd


pi=np.pi        
# plt.style.use(['science','no-latex'])
# plt.rcParams['figure.dpi'] = 150
# plt.rcParams['axes.facecolor'] = 'white'
# plt.rcParams['axes.edgecolor'] = 'black'
# plt.rcParams['axes.grid'] = False
# plt.rcParams['figure.frameon'] = True
# plt.rcParams.update(plt.rcParamsDefault)
# for a complete set of parameters "print(plt.rcParams)"
sns.set_style('ticks')
plt.rcParams['text.usetex'] = False
plt.rcParams['font.family'] =  'Arial'
plt.rcParams['font.size'] = 16
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["xtick.major.top"] = True
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.titlesize'] = 16
plt.rcParams["xtick.major.bottom"] = True
plt.rcParams["xtick.top"] = True
plt.rcParams["xtick.bottom"] = True
plt.rcParams["ytick.left"] = True
plt.rcParams["ytick.right"] = True
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams["ytick.major.right"] = True
plt.rcParams["ytick.labelright"] = False
plt.rcParams["ytick.minor.visible"] = False

device_name = "WM1"
project =  'dynamical-decoupling'

class qubit():
    
    #%% INITIALIZATION
    default_qb_pars = {
                    "qb_LO":                     3.829e9,
                    "qb_freq":                   3.879e9,
                    "qb_IF":                     0,
                    "qb_mixer_offsets":          [0,0], # I,Q
                    "qb_mixer_imbalance":        [0,0], # gain,phase
                    "pi_len":                       64, # needs to be multiple of 4
                    "pi_half_amp":                  0.2,
                    "pi_amp":                       0.45,
                    "gauss_len":                    64,
                    "gauss_amp":                    0.45,
                    "rr_LO":                        6.70438e9,
                    "rr_freq":                      6.4749e9,
                    'rr_IF':                        0,
                    "rr_mixer_offsets":             [0,0],
                    "rr_mixer_imbalance":           [0,0],
                    "amp_r":                        0.45,
                    'rr_atten':                     25,
                    "rr_resettime":                 0.5e-6,
                    'cav_resp_time':                0.5e-6,
                    'T1':                           50e-6,
                    'T2R':                          20e-6,
                    'T2E':                          30e-6,
                    "IQ_rotation":                  -0/180*np.pi, # phase rotation applied to IQ data
                    "analog_input_offsets":         [0,0],
                    }
    
    default_exp_pars = {
                    "Tmax":                         10e-6,
                    'active_reset':                 False,
                    "satur_dur":                    20e-6,
                    "tof":                          260, # time of flight in ns
                    }
    
    def __init__(self, qb):
        # load pars from json, OR create new json file
        self.name = qb
        
        try:
            print(f'Loading parameter {qb}_pars JSON file')
            with open(f'{qb}_pars.json', 'r') as openfile:
                self.qb_pars = json.load(openfile)
            # with open(f'{qb}_exp_pars.json', 'r') as openfile:
            #     self.exp_pars = json.load(openfile)

        except FileNotFoundError:
            print('Parameter file not found; loading parameters from template')
            self.qb_pars = self.default_qb_pars
            self.exp_pars = self.default_exp_pars

        self.inst_init()
        self.zi_init()

    def zi_init(self,qa_delay=900):
        '''
        Applies initial settings to zurich instrument devices

        qa_delay:        delay in time samples (1/1.8GHz) between trigger receive and integration
        '''
        awg_setting = [
            ['/dev8233/system/clocks/referenceclock/source', 1], # sets clock reference to external 10MHz
            ['/dev8233/system/awg/channelgrouping', 0], #sets grouping mode to 2x2
            ['/dev8233/sigouts/*/on', 1], #turn on outputs 1 and 2
            ['/dev8233/sigouts/*/range', 2], # sets range of awg
            ['/dev8233/system/awg/oscillatorcontrol', 1], #enables oscillator control via AWG (needed so we can resetOscPhase via sequence program)
            ['/dev8233/oscs/0/freq', self.qb_pars['qb_IF']], # sets the oscillator freq
            ['/dev8233/awgs/0/time', 1], # set AWG sampling rate to 1.2GHz
            ['/dev8233/sigouts/0/offset', self.qb_pars['qb_mixer_offsets'][0]],
            ['/dev8233/sigouts/1/offset', self.qb_pars['qb_mixer_offsets'][1]],
            ['/dev8233/awgs/0/outputs/0/modulation/mode', 3], # Output 1 modulated with Sine 1
            ['/dev8233/awgs/0/outputs/1/modulation/mode', 4], # Output 2 modulated with Sine 2
            ['/dev8233/sines/0/phaseshift', 0],   # Sine 1 phase
            ['/dev8233/sines/1/phaseshift' , 90],   # Sine 2 phase
            ['/dev8233/triggers/out/0/source', 4] # sets the marker 1 channel output
        ]
        print('Applying initial settings to HDAWG')
        self.IQ_imbalance(g=self.qb_pars['qb_mixer_imbalance'][0], phi=self.qb_pars['qb_mixer_imbalance'][1])
        self.awg.set(awg_setting)
        self.awg.sync()

        qa_setting = [
            # daq.setInt('/dev2528/awgs/0/outputs/0/mode', 1) # AWG set to modulation mode
            ['/dev2528/system/extclk', 1], # sets clock reference to external 10MHz
            ['/dev2528/sigouts/*/range',0.15], #sets output range to 150mV
            ['/dev2528/sigouts/*/on', 1], #turn on outputs 1 and 2
            ['/dev2528/sigins/0/range', 0.8],  #sets input range to 500mV
            ['/dev2528/sigouts/0/offset', self.qb_pars['rr_mixer_offsets'][0]],
            ['/dev2528/sigouts/1/offset', self.qb_pars['rr_mixer_offsets'][1]],
            ['/dev2528/system/jumbo', 1], # enables jumbo frames for faster connection | make sure you enable jumbo frames in windows network settings (how-to link:https://elements.tv/blog/achieve-maximum-ethernet-performance/)
            ['/dev2528/qas/0/integration/sources/0', 0], #sets channel mapping
            ['/dev2528/awgs/0/outputs/0/modulation/mode' , 0],
            ['/dev2528/qas/0/delay',round(qa_delay*1.8e9)], # sets delay between trigger receive and start of integration
            ['/dev2528/qas/0/integration/mode', 0], # 0 for standard (4096 samples max), 1 for spectroscopy
            ['/dev2528/oscs/0/freq', self.qb_pars['rr_IF']], # sets the oscillator freq
            ['/dev2528/qas/0/integration/length', 4096],
            ['/dev2528/qas/0/integration/trigger/channel', 7], # 0 for trigger input ch 1, 7 for internal triggering (QA AWG -> QA)]
            ['/dev2528/qas/0/monitor/trigger/channel', 7],
            ['/dev2528/qas/0/result/mode',0], #sets averaging to cyclic
            ['/dev2528/qas/0/result/reset', 1],
            ['/dev2528/awgs/0/auxtriggers/0/channel', 2], # sets the digital trigger signal 1 input channel to physical input channel 3
        ]
        print('Applying initial settings to UHFQA')
        self.qa.set(qa_setting)
        self.qa.sync()

    #%% setup_exp_funcs
    def setup_active_reset(self):
        """
        Sets up active reset. The discrimination threshold has to be previously established via a Rabi Measurement.

        Args:
            threshold (int, optional): The value that is going to be used to discriminate the measurement results. Defaults to 5.

        Returns:
            None.

        """
        # Configure AWG settings
        # select trigger sources
        self.awg.setInt('/dev8233/awgs/0/auxtriggers/0/channel', 0)
        self.awg.setInt('/dev8233/awgs/0/auxtriggers/1/channel', 1)
        # Select trigger slope. First trigger is QA Result Trigger (rise), second is QA Result (level)
        self.awg.setInt('/dev8233/awgs/0/auxtriggers/0/slope', 1)
        self.awg.setInt('/dev8233/awgs/0/auxtriggers/1/slope', 0)
        # sets trigger level
        self.awg.setDouble('/dev8233/triggers/in/0/level', 0.3)
        self.awg.setDouble('/dev8233/triggers/in/1/level', 0.3)

        #Configure QA settings
        # select trigger sources
        self.qa.setInt('/dev2528/triggers/out/0/source', 74) # QA result trigger. The QA sends a trigger to HDAWG Ch. 1 when measurement is done.
        self.qa.setInt('/dev2528/triggers/out/1/source', 64) # QA result. Sends trigger to HDAWG based on the measurement result.
        # set trigger mode to output ("drive")
        self.qa.setInt('/dev2528/triggers/out/0/drive', 1)
        self.qa.setInt('/dev2528/triggers/out/1/drive', 1)
        # set trigger levels to 3 V
        # self.qasetDouble('/dev2528/triggers/in/0/level', 3)
        # self.qasetDouble('/dev2528/triggers/in/1/level', 3)
        # sets QA result threshold
        self.qa.setDouble('/dev2528/qas/0/thresholds/0/level', self.exp_pars['threshold'])
        self.qa.sync()

    def pulsed_spec_setup(self):
        '''Setup HDAWG and UHFQA for pulsed spectroscopy'''
        self.awg.set('/dev8233/awgs/0/outputs/*/modulation/mode', 0)
        self.awg.setInt('/dev8233/awgs/0/time',2) # sets AWG sampling rate to 600 MHz
        self.setup_awg()
        

    def single_shot_setup(self):
        '''Setup HDAWG for single shot experiment'''
        print('-------------Setting HDAWG sequence-------------')
        self.awg.setInt('/dev8233/awgs/0/time',0) # sets AWG sampling rate to 2.4 GHz
        self.sequence = Sequence()
        self.sequence.code = self.single_shot_sequence()
        self.sequence.constants['qubit_reset_time'] = self.roundToBase(self.exp_pars['qubit_reset_time']*1.17e6)
        self.sequence.constants['num_samples'] = self.exp_pars['num_samples']
        N = self.qb_pars['pi_len']
        pi_pulse = self.qb_pars['pi_amp']*gaussian(N,N/5)
        self.sequence.waveforms = Waveforms()
        self.sequence.waveforms[1] = (Wave(pi_pulse, name="w_pi_I", output=OutputType.OUT1|OutputType.OUT2),
            Wave(np.zeros(N), name="w_pi_Q", output=OutputType.OUT1|OutputType.OUT2))
        self.upload_to_awg()
        

    def ssb_setup(self,qb_IF=50e6):

        self.awg.setDouble('/dev8233/oscs/0/freq', qb_IF)
        self.awg.setDouble('/dev8233/sines/1/phaseshift', 90)
        self.awg.setInt('/dev8233/awgs/0/outputs/0/modulation/mode', 3)
        self.awg.setInt('/dev8233/awgs/0/outputs/1/modulation/mode', 4)
        self.awg.setDouble('/dev8233/sigouts/0/range', 0.6)
        self.awg.setDouble('/dev8233/sigouts/1/range', 0.6)

        
    def single_shot(self,save_data=True):
        '''
        DESCRIPTION: Executes single shot experiment.

        '''
        result_length=2*self.exp_pars['num_samples']
        self.exp = 'single-shot'

        readout_pulse_length = 2.3e-6 + self.qb_pars['cav_resp_time'] + 2e-6

        
        self.single_shot_setup()
        # setup QA
        self.setup_qa_awg(readout_pulse_length) # setup QA AWG for readout
        time.sleep(0.1)

        sweep_data, paths = self.create_sweep_data_dict()
        data_pi = []
        data_OFF = []

       
        self.qa_result_reset()
        self.enable_awg(self.awg, enable=0)
        self.config_qa(result_length,source=7)
        self.qa.sync()

        self.qa_result_enable()
        self.enable_awg(self.qa) # start the readout sequence
        

        print('Start measurement')
        bt = time.time()
        self.enable_awg(self.awg,enable=1)
        data = self.acquisition_poll(paths, result_length, timeout = 5*self.exp_pars['num_samples']*self.exp_pars['qubit_reset_time']) # transfers data from the QA result to the API for this frequency point
        et = time.time()
        duration = et-bt
        print(f'Measurement time: %.1f s'%duration)
        # seperate OFF/ON data and average
        data_OFF = np.append(data_OFF, [data[paths[0]][k] for k in self.even(len(data[paths[0]]))])/4096
        data_pi =  np.append(data_pi, [data[paths[0]][k] for k in self.odd(len(data[paths[0]]))])/4096

        self.enable_awg(self.awg,enable=0)
        self.stop_result_unit( paths)
        self.enable_awg(self.qa, enable = 0)

        data = dict()
        data['I'] = np.real(data_OFF)
        data['Q'] = np.imag(data_OFF)
        data['Iexc'] = np.real(data_pi)
        data['Qexc'] = np.imag(data_pi)
        
        if save_data:
            self.save_data(project,device_name,data=np.vstack((data_OFF,data_pi)))
            
        return data

    def spectroscopy(self,element,freqs,save_data=True):
        '''
        DESCRIPTION: Executes qubit spectroscopy.

        '''
        if element == 'resonator':
            instfuncs.set_rr_output(False)
            
        self.exp = 'spectroscopy'
        self.exp_pars['fsAWG'] = 600e6
        if self.exp_pars['on_off']:
            result_length = 2*self.exp_pars['n_avg']
        else:
            result_length = self.exp_pars['n_avg']
            
        instfuncs.set_attenuator(self.exp_pars['rr_attenuation'])

        readout_pulse_length = 2.3e-6 + self.qb_pars['cav_resp_time'] + 2e-6

        self.pulsed_spec_setup()
        self.setup_qa_awg(readout_pulse_length)
        self.config_qa(result_length=result_length,delay=self.qb_pars['cav_resp_time'])
        
        print('Start measurement')
        sweep_data, paths = self.create_sweep_data_dict()
        data_ON = []
        data_OFF = []

        self.enable_awg(self.qa) # start the readout sequence
        
        bt = time.time()
        
        with tqdm(total = len(freqs)) as progress_bar:
            for f in freqs:
                instfuncs.set_qb_LO(f*1e9,sweep=True)
                self.qa_result_reset()
                self.qa_result_enable()
                self.enable_awg(self.awg) #runs the drive sequence
                data = self.acquisition_poll(paths, result_length, timeout = 10) # transfers data from the QA result to the API for this frequency point
                # seperate OFF/ON data and average
                if self.exp_pars['on_off']:
                    data_OFF = np.append(data_OFF, np.mean([data[paths[0]][k] for k in self.even(len(data[paths[0]]))]))
                else:
                    pass
                data_ON =  np.append(data_ON, np.mean([data[paths[0]][k] for k in self.odd(len(data[paths[0]]))]))
                progress_bar.update(1)

        et = time.time()
        duration = et-bt
        print(f'Measurement time: {duration:.1f} seconds')

        if self.exp_pars['on_off']:
            data = (data_ON-data_OFF)/4096
        else:
            data = data_ON/4096
        I_data= data.real
        Q_data = data.imag

        power_data = np.abs(I_data*I_data.conjugate()+Q_data*Q_data.conjugate())

        self.enable_awg(self.awg,enable=0)
        self.stop_result_unit(paths)
        self.enable_awg(self.qa, enable = 0)
        
        if save_data:
            self.save_data(project,device_name,data=np.vstack((freqs,I_data,Q_data)))
        
        if element == 'resonator':
            instfuncs.set_rr_output(True)
            
        return power_data,I_data,Q_data

    def pulsed_exp(self,exp='rabi',verbose=1,check_mixers=False,save_data=True):

        source = 2
        readout_pulse_length = 2.3e-6 + self.qb_pars['cav_resp_time'] + 0.5e-6

        self.update_qb_value('qb_IF', self.exp_pars['qubit_drive_freq']-self.qb_pars['qb_LO'])
        self.awg.setDouble('/dev8233/oscs/0/freq', self.qb_pars['qb_IF'])
        # if check_mixers:
        
        # stops AWGs and reset the QA
        self.enable_awg(self.awg,enable=0)
        self.enable_awg(self.qa,enable=0)
        self.qa_result_reset()
        
        self.exp = exp

        #setup HDAWG
        if self.exp == 'p-rabi':
            self.x0 = self.exp_pars['x0']
            self.xmax = self.exp_pars['xmax']
            self.dx = self.exp_pars['dx']
            x_array = np.arange(self.x0,self.xmax+self.dx/2,self.dx)
            self.n_steps = len(x_array)
        else:
            self.n_points,self.n_steps,self.x0,self.xmax,self.dx = self.calc_steps(verbose=verbose)
            
        
        if self.qb_pars['qb_IF'] != 0:
            self.awg.set('/dev8233/awgs/0/outputs/0/modulation/mode', 3)
            self.awg.set('/dev8233/awgs/0/outputs/1/modulation/mode', 4)
        else:
            self.awg.set('/dev8233/awgs/0/outputs/*/modulation/mode', 0)
            
        self.awg.setInt('/dev8233/awgs/0/time',int(2.4e9/self.exp_pars['fsAWG']-1)) # sets AWG sampling rate to 600 MHz
        self.setup_awg()
        
        # setup QA
        self.setup_qa_awg(readout_pulse_length=readout_pulse_length) # setup QA AWG for readout
        
        # setup QA Result unit
        self.config_qa(result_length=self.n_steps,source=source) # configure UHFQA result unit, source = 2 means data is rotated
        self.qa.sync()
        
        if self.exp_pars['active_reset'] == True:
            self.setup_active_reset(threshold=threshold)

        exp_dur = self.calc_timeout()
        print('Estimated Measurement Time (with/without active reset): %.3f/%.3f sec'%(int(1/8*exp_dur),exp_dur))

        
        if self.exp_pars['active_reset'] == True:
            timeout = 0.2*1.2*exp_dur
        else:
            timeout = 1.2*exp_dur

        sweep_data, paths = self.create_sweep_data_dict()

        self.enable_awg(self.qa,enable=1) # start the readout sequence

        self.qa_result_enable() # arm the qa

        str_meas = time.time()
        self.enable_awg(self.awg,enable=1) #runs the drive sequence
        data = self.acquisition_poll(paths, num_samples = self.n_steps, timeout = 2*timeout) # retrieve data from UHFQA

        for path, samples in data.items():
            sweep_data[path] = np.append(sweep_data[path], samples)

        # reset QA result unit and stop AWGs
        self.stop_result_unit(paths)
        self.enable_awg(self.awg, enable = 0)
        self.enable_awg(self.qa, enable = 0)


        end_meas = time.time()
        print('\nmeasurement duration: %.1f s' %(end_meas-str_meas))
        
        x_array = self.get_xdata_frm_ct()
        if self.exp != 'p-rabi':
            x_array = x_array/self.exp_pars['fsAWG']
        elif self.exp == 'p-rabi':
            x_array = self.get_xdata_frm_ct()
            
        data = sweep_data[paths[0]][0:self.n_steps]/4096
        
        if source == 2 or source == 1:
            results = data
        elif source == 7:
            I = data.real
            Q = data.imag
            results = [[I],[Q]]

        if save_data:
            self.save_data(project,device_name,data=np.vstack((x_array,results)))
            
        return x_array[0],results,self.n_steps
        
        # self.create_wfms(sequence=sequence, n_points=n_points, Tmax=Tmax)


        # elif setup[0] == 2:
        #     bt = time.time()
        #     # replace waveforms, don't recompile program
        #     n_points,n_steps,pulse_length_increment,pulse_length_start = self.calc_n_steps(sequence=sequence,fsAWG=fs,piWidth_Y=piWidth_Y,
        #                                                                            stepSize=stepSize,Tmax=Tmax,verbose=verbose)
        #     if B0 == 0:
        #         noise_instance = np.zeros(n_points)
        #     if mu != 0 or sigma != 0:
        #         white_noise = np.random.normal(mu, sigma, n_points)
        #         waveforms_native = convert_awg_waveform(wave1=noise_instance,wave2=white_noise)
        #     else:
        #         waveforms_native = convert_awg_waveform(wave1=noise_instance)
        #     path = '/dev8233/awgs/0/waveform/waves/0'
        #     self.awg.setVector(path,waveforms_native)
        #     self.awg.sync()
        #     et = time.time()
        #     print('Replacing waveforms took: %.1f ms'%(1e3*(et-bt)))

    #%% setup_HDAWG
    def setup_awg(self):
        
        # Initialize sequence class
        self.sequence = Sequence()
        # setup sequence code
        self.sequence.code = self.gen_seq_code()
        # setup constants
        self.setup_seq_pars()
        # setup waveforms
        self.setup_waveforms()
        if self.exp == 'p-rabi':
            self.make_ct(sweep_var='amp')
        elif self.exp != 'spectroscopy': 
            # setup command table    
            self.make_ct(sweep_var='time')
        
        # upload everything to awg
        self.upload_to_awg()
            
    def gen_seq_code(self):
        
        if self.exp == 'spectroscopy':
            code = self.spec_sequence()
        elif self.exp == 'rabi':
            code = self.rabi_sequence()
        elif self.exp == 'p-rabi':
            code = self.power_rabi_sequence()
        elif self.exp == 'T1':
            code = self.T1_sequence()
        elif self.exp == 'ramsey':
            code = self.ramsey_sequence()
        elif self.exp == 'echo':
            code = self.echo_sequence()

            
        
        return code    
    
    def setup_seq_pars(self):
        
        # i = 0
        # for key,value in self.exp_pars.items():
        #     if (i > 1 and self.exp == 'spectroscopy') or i > 2:
        #         break
        #     else:
        #         if key == 'qubit_reset_time':
        #             value = self.roundToBase(value*self.exp_pars['fsAWG']) # converts the reset time from us to num of samples
        #         else:
        #             pass
        #         self.sequence.constants[key] = value
        #     i += 1
        self.sequence.constants['n_avg'] = self.exp_pars['n_avg']
        self.sequence.constants['qubit_reset_time'] = self.roundToBase(self.exp_pars['qubit_reset_time']*1.17e6)
        if self.exp != 'spectroscopy':
            self.sequence.constants['n_steps'] = self.n_steps
     
    def setup_waveforms(self):
        
        self.sequence.waveforms = Waveforms()
        
        if self.exp == 'spectroscopy':
            
            qubit_drive_dur = self.roundToBase(self.exp_pars['satur_dur']*self.exp_pars['fsAWG'])
            const_pulse = 2*self.exp_pars['amp_q'] * np.ones(qubit_drive_dur)
            self.sequence.waveforms[0] = (Wave(const_pulse, name="w_const", output=OutputType.OUT1),
                Wave(np.zeros(qubit_drive_dur), name="w_zero", output=OutputType.OUT2))
            
            
        elif self.exp == 'rabi':
            
            N = 64
            # gauss_pulse = self.calc_amp('dev8233', amp=self.exp_pars['amp_q'])*gaussian(N,N/5)
            gauss_pulse = self.exp_pars['amp_q']*gaussian(N,N/5)
            gauss_rise = gauss_pulse[:int(N/2)]
            gauss_fall = gauss_pulse[int(N/2):]
            
            self.sequence.waveforms[0] = (Wave(gauss_rise, name="w_gauss_rise_I", output=OutputType.OUT1|OutputType.OUT2),
                Wave(np.zeros(len(gauss_rise)), name="w_zero1", output=OutputType.OUT1|OutputType.OUT2))
                
            self.sequence.waveforms[1] = (Wave(gauss_fall, name="w_gauss_fall_I", output=OutputType.OUT1|OutputType.OUT2),
                Wave(np.zeros(len(gauss_fall)), name="w_zero2", output=OutputType.OUT1|OutputType.OUT2))
            
            self.sequence.waveforms[2] = (Wave(self.exp_pars['amp_q']*np.ones(self.n_points), name="w_const_I", output=OutputType.OUT1|OutputType.OUT2),
                Wave(np.zeros(self.n_points), name="w_const_Q", output=OutputType.OUT1|OutputType.OUT2))
            
        elif self.exp == 'p-rabi':
            N = self.qb_pars['gauss_len']
            gauss_pulse = gaussian(N,N/5)
            self.sequence.waveforms[0] = (Wave(gauss_pulse, name="w_gauss", output=OutputType.OUT1|OutputType.OUT2),
                Wave(np.zeros(N), name="w_zero", output=OutputType.OUT1|OutputType.OUT2))
            
        elif self.exp == 'T1':
            
            N = self.qb_pars['pi_len']
            pi_pulse = self.calc_amp('dev8233', amp=self.qb_pars['pi_amp'])*gaussian(N,N/5)
            
            self.sequence.waveforms[0] = (Wave(pi_pulse, name="w_pi_I", output=OutputType.OUT1|OutputType.OUT2),
                Wave(np.zeros(N), name="w_pi_Q", output=OutputType.OUT1|OutputType.OUT2))
            self.sequence.waveforms[1] = (Wave(np.zeros(self.n_points), name="w_zero_I", output=OutputType.OUT1|OutputType.OUT2),
                Wave(np.zeros(self.n_points), name="w_zero_Q", output=OutputType.OUT1|OutputType.OUT2))
                
        elif self.exp == 'ramsey':
            N = self.qb_pars['pi_len']
            pi2_pulse = self.calc_amp('dev8233', amp=self.qb_pars['pi_half_amp']) * gaussian(N,N/5)
            self.sequence.waveforms[0] = (Wave(pi2_pulse, name="w_pi2_I", output=OutputType.OUT1|OutputType.OUT2),
                Wave(np.zeros(N), name="w_pi2_Q", output=OutputType.OUT1|OutputType.OUT2))
            self.sequence.waveforms[1] = (Wave(np.zeros(self.n_points), name="w_zero_I", output=OutputType.OUT1|OutputType.OUT2),
                Wave(np.zeros(self.n_points), name="w_zero_Q", output=OutputType.OUT1|OutputType.OUT2))
            
        elif self.exp == 'echo':
            N = self.qb_pars['pi_len']
            pi2_pulse = self.calc_amp('dev8233', amp=self.qb_pars['pi_half_amp']) * gaussian(N,N/5)
            pi_pulse = self.calc_amp('dev8233', amp=self.qb_pars['pi_amp'])*gaussian(N,N/5)
            
            self.sequence.waveforms[0] = (Wave(pi2_pulse, name="w_pi2_I", output=OutputType.OUT1|OutputType.OUT2),
                Wave(np.zeros(N), name="w_pi2_Q", output=OutputType.OUT1|OutputType.OUT2))
            self.sequence.waveforms[1] = (Wave(pi_pulse, name="w_pi_I", output=OutputType.OUT1|OutputType.OUT2),
                Wave(np.zeros(N), name="w_pi_Q", output=OutputType.OUT1|OutputType.OUT2))
            self.sequence.waveforms[2] = (Wave(np.zeros(self.n_points), name="w_zero_I", output=OutputType.OUT1|OutputType.OUT2),
                Wave(np.zeros(self.n_points), name="w_zero_Q", output=OutputType.OUT1|OutputType.OUT2))
            
        elif self.exp == 'single_shot' or self.exp_pars['active_reset'] == True:
            N = self.qb_pars['pi_len']
            pi_pulse = self.qb_pars['pi_amp']*gaussian(N,N/5)
            self.sequence.waveforms[1] = (Wave(pi_pulse, name="w_pi_I", output=OutputType.OUT1|OutputType.OUT2),
                Wave(np.zeros(N), name="w_pi_Q", output=OutputType.OUT1|OutputType.OUT2))
            
    def setup_mixer_calib(self,inst,amp = 1):
        self.awg.setInt('/dev8233/awgs/0/time',11) # sets AWG sampling rate to 600 MHz
        self.sequence = Sequence()
        self.sequence.code = self.mixer_calib_sequence()
        self.sequence.constants['amp'] = amp
        self.sequence.constants['N'] = self.roundToBase(0.1*1.17e6)
        if inst == self.awg:
            with self.hdawg_core.set_transaction():
                    self.hdawg_core.awgs[0].load_sequencer_program(self.sequence)
        elif inst == self.qa:
            with self.hdawg_core.set_transaction():
                self.qa_awg_core.awgs[0].load_sequencer_program(self.sequence)
            
        
    def spec_sequence(self):

        
        if self.exp_pars['on_off']:
            awg_program = '''       
    
            // Beginning of the core sequencer program executed on the HDAWG at run time
            wave w_marker = 2*marker(512,1);
    
            repeat(n_avg) {
                // OFF Measurement
                playWave(1,w_marker);
                playZero(qubit_reset_time,AWG_RATE_1P2MHZ);
                // ON measurement
                playWave(1,w_const,2,w_zero);
                playWave(1,w_marker);
                playZero(qubit_reset_time,AWG_RATE_1P2MHZ);
                        }'''
        else:
            awg_program = '''       
    
            // Beginning of the core sequencer program executed on the HDAWG at run time
            wave w_marker = 2*marker(512,1);
    
            repeat(n_avg) {
                playWave(1,w_const,2,w_zero);
                playWave(1,w_marker);
                playZero(qubit_reset_time,AWG_RATE_1P2MHZ);
                        }'''
            
        return awg_program

    def rabi_sequence(self):
        
        awg_program = '''
        resetOscPhase();
        
        wave w_marker = 2*marker(512,1);
        var i;
        // Beginning of the core sequencer program executed on the HDAWG at run time
        repeat(n_avg){
            for (i=0; i<n_steps; i++) {
                    playWave(1,2,w_gauss_rise_I,1,2,w_zero1);
                    executeTableEntry(i);
                    playWave(1,2,w_gauss_fall_I,1,2,w_zero2);
                    playWave(1,w_marker);
                    playZero(qubit_reset_time,AWG_RATE_1P2MHZ);
          }
        }'''

        return awg_program
    
    def power_rabi_sequence(self):
        
        awg_program = '''
        resetOscPhase();
        
        wave w_marker = 2*marker(512,1);
        var i;
        // Beginning of the core sequencer program executed on the HDAWG at run time
        repeat(n_avg){
            for (i=0; i<n_steps; i++) {
                    executeTableEntry(i);
                    playWave(1,w_marker);
                    playZero(qubit_reset_time,AWG_RATE_1P2MHZ);
          }
        }'''

        return awg_program
    
    def T1_sequence(self,active_reset=0):
        awg_program = '''
        resetOscPhase();
        wave w_marker = 2*marker(512,1);
        var i;
        // Beginning of the core sequencer program executed on the HDAWG at run time
        repeat(n_avg){
            for (i=0; i<n_steps; i++) {
                    playWave(1,2,w_pi_I,1,2,w_pi_Q);
                    executeTableEntry(i);
                    playWave(1,w_marker);
                    playZero(qubit_reset_time,AWG_RATE_1P2MHZ);
          }
        }'''

        return awg_program
    
    def ramsey_sequence(self,nAverages=1000,n_steps=100,active_reset=0):
        
        awg_program = '''
            resetOscPhase();
            
            wave w_marker = 2*marker(512,1);
            var i;
            // Beginning of the core sequencer program executed on the HDAWG at run time
            repeat(n_avg){
                for (i=0; i<n_steps; i++) {
                        playWave(1,2,w_pi2_I,1,2,w_pi2_Q);
                        executeTableEntry(i);
                        playWave(1,2,w_pi2_I,1,2,w_pi2_Q);
                        playWave(1,w_marker);
                        playZero(qubit_reset_time,AWG_RATE_1P2MHZ);
              }
            }'''

        return awg_program

    def echo_sequence(self,nAverages=1000,n_steps=100,active_reset=0):

        awg_program = '''
        resetOscPhase();
        wave w_marker = 2*marker(512,1);
        var i;
        // Beginning of the core sequencer program executed on the HDAWG at run time
        repeat(n_avg){
            for (i=0; i<n_steps; i++) {
                    playWave(1,2,w_pi2_I,1,2,w_pi2_Q);
                    executeTableEntry(i);
                    playWave(1,2,w_pi_I,1,2,w_pi_Q);
                    executeTableEntry(i);
                    playWave(1,2,w_pi2_I,1,2,w_pi2_Q);
                    playWave(1,w_marker);
                    playZero(qubit_reset_time,AWG_RATE_1P2MHZ);
          }
        }'''

        return awg_program
    
    def single_shot_sequence(self):
        awg_program = '''
        resetOscPhase();
        
        wave w_marker = 2*marker(512,1);
        // Beginning of the core sequencer program executed on the HDAWG at run time
        repeat(num_samples) {
            // OFF Measurement
            playWave(1,w_marker);
            playZero(qubit_reset_time,AWG_RATE_1P2MHZ);
            playWave(1,2,w_pi_I,1,2,w_pi_Q);
            playWave(1,w_marker);
            playZero(qubit_reset_time,AWG_RATE_1P2MHZ);
            }'''

        return awg_program
    
    def mixer_calib_sequence(self):
        
        awg_program = '''
            resetOscPhase();
            wave w_const = amp*ones(N);
            wave w_zeros = zeros(N);

            while (true) {
                playWave(1,2,w_const,1,2,w_zeros);
                waitWave();
                }
        '''

        return awg_program
    
    def reset_sequence(self):

        awg_program = '''
            waitDigTrigger(1);
            wait(1);
            wait(2000);
            if (getDigTrigger(2) == 0) {
                wait(10);
            } else {
                playWave(1,2,w_pi_I,1,2,w_pi_Q);
                }
            wait(10000);
          '''
          
        return awg_program
         
    #%% setup_QA
    def setup_qa_awg(self,readout_pulse_length=5e-6):
        
        self.qa.setInt('/dev2528/awgs/0/time',2)
        
        self.qa_sequence = Sequence()
        self.qa_sequence.code = """
        const cav_resp_time = delay;

        // Readout pulse
        wave readoutPulse = 0.67*ones(L);

        repeat(N) {
            waitDigTrigger(1,1);
            startQA();
            playWave(1,readoutPulse);
        }
        """
        
        self.qa_sequence.constants['delay'] = self.roundToBase(self.qb_pars['cav_resp_time']*450e6)
        self.qa_sequence.constants['L'] = self.roundToBase(readout_pulse_length*450e6)
        self.qa_sequence.constants['N'] = self.exp_pars['n_avg']*self.n_steps 
        # self.qa_sequence.constants['amp'] = self.calc_amp('dev2528', amp=0.1)
        
        with self.qa_awg_core.set_transaction():
            self.qa_awg_core.awgs[0].load_sequencer_program(self.qa_sequence)
            
        self.qa.setInt('/dev2528/awgs/0/auxtriggers/0/channel', 2) # sets the source of digital trigger 1 to be the signal at trigger input 3 (back panel)
        self.qa.setDouble('/dev2528/triggers/in/2/level', 0.1)
        self.qa.setInt('/dev2528/awgs/0/auxtriggers/0/slope', 1)
        

    def config_qa(self,result_length=1,delay=300e-9,source=7):
        # print('-------------Configuring QA-------------\n')
        base_rate=1.8e9
        bypass_crosstalk=0
        # set modulation frequency of QA AWG to some IF and adjust input range for better resolution
        self.qa.setDouble('/dev2528/sigins/0/range',0.5)
        # self.qa.setInt('/dev2528/oscs/0/freq'.format(device),int(rr_IF))
        # self.qa.setInt('/dev2528/awgs/0/outputs/0/mode', 1) # AWG set to modulation mode

        # QA setup Settings
        self.qa.setInt('/dev2528/qas/0/integration/sources/0', 0)
        self.qa.setInt('/dev2528/qas/0/delay',round(delay*base_rate))
        # if sequence =='spec':
        #     self.qa.setInt('/dev2528/qas/0/integration/mode'.format(device), 0) # 0 for standard (4096 samples max), 1 for spectroscopy
        # elif sequence =='pulse':
        #     self.qa.setInt('/dev2528/qas/0/integration/mode'.format(device), 0)
        self.qa.setInt('/dev2528/qas/0/integration/mode', 0) # 0 for standard (4096 samples max), 1 for spectroscopy
        self.qa.setInt('/dev2528/qas/0/integration/length', 4096)
        self.qa.setInt('/dev2528/qas/0/bypass/crosstalk', bypass_crosstalk)   #No crosstalk matrix
        self.qa.setInt('/dev2528/qas/0/bypass/deskew', 1)   #No crosstalk matrix
        # self.qa.setInt('/dev2528/qas/0/bypass/rotation'.format(device), 1)   #No rotation
        # x = np.linspace(0,integration_length,round(integration_length*base_rate))
        # weights_I = np.sin(2*np.pi*rr_IF*x)
        # weights_Q = np.zeros(round(integration_length*base_rate))
        weights_I = weights_Q = np.ones(4096)
        # weights_Q = np.zeros(round(integration_length*base_rate))
        self.qa.setVector('/dev2528/qas/0/integration/weights/0/real', weights_I)
        self.qa.setVector('/dev2528/qas/0/integration/weights/0/imag', weights_Q)
        self.qa.sync()
        self.qa.setInt('/dev2528/qas/0/integration/trigger/channel', 7); # 0 for trigger input ch 1, 7 for internal triggering (QA AWG -> QA)

        # QA Monitor Settings
        self.qa.setInt('/dev2528/qas/0/monitor/trigger/channel', 7)
        if self.exp == 'spectroscopy':
            self.qa.setInt('/dev2528/qas/0/monitor/averages',1)
            self.qa.setInt('/dev2528/qas/0/result/averages', 1)
        else:
            self.qa.setInt('/dev2528/qas/0/monitor/averages',self.exp_pars['n_avg'])
            self.qa.setInt('/dev2528/qas/0/result/averages', self.exp_pars['n_avg'])
        self.qa.setInt('/dev2528/qas/0/monitor/length', 4096)
        # configure triggering (0=trigger input 1 7 for internal trigger)

        # QA Result Settings
        self.qa.setInt('/dev2528/qas/0/result/length', result_length)
        
        self.qa.setInt('/dev2528/qas/0/result/source', source) # 2 -> source = rotation | 7 = integration
        self.qa.setInt('/dev2528/qas/0/result/reset', 1)
        self.qa.setInt('/dev2528/qas/0/result/enable', 1)
        self.qa.setInt('/dev2528/qas/0/result/mode',0) # cyclic averaging
        self.qa.sync()
        
    def qa_result_reset(self,reset = 1):
        '''
        result reset of QA

        '''

        self.qa.setInt('/dev2528/qas/0/result/reset', reset)

    def qa_result_enable(self,enable = 1):
        '''
        enable QA result
        '''

        self.qa.setInt('/dev2528/qas/0/result/enable', enable)

    def create_sweep_data_dict(self):
        '''
        create sweep data dictionary for ch1 and ch2

        device:         device ID
        '''

        channels = [0, 1]
        # Subscribe to result waves
        paths = []
        for ch in channels:
            path = f'/dev2528/qas/0/result/data/{ch}/wave'
            paths.append(path)
        self.qa.subscribe(paths)
        sweep_data = dict()

        for path in paths:
            sweep_data[path] = np.array([])
        return sweep_data, paths

    def acquisition_poll(self,paths, num_samples, timeout):
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

            dataset = self.qa.poll(poll_length, poll_timeout, poll_flags, poll_return_flat_dict)

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
    
    def stop_result_unit(self, paths):
        '''
        stop QA result unit,

        daq:             dag ID
        device:          device ID
        paths:           data paths
        '''
        self.qa.unsubscribe(paths)
        self.qa.setInt('/dev2528/qas/0/result/enable', 0)
    #%% signal_funcs
    def pull_wfm(self,sweep_name,nu,tauk,sequence='ramsey'):
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

    def create_wfms(self,sequence="ramsey",n_points=1000,Tmax=5e-6):
        '''Generates all waveform text files to be used by the HDAWG sequence file'''
        # # create RTN noise or pull instance from file (only for parameter sweeps)
        # if B0 != 0 and len(noise_instance) == 0:
        #     print('Generating New Waveform')
        #     t = np.linspace(0,Tmax,n_points)
        #     if wk == 0:
        #         qubit_free_evol = B0 * np.cos(2*np.pi*nu*1e3*t+phi*2*np.pi*np.random.rand()) * gen_tel_noise(n_points, tauk, dt=Tmax/n_points)
        #     else:
        #         qubit_free_evol = B0 * gen_WK_sig(fs=2*n_points/Tmax, nu=nu*1e3, tauk=tauk*1e-6, Tmax=Tmax)
        # elif B0 != 0 and len(noise_instance) > 0:
        #     print('Loading Waveform from csv File')
        #     qubit_free_evol = noise_instance
        # else:
        #     qubit_free_evol = np.zeros(n_points)

        # qubit_free_evol = qubit_free_evol[...,None] # transpose waveforms such that they are columns (necessary such that they are readable by AWG seqc files)

        # # create white noise instance
        # if mu != 0 or sigma != 0:
        #     if sweep == 0:
        #         white_noise = np.random.normal(loc=mu, scale=sigma, size=n_points)
        #     elif sweep == 1:
        #         white_noise = white_noise_instance
        
        #     white_noise = white_noise[...,None]
        #     wfm_arr = np.hstack((qubit_free_evol,white_noise)) # stack waveform columns horizontally
        # else:
        #     wfm_arr = qubit_free_evol
        if sequence == 'rabi':
            fsAWG = 2.4e9
            gauss_pulse = gaussian(self.roundToBase(self.qb_pars['gauss_len']*self.exp_pars['fsAWG']), self.roundToBase(self.pars['gauss_len']/5*self.exp_pars['fsAWG']))
            drive_pulse = np.concatenate((gauss_pulse[:len(gauss_pulse)/2],np.ones(self.roundToBase(Tmax*fsAWG)),gauss_pulse[len(gauss_pulse)/2+1:]))
            # save wfm to file
            self.make_wfm_file(filename='drive_pulse', wfm_data=drive_pulse)
        # elif sequence == 'T1':
            
        # elif sequence == 'ramsey':
            

    def gen_WK_sig(self,fs,nu,tauk,Tmax):

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

    def calc_autocorr(self,sig):
        '''Calculates the autocorrelation of the given signal'''
        return sm.tsa.acf(sig,nlags=len(sig))

    def gen_noise_realizations(self,par1_arr=np.linspace(0,10,100),par2_arr=[0],numRealizations=3,n_points=1000,T_max=5e-6,sweep_count=1,
                               meas_device='CandleQubit_6',sequence='ramsey',wk=False,plot=False):
        """
        Generates noise realizations and saves them to a csv file for parameter sweep

        Args:
            par1_arr (array, optional): array of tauk values. Defaults to np.linspace(0,10,100).
            par2_arr (array, optional): arary of nu values. Defaults to [0].
            numRealizations (TYPE, optional): DESCRIPTION. Defaults to 3.
            n_points (int, optional): DESCRIPTION. Defaults to 1000.
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
        t = np.linspace(0,T_max,n_points)
        parent_dir = 'E:\\generalized-markovian-noise\\%s\\sweep_data\\%s\\'%(meas_device,sequence)
        directory = 'sweep_%03d\\noise_instances'%(sweep_count)
        path = os.path.join(parent_dir,directory)
        noise_arr = np.zeros((numRealizations,n_points))
        for i in range(numPoints_par2):
            for k in range(numPoints_par1):
                filename = "nu_%d_Hz_tau_%d_ns.csv" % (round(par2_arr[i]*1e3),round(par1_arr[k]*1e3))
                with open(os.path.join(path,filename),"w",newline="") as datafile:
                    writer = csv.writer(datafile)
                    for j in range(numRealizations):
                        if len(par2_arr) > 1 or par2_arr != 0:
                            if wk:
                                noise_arr[j,:] = np.cos(2*np.pi*par2_arr[i]*1e3*t + phi*2*np.pi*np.random.rand()) * gen_tel_noise(n_points, par1_arr[k], dt = T_max/n_points)
                            elif wk:
                                noise,psd,freqs,autocorr = gen_WK_sig(fs=2*n_points/T_max, nu=par2_arr[i]*1e3, tauk=par1_arr[k]*1e-6, Tmax=1e-3)
                                noise_arr[j,:] = noise[:n_points]/max(noise[:n_points])
                        elif len(par2_arr) <= 1 and par2_arr[0] == 0:
                            noise_arr[j,:] = gen_tel_noise(n_points, par1_arr[k], dt = T_max/n_points)

                    writer.writerows(noise_arr)

        if plot:
            fig = plt.figure(figsize=(4,8),dpi=300)
            ax1 = fig.add_subplot(3,1,1) # noise realization plot
            ax2 = fig.add_subplot(3,1,2) # mean autocorrelation plot
            ax3 = fig.add_subplot(3,1,3) # PSD plot (real & imag)

            ac = np.zeros(n_points)
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

    def gen_tel_noise(self,numPoints,tau,dt):
        '''Generates a single instance of telegraph noise'''
        signal = np.ones(numPoints)*(-1)**np.random.randint(0,2)
        for i in range(1,numPoints-1):
            if np.random.rand() < 1/(2*tau*1e-6/dt)*np.exp(-1/(2*tau*1e-6/dt)):
                signal[i+1] = - signal[i]
            else:
                signal[i+1] = signal[i]
        return signal

    def qa_monitor_avg(self,length,averages):

        settings = [
            ("qas/0/monitor/enable", 0),
            ("qas/0/monitor/length", length),
            ("qas/0/monitor/averages", averages),
            ("qas/0/monitor/enable", 1),
            ("qas/0/monitor/reset", 1),
            ('qas/0/monitor/trigger/channel', 7)
        ]
        self.qa.set([(f"/{'dev2528'}/{node}", value) for node, value in settings])
        # Signals to measure

        paths = []
        for channel in range(2):
            path = f"/{'dev2528':s}/qas/0/monitor/inputs/{channel:d}/wave"
            paths.append(path)
        self.qa.setInt('/dev2528/qas/0/monitor/reset', 1)
        self.qa.setInt('/dev2528/qas/0/monitor/enable', 1)

        self.qa.sync()
        time.sleep(0.2)
        # self.qa.setInt('dev2528/qas/0/monitor/trigger/channel',1)
        self.qa.subscribe(paths)

        # Perform acquisition
        print("Acquiring data...")
        data = self.acquisition_poll(paths, length,timeout=60)
        self.qa.setInt('/dev2528/qas/0/monitor/enable', 0)
        fig = plt.figure(figsize=(12,6))
        plt.plot(data[paths[0]])
        plt.plot(data[paths[1]])
        print(len(data[paths[0]]))
        print(len(data[paths[1]]))
        plt.title(f'Input signals after {averages:d} averages')
        plt.xlabel('n_points')
        plt.ylabel('Amp (V)')
        plt.grid()
        plt.show()

        return data

    def scope_meas(length=8192,nAverages=128,samp_rate=1.8e9,trigLevel=0.1):

        '''Executes a measurement with the UHFQA'''
        #setup and initialize scope
        scope = self.config_scope(scope_length=length,scope_avg_weight=1,scope_mode=0)
        scope = self.qa.scopeModule()
        self.qa.setInt('/dev2528/scopes/0/channel', 3)# 1: ch1; 2(DIG):ch2; 3: ch1 and ch2(DIG)
        self.qa.setInt('/dev2528/scopes/0/channels/0/inputselect', 0) # 0: sigin 1; 1: sigin 2
        self.qa.setDouble('/dev2528/scopes/0/length', length)
        scope.set('scopeModule/averager/weight', 1)
        scope.set('scopeModule/mode', 2)
        self.qa.setDouble('dev2528/scopes/0/length', length*nAverages)
        self.qa.setInt('/dev2528/scopes/0/single', 1)
        self.qa.setInt('/dev2528/scopes/0/trigchannel', 0)
        self.qa.setInt('/dev2528/scopes/0/trigenable', 1)
        self.qa.setDouble('/dev2528/scopes/0/trigholdoff', 50e-6) #In units of second. Defines the time before the trigger is rearmed after a recording event
        self.qa.setDouble('/dev2528/scopes/0/triglevel', trigLevel) #in Volts
        # self.qa.setInt('/dev2528/scopes/0/segments/enable', 1)
        self.qa.setInt('/dev2528/scopes/0/time', int(np.log2(1.8e9/samp_rate))) #set sampling rate
        self.qa.setInt('/dev2528/scopes/0/channel', 3) # enables both signal inputs
        # self.qa.setInt('/dev2528/scopes/0/segments/count',nAverages)
        self.qa.setDouble('/dev2528/sigins/0/range',0.4)
        self.qa.sync()

        self.restart_avg_scope(scope)
        self.enable_scope(enable=1)
        self.subscrib_scope(scope,'dev2528')
        self.execute_scope(scope)

        self.enable_awg(self.awg,enable=1)

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

        self.enable_scope(enable=0)
    #    scope.set('scopeModule/clearhistory', 0)
        self.finish_scope(scope)

        # return scope,avgCh1Data,avgCh2Data
        return scope,ch1Data,ch2Data

    def calc_timeout(self):
        '''Calculates timeout for experiment. If measurement takes more than timeout seconds, then the program stops and gives an error'''
        t = (self.n_steps*(self.dx/self.exp_pars['fsAWG']+self.exp_pars['qubit_reset_time']))*self.exp_pars['n_avg']
        return t

    def init_arrays(numRealizations=128,interval=2,n_pointsBackground=200,n_points=200):
        '''Initializes arrays to be used for storing exp data during parameter sweep'''
        bData = np.zeros((int(numRealizations/interval),n_pointsBackground),dtype=float)
        data = np.zeros((numRealizations,n_points),dtype=float)
        return bData,data

    # def set_AWG_output_amplitude(range)

    def calc_sweep_time(par1,par2,measTimeBackground=1,measTime=25,nMeasBackground=100,nMeas=100):
        return (measTimeBackground*nMeasBackground+measTime*nMeas)*len(par1)*len(par2)


    def create_echo_wfms(fs=1.2e9,mu=0,sigma=0,B0=0,n_points=1024,Tmax=5e-6,amp=0,pi2Width=50e-9,n_steps=101,pulse_length_increment=32):

        '''
        DESCRIPTION: Generates a series of waveforms to be uploaded into the AWG. The output is a series of csv files.
        '''

        start = time.time()
        ACpre = mu*np.ones(roundToBase(1500e-9*fs))
        pi2 = amp*np.ones(int(fs*pi2Width))
        pi2pre = 0 * ACpre
        ac_noise = np.random.normal(mu, sigma, 2*n_points)
        tel_noise = np.zeros(2*n_points)
        for i in range(n_steps):
            ch1_wfm = np.concatenate((pi2pre,pi2,tel_noise[0:i*pulse_length_increment],pi2,pi2,tel_noise[i*pulse_length_increment:2*i*pulse_length_increment],pi2,pi2pre))
            ch2_wfm = np.concatenate((ACpre,mu*pi2/amp,ac_noise[0:i*pulse_length_increment],mu*pi2/amp,mu*pi2/amp,ac_noise[i*pulse_length_increment:2*i*pulse_length_increment],mu*pi2/amp,ACpre))
            ch1_wfm = ch1_wfm[...,None]
            ch2_wfm = ch2_wfm[...,None]
            wfm_2D_arr = np.hstack((ch1_wfm,ch2_wfm))
            np.savetxt("C:/Users/LFL/Documents/Zurich Instruments/LabOne/WebServer/awg/waves/"+"echo_"+"wfm_%03d"%(i)+".csv", wfm_2D_arr, delimiter = ",")

        end = time.time()
        print('Generating echo Waveforms took %.1f' %(end-start))

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

    #%% command_table_funcs
    # Please refer to https://docs.zhinst.com/hdawg/commandtable/v2/schema for other settings
    def make_ct(self,sweep_var='time'):
        
        self.init_ct()
        
        if self.exp == 'rabi':
            wfm_index = 2
        elif self.exp == 'p-rabi':
            wfm_index = 0
        elif self.exp == 'T1':
            wfm_index = 1
        elif self.exp == 'ramsey':
            wfm_index = 1
        if self.exp == 'echo':
            wfm_index = 2
            
        if sweep_var == 'time':
            self.ct_sweep_length(wfm_index)
        elif sweep_var == 'amp':
            self.ct_sweep_amp(wfm_index)
        
    def ct_sweep_length(self,wfm_index):
            
        for i in range(self.n_steps):
            wfm_length = self.x0 + i * self.dx
            self.ct.table[i].waveform.index = wfm_index
            self.ct.table[i].waveform.length = wfm_length
            
    def ct_sweep_amp(self,wfm_index):
        
        for i in range(self.n_steps):
            amp = self.exp_pars['x0'] + i * self.exp_pars['dx']
            self.ct.table[i].waveform.index = wfm_index
            self.ct.table[i].amplitude0.value = amp
            self.ct.table[i].amplitude1.value = amp

    def init_ct(self):
        
        ct_schema = self.hdawg_core.awgs[0].commandtable.load_validation_schema()
        self.ct = CommandTable(ct_schema)
        
    #%% mixer_opt_funcs
    def get_power(self,fc=4e9,span=0.5e6,threshold=-50,config=False,plot=False,output=False):
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
            self.sa.setValue('Span',span)
            self.sa.setValue('Center frequency', fc)
            self.sa.setValue('Threshold',threshold)
            
        # signal = sa_get_sweep_64f(self.sa)['max']
        signal = self.sa.getValue('Signal')['y']
        power = np.max(signal)
        
       
        if plot:
            freqs = np.linspace(fc-span/2,fc+span/2,num=len(signal))
            self.power_plot(freqs, signal, power, fc=fc)
            if output:
                print(f'{power} dBm at {fc/1e9} GHz')
        return power
    
    def meas_on_power(self,inst='awg',amp=0.1):
        if meas_on_power:
            if inst == 'awg':
                self.setup_mixer_calib('awg',amp)
            elif inst == 'qa':
                self.setup_mixer_calib('qa',amp)
                
        self.get_power(fc=self.qb_pars['qb_LO']+self.qb_pars['qb_IF'],threshold=0,config=True)
        
    def config_sa(self,fc,span=0.5e6,threshold=-30):
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

            self.sa.setValue('Span',span)
            self.sa.setValue('Center frequency', fc)
            self.sa.setValue('Threshold',threshold)

    def min_leak(self,inst,mode='fine',mixer='qubit',threshold=-50,measON=False,plot=False):
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
        
        if inst == self.awg:
            device = 'dev8233'
        elif inst == self.qa:
            device = 'dev2528'
            atten = self.qb_pars['rr_atten']
        
        start = time.time()
        if mode == 'coarse':
            span = 10e-3
            dV = 1e-3
        elif mode == 'fine':
            span = 2e-3
            dV = 0.1e-3

        # generate arrays for optimization parameters
        vStart = np.zeros(2)
        for i in range(len(vStart)):
            vStart[i] = inst.get(f'/{device}/sigouts/{i}/offset')[f'{device}']['sigouts'][f'{i}']['offset']['value']
            inst.sync()
        VoltRange1 = np.arange(vStart[0]-span/2,vStart[0]+span/2,dV)
        VoltRange2 = np.arange(vStart[1]-span/2,vStart[1]+span/2,dV)

        # vStart[i] = inst.get(f'/{device}/sigouts/{channels[i]}/offset')[f'{device}']['sigouts'][f'{channels[i]}']['offset']['value']
        inst.sync()

        VoltRange1 = np.arange(vStart[0]-span/2,vStart[0]+span/2,dV)
        VoltRange2 = np.arange(vStart[1]-span/2,vStart[1]+span/2,dV)

        L1 = len(VoltRange1)
        L2 = len(VoltRange2)
        power_data = np.zeros((L1,L2))
        
        if mixer == 'qubit':
            f_LO = self.qb_pars['qb_LO']
        elif mixer == 'resonator':
            f_LO = self.qb_pars['rr_LO']
            self.update_qb_value('rr_atten', 0)
            
        self.config_sa(fc=f_LO,threshold=threshold)
            
        # Sweep individual channel voltages and find leakage
        with tqdm(total = L1*L2) as progress_bar:
            for i,V1 in enumerate((VoltRange1)):
                for j,V2 in enumerate((VoltRange2)):
                    inst.set(f'/{device}/sigouts/0/offset',V1)
                    inst.set(f'/{device}/sigouts/1/offset',V2)
                    inst.sync()
                    power_data[i,j] = self.get_power(fc=f_LO,plot=False,config=False)
                    progress_bar.update(1)

        # find index of voltage corresponding to minimum LO leakage
        argmin = np.unravel_index(np.argmin(power_data), power_data.shape)

        opt_I = VoltRange1[argmin[0]]
        opt_Q = VoltRange2[argmin[1]]
        # set voltages to optimal values
        inst.set(f'/{device}/sigouts/0/offset',opt_I)
        inst.set(f'/{device}/sigouts/1/offset',opt_Q)
        if inst == self.awg:
            self.update_qb_value('qb_mixer_offsets', [opt_I,opt_Q])
        elif inst == self.qa:
            self.update_qb_value('rr_mixer_offsets', [opt_I,opt_Q])
        inst.sync()
        print(f'optimal I_offset = {round(opt_I*1e3,1)} mV, optimal Q_offset = {round(1e3*opt_Q,1)} mV')

        end = time.time()
        print('%s mixer Optimization took %.1f seconds'%(mixer,(end-start)))

        # get LO leakage for optimal DC values
        OFF_power = self.get_power(fc=f_LO,threshold=threshold,plot=False)

        # if measON:
        #     offset = inst.get(f'/{device}/sigouts/0/offset')[f'{device}']['sigouts'][0]['offset']['value'][0]
        #     #get ON power
        #     inst.set(f'/{device}/sigouts/0/offset', amp)
        #     inst.sync()
        #     ON_power = self.get_power(fc=f_LO,threshold=0)
        #     inst.set(f'/{device}/sigouts/0/offset', offset)
        #     inst.sync()
        # else:
        #     pass

        if device == 'dev8233':
             element = 'qubit'
        elif device == 'dev2528':
             element = 'readout'

        if plot:
            self.plot_mixer_opt(VoltRange1, VoltRange2, power_data,cal='LO',mixer=mixer,fc=f_LO)
        
        if mixer == 'resonator':
            self.update_qb_value('rr_atten', atten)
        
    def suppr_image(self,inst,mode='fine',mixer='qubit',threshold=-50,f_LO=3.875e9,f_IF=50e6,amp=1,
                    sb='lsb',gen=0,plot=True):
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
                gen (int): The oscillator used for modulation.
                plot (boolean): Whether or not to plot the image power as a function the parameters.
            """

            if inst == self.awg:
                device = 'dev8233'
            elif inst == self.qa:
                device = 'dev2528'
                atten = self.qb_pars['rr_atten']
                
            if sb == 'lsb':
                f_im = f_LO - f_IF
            elif sb == 'usb':
                f_im = f_LO + f_IF

            start = time.time()
            if mode == 'coarse':
                span_amp= 10.001e-2
                da = 20e-3
                span_phi = 3.1
                dp = 0.5
            elif mode == 'fine':
                span_amp = 5.001e-3
                da = 5e-4
                span_phi = 0.5
                dp = 0.015

            # get current values of phase and amplitude
            a0 = self.qb_pars['qb_mixer_imbalance'][0]
            p0 = self.qb_pars['qb_mixer_imbalance'][1]
            # generate arrays for optimization parameters based on current values of phi and a used
            phiArr = np.arange(p0-span_phi/2,p0+span_phi/2,dp)
            ampArr = np.arange(a0-span_amp/2,a0+span_amp/2,da)

            # upload and run AWG sequence program
            self.setup_mixer_calib(inst,amp=amp)
            self.awg.set('/dev8233/oscs/0/freq',f_IF)
            self.enable_awg(inst,enable=1)
            
            
            L1 = len(ampArr)
            L2 = len(phiArr)
            power_data = np.zeros((L1,L2))

            self.config_sa(fc=f_im,threshold=threshold)

            # Sweep individual channel voltages and find leakage
            with tqdm(total = L1*L2) as progress_bar:
                for i,g in enumerate((ampArr)):
                    for j,phi in enumerate((phiArr)):
                        self.IQ_imbalance(g, phi)
                        inst.sync()
                        power_data[i,j] = self.get_power(fc=f_im,threshold=threshold)
                        progress_bar.update(1)
            
            self.enable_awg(inst,enable=0)
            # find index of voltage corresponding to minimum LO leakage
            argmin = np.unravel_index(np.argmin(power_data), power_data.shape)
            
            opt_amp = ampArr[argmin[0]]
            opt_phi = phiArr[argmin[1]]
            # set voltages to optimal values
            self.IQ_imbalance(g=opt_amp, phi=opt_phi)
            self.update_qb_value('qb_mixer_imbalance', [opt_amp,opt_phi])
            inst.sync()
            print(f'g = {round(opt_amp*1e2,1)}, phi = {round(opt_phi,3)}')

            end = time.time()
            print('%s mixer Optimization took %.1f seconds'%(mixer,(end-start)))

            if plot:
                self.plot_mixer_opt(ampArr*1e2,phiArr,  power_data,cal='SB',mixer='qubit',fc=f_im)

    #%% utilities
    def inst_init(self):
        
        session = Session('localhost')
        self.hdawg_core = session.connect_device('DEV8233')
        self.qa_awg_core = session.connect_device('DEV2528')
        self.awg,device_id,_ = create_api_session('dev8233',api_level=6)
        self.qa,device_id,_ = create_api_session('dev2528',api_level=6)
        
        instfuncs.set_qb_LO(self.qb_pars['qb_LO'])
        instfuncs.set_rr_LO(self.qb_pars['rr_LO'])
        instfuncs.set_attenuator(self.qb_pars['rr_atten'])
        
        self.sa = instfuncs.init_sa()
        
    def enable_awg(self,inst, enable=1):
        '''
        enable/disable AWG

        daq:             instrument (awg or qa)
        device:          device ID
        enable:     0:   disable AWG
                    1:   enable AWG
        '''
        if inst == self.awg:
            inst.syncSetInt('/dev8233/awgs/0/enable', enable)
        elif inst == self.qa:
            inst.syncSetInt('/dev2528/awgs/0/enable', enable)
        
    def upload_to_awg(self):
        
        
        with self.hdawg_core.set_transaction():
            # upload sequence file
            try:
                self.hdawg_core.awgs[0].load_sequencer_program(self.sequence)
            except:
                print(self.sequence)
            # upload waveforms
            self.hdawg_core.awgs[0].write_to_waveform_memory(self.sequence.waveforms)
            # validate waveforms
            # waveforms = self.hdawg_core.awgs[0].read_from_waveform_memory()
            # waveforms.validate(self.hdawg_core.awgs[0].waveform.descriptors())
            # upload command table if applicable
            if self.exp != 'spectroscopy' and self.exp != 'single-shot':
                self.hdawg_core.awgs[0].commandtable.upload_to_device(self.ct)
            
        self.awg.sync()
        
    def get_xdata_frm_ct(self):
        
        x_array = np.zeros((1,self.n_steps))
        ct = self.hdawg_core.awgs[0].commandtable.load_from_device()
        for i in range(self.n_steps):
            if self.exp != 'p-rabi':
                value = ct.table[i].waveform.length
                x_array[0][i] = int(value)
            elif self.exp == 'p-rabi':
                value = ct.table[i].amplitude0.value
                x_array[0][i] = value
            else:
                raise Exception(f'No command table used for {self.exp}')
            
         
        return x_array
            
    # def create_and_compile_awg(self,inst,device_id, awg_program, seqr_index= 0, timeout=1,verbose=0):
    #     """
    #     Compiles and uploads the sequence file into the AWG

    #     Args:
    #         inst: instrument to set awg sequence file to (awg or qa)
    #         awg_program (string): Sequence file.
    #         seqr_index (int, optional): Which AWG to upload the sequence file to. Defaults to 0.
    #         timeout (float, optional): How long to wait for sequence file upload before time out. Defaults to 1.
    #         verbose (TYPE, optional): DESCRIPTION. Defaults to 0.

    #     Raises:
    #         Exception: DESCRIPTION.

    #     Returns:
    #         None.

    #     """

    #     awgModule = inst.awgModule()
    #     awgModule.set('device', device_id)
    #     awgModule.set('index', seqr_index)
    #     awgModule.execute()
    #     """Compile and upload awg_program as .elf file"""
    #     if verbose==0:
    #         # print("Starting compilation.")
    #         awgModule.set('compiler/sourcestring', awg_program)
    #         compilerStatus = -1
    #         while compilerStatus == -1:
    #             compilerStatus = awgModule.getInt('compiler/status')
    #             time.sleep(0.1)
    #         compilerStatusString = awgModule.getString('compiler/statusstring')
    #         # print(f"Compiler messages:\n--------------\n{compilerStatusString}\n--------------")
    #         if compilerStatus == 1: # compilation failed
    #             print(awg_program)
    #             raise Exception("Compilation failed.")
    #         # if compilerStatus == 0:
    #         #     print("Compilation successful with no warnings.")
    #         # if compilerStatus == 2:
    #         #     print("Compilation successful with warnings.")
    #         # print("Waiting for the upload to the instrument.")
    #         elfProgress = 0
    #         elfStatus = 0
    #         lastElfProgressPrc = None
    #         while (elfProgress < 1.0) and (elfStatus != 1):
    #             elfProgress = awgModule.getDouble('progress')
    #             elfStatus = awgModule.getInt('elf/status')
    #             elfProgressPrc = round(elfProgress * 100);
    #             if elfProgressPrc != lastElfProgressPrc:
    #                 # print(f'Upload progress: {elfProgressPrc:2.0f}%')
    #                 lastElfProgressPrc = elfProgressPrc
    #             time.sleep(0.1)
    #     else:
    #         print("Starting compilation.")
    #         awgModule.set('compiler/sourcestring', awg_program)
    #         compilerStatus = -1
    #         while compilerStatus == -1:
    #             compilerStatus = awgModule.getInt('compiler/status')
    #             time.sleep(0.1)
    #         compilerStatusString = awgModule.getString('compiler/statusstring')
    #         print(f"Compiler messages:\n--------------\n{compilerStatusString}\n--------------")
    #         if compilerStatus == 1: # compilation failed
    #             raise Exception("Compilation failed.")
    #         if compilerStatus == 0:
    #             print("Compilation successful with no warnings.")
    #         if compilerStatus == 2:
    #             print("Compilation successful with warnings.")
    #         print("Waiting for the upload to the instrument.")
    #         elfProgress = 0
    #         elfStatus = 0
    #         lastElfProgressPrc = None
    #         while (elfProgress < 1.0) and (elfStatus != 1):
    #             elfProgress = awgModule.getDouble('progress')
    #             elfStatus = awgModule.getInt('elf/status')
    #             elfProgressPrc = round(elfProgress * 100);
    #             if elfProgressPrc != lastElfProgressPrc:
    #                 print(f'Upload progress: {elfProgressPrc:2.0f}%')
    #                 lastElfProgressPrc = elfProgressPrc
    #             time.sleep(0.1)
    #     if elfStatus == 0 and verbose == 1:
    #         print("Upload to the instrument successful.")
    #     if elfStatus == 1:
    #         raise Exception("Upload to the instrument failed.")


    def IQ_imbalance(self,g,phi,amp=0.9):
        """
        Calculates the correction matrix (C**-1) and applies it to the AWG

        Args:
            awg: awg instance to apply correction
            g (TYPE): amplitude imbalance.
            phi (TYPE): phase imbalance in degrees

        Returns:
            list: correction matrix.

        """
        # phi = phi*pi/180
        # c = np.cos(phi)
        # s = np.sin(phi)
        # N = amp / ((1-g**2)*(2*c**2-1))
        # corr = np.array([float(N * x) for x in [(1-g)*c, (1+g)*s, (1-g)*s, (1+g)*c]]).reshape(2,2)
        # # print(corr)

        # for i in range(2):
        #     for j in range(2):
        #         self.awg.set(f'/dev8233/awgs/0/outputs/{j}/gains/{i}', corr[i,j])
        self.awg.set('/dev8233/awgs/0/outputs/0/gains/0', amp*(1+g/2))
        self.awg.set('/dev8233/awgs/0/outputs/1/gains/0', -amp*(1-g/2))
        self.awg.set('/dev8233/sines/1/phaseshift',90+phi)
        
                
    def save_data(self,project,device_name,par_dict={},data=[]):
        
        dir_path = f'D:\\{project}\\{device_name}\\{self.exp}-data'

        if os.path.exists(dir_path):
            pass
        else:
            print(f'Directory not found; making new directory at {dir_path}')
            os.makedirs(dir_path)
            
        try:
            latest_file = max(glob.glob(os.path.join(dir_path, '*.csv')), key=os.path.getmtime)
            self.iteration = int(re.findall(r'\d+', latest_file)[1]) + 1
        except:
            self.iteration = 1
            
        filename = f'\\{self.exp}'+f'_data_{self.iteration}.csv'

        with open(dir_path+filename,"w",newline="") as datafile:
            writer = csv.writer(datafile)
            writer.writerow(self.qb_pars.keys())
            writer.writerow(self.qb_pars.values())
            writer.writerow(self.exp_pars.keys())
            writer.writerow(self.exp_pars.values())
            for i in range(data.shape[0]):
                writer.writerow(data[i][:])
    
    def calc_steps(self,verbose=1):
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
            n_points (int): number of waveform points.
            n_steps (int): number of sequence steps.
            t0 (int): sequence step size in units of samples 1/fsAWG.
            dt (int): starting sequence point in units of samples.

        """
        t0 = self.roundToBase(self.exp_pars['fsAWG']*self.exp_pars['x0'])
        dt = self.roundToBase(self.exp_pars['fsAWG']*self.exp_pars['dx'])
        n_points = self.roundToBase(self.exp_pars['xmax']*self.exp_pars['fsAWG']) # this ensures there is an integer number of time points
        n_steps = int((n_points-t0)/dt) # 1 is added to include the first point
        tmax = dt*n_steps
       
        if verbose == 1:
            print("dt is %.1f ns (%d pts) ==> f_s = %.1f MHz \nn_points = %d | n_steps is %d | Pulse length start = %.1f ns (%d pts)" %(dt/self.exp_pars['fsAWG']*1e9,dt,1e-6*self.exp_pars['fsAWG']/dt,n_points,n_steps,t0*1e9/self.exp_pars['fsAWG'],t0))
        else:
            pass
        
        if n_steps > 1024:
            raise Exception('Error: The maximum number of steps is 1024')

        return n_points,n_steps,t0,tmax,dt
    
    def update_qb_value(self,key,value):
        
        if key == 'gauss_len':
            value = self.roundToBase(value)
        else:
            pass
        
        print(f'Updating {key} to {value}')
        self.qb_pars[key] = value
        # self.make_config(self.qb_pars)

        if key == 'qb_LO':
            instfuncs.set_qb_LO(value)
        elif key == 'rr_LO':
            instfuncs.set_rr_LO(value)
        elif key == 'rr_atten':
            instfuncs.set_attenuator(value)
        elif key == 'qb_freq':
            self.update_qb_value('qb_IF',self.qb_pars['qb_freq']-self.qb_pars['qb_LO'])
            self.awg.set('/dev8233/oscs/0/freq',self.qb_pars['qb_IF'])
        elif key == 'qb_mixer_imbalance':
            self.IQ_imbalance(g=value[0], phi=value[1])
        
        self.write_pars()
        
    def update_pi(self,pi_amp):
        self.update_qb_value('pi_len', self.qb_pars['gauss_len'])
        self.update_qb_value('pi_amp',pi_amp)
        self.update_qb_value('pi_half_amp',pi_amp/2)
        
    def update_exp_value(self,key,value):
        print(f'Updating {key} to {value}')
        self.exp_pars[key] = value
        # self.make_config(self.exp_pars)

    def write_pars(self):
        with open(f'{self.name}_pars.json', "w") as outfile:
            json.dump(self.qb_pars, outfile)
    
    def remove_key(self, key):
        print(f'Removing {key} from pars')
        del self.pars[key]
        self.write_pars()

    def add_key(self, key, value):
        print(f'Adding {key} = {value} to pars')
        self.pars[key] = value
        self.write_pars()

    def roundToBase(self,n_points,base=16):
        '''Make the AWG happy by uploading a wfm whose points are multiple of 16'''
        y = int(base*round(n_points/base))
        if y==0:
            y = int(base*round(n_points/base+1))
            
        return y
    
    def calc_amp(self,device,amp):
        
        awg_range = self.awg.get(f'/{device}/sigouts/0/range',flat=True)[f'/{device}/sigouts/0/range']['value'][0]
        if 2*amp*awg_range > awg_range:
            raise ValueError(f'Waveform amplitude [{amp}] exceeds {device} range [{awg_range}]')
        else:
            amp = awg_range*amp/2
            
        return amp

    def odd(self,n):
        return range(1,n,2)

    def even(self,n):
        return range(0,n,2)
    
    def make_config(self, awg,pars):
        # gauss_wf_4ns = self.delayed_gauss()

        self.config = {


            "controllers": {
                "awg": {

            "waveforms": {
                "zero_wf": {"type": "constant", "sample": 0.0},
                "const_wf": {"type": "constant", "sample": pars['amp_q']},
                "const_wf_rr": {"type": "constant", "sample": pars['amp_r']},
                "gaussian_wf": {"type": "arbitrary", "samples": [float(arg) for arg in pars['gauss_amp'] * gaussian(pars['gauss_len'], pars['gauss_len']/5)]},
                # "gaussian_4ns_wf": {"type": "arbitrary", "samples": gauss_wf_4ns},
                "ro_wf1": {"type": "constant", "sample": pars['amp_r']},
                "pi_wf_i1": {"type": "arbitrary", "samples": [float(arg) for arg in pars['pi_amp'] * gaussian(pars['pi_len'], pars['pi_len']/5)]},
                "pi_wf_q1": {"type": "constant", "sample": 0.0},
                "pi_half_wf_i1": {"type": "arbitrary", "samples": [float(arg) for arg in pars['pi_half_amp'] * gaussian(pars['pi_half_len'], pars['pi_half_len']/5)]},
                "pi_half_wf_q1": {"type": "constant", "sample": 0.0},
                "arb_wfm": {"type": "arbitrary", "samples": [0.2]*10+[0.3]*10+[0.25]*20},
            },}}}

            
    # def create_wfms(self,awg):
        
        
        
    def set_config(self,awg,qa):
        
        self.IQ_imbalance(awg, g=self.pars['qubit_mixer_imbalance'][0], phi=self.pars['qubit_mixer_imbalance'][1])
        
        awg_setting = [
            ['/dev8233/oscs/0/freq', self.pars['qubit_IF']], # sets the oscillator freq
            ['/dev8233/sigouts/0/offset', self.pars['qubit_mixer_offsets'][0]],
            ['/dev8233/sigouts/1/offset', self.pars['qubit_mixer_offsets'][1]],
        ]
        print('Updating settings on HDAWG')
        awg.set(awg_setting)
        awg.sync()
    
        qa_setting = [
            ['/dev2528/sigouts/0/offset', self.pars['rr_mixer_offsets'][0]],
            ['/dev2528/sigouts/1/offset', self.pars['rr_mixer_offsets'][1]],
        ] 
        print('Updating settings on UHFQA')
        qa.set(qa_setting)
        qa.sync()
        
    def make_wfm_file(self,filename,wfm_data):
        
        path = "C:/Users/LFL/Documents/Zurich Instruments/LabOne/WebServer/awg/waves/"
        np.savetxt(path+filename+"_wfm"+".csv",wfm_data, delimiter = ",") # save file where it can be called by the AWG sequence program
        
    #%% Plot functions
    # from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    # from matplotlib.colors import LightSource
    # set some deafault
    # plt.rcParams.update(plt.rcParamsDefault)
    # for a complete set of parameters "print(plt.rcParams)"
    # sns.set_style('ticks')
    # plt.rcParams['text.usetex'] = False
    # plt.rcParams['font.family'] =  'Arial'
    # plt.rcParams['font.size'] = 16
    # plt.grid(False)
    # plt.rcParams["xtick.direction"] = "in"
    # plt.rcParams["ytick.direction"] = "in"
    # plt.rcParams["xtick.major.top"] = True
    # plt.rcParams['xtick.labelsize'] = 14
    # plt.rcParams['axes.labelsize'] = 16
    # plt.rcParams['axes.titlesize'] = 16
    # plt.rcParams["xtick.major.bottom"] = True
    # plt.rcParams["xtick.top"] = True
    # plt.rcParams["xtick.bottom"] = True
    # plt.rcParams["ytick.left"] = True
    # plt.rcParams["ytick.right"] = True
    # plt.rcParams['ytick.labelsize'] = 14
    # plt.rcParams["ytick.major.right"] = True
    # plt.rcParams["ytick.labelright"] = False
    # plt.rcParams["ytick.minor.visible"] = False

    def power_plot(self,freqs,signal,power,fc):
        plt.plot(freqs*1e-6, signal,'-')
        plt.xlabel('Frequency [MHz]')
        plt.ylabel('Power [dBm]')
        plt.show()

    def tof_plot(self,adc1,adc2):
        plt.figure()
        plt.title('time-of-flight calibration analysis')
        plt.plot(adc1)
        plt.plot(adc2)
        plt.legend(["adc1", "adc2"])
        plt.show()

    def qb_spec_plot(self,freq,I,Q,attenuation=-30,find_peaks=False):

        I = I*1e3
        Q = Q*1e3
        mag = np.abs(I+1j*Q)

        phase = np.unwrap(np.angle(I+1j*Q))
        if find_peaks:
            sigma = np.std(mag)
            print(f'Peak threshold at {np.mean(mag)+3*sigma:.1f}')
            peaks,_ = scy.signal.find_peaks(mag,height=np.mean(mag)+3*sigma,distance=200,width=3)
            try:
                for i in peaks:
                    print(f'Peaks at: {round(freq[i],5)} GHz\n')
            except:
                print('Peaks not found or do not exist.')

        fig = plt.figure(figsize=(5,4))

        # # I data
        # ax1 = fig.add_subplot(221)
        # ax1.plot(freq,I,'-o', markersize = 3, c='C0')
        # ax1.set_xlabel('Frequency (GHz)')
        # ax1.set_ylabel('I (mV)')
        # # Q data
        # ax1 = fig.add_subplot(222)
        # ax1.plot(freq,Q,'-o', markersize = 3, c='C0')
        # ax1.set_xlabel('Frequency (GHz)')
        # ax1.set_ylabel('Q (mV)')
        # Power data
        ax1 = fig.add_subplot()
        ax1.plot(freq,mag,'-o', markersize = 3, c='C0')
        ax1.set_xlabel('Frequency (GHz)')
        ax1.set_ylabel('Magnitude (mV)')
        # phase data
        # ax1 = fig.add_subplot(212)
        # ax1.plot(freq,phase,'-o', markersize = 3, c='C0')
        # ax1.set_xlabel('Frequency (GHz)')
        # ax1.set_ylabel('Phase (rad)')

        if len(peaks) == 2:
            txt = '$\omega_{01}$ = %.4f GHz\n$\omega_{12}$ = %.4f GHz\n$\\alpha$ = %.1f MHz\n$A_{qb}$ = %.2f V\nReadout Attn = %d dBm\n$f_r$ = %.4f GHz'%(freq[peaks[1]],freq[peaks[0]],(freq[peaks[0]]-freq[peaks[1]])*1e3,self.exp_pars['amp_q'],self.exp_pars['rr_attenuation'],(self.qb_pars['rr_LO']+self.qb_pars['rr_IF'])*1e-9)
        elif len(peaks) == 1:
            txt = '$\omega_{01}$ = %.4f GHz\n$A_{qb}$ = %.1f mV\nReadout Attn = %d dBm\n$f_r$ = %.4f GHz'%(freq[peaks[0]],self.exp_pars['amp_q']*1e3,self.exp_pars['rr_attenuation'],(self.qb_pars['rr_LO']+self.qb_pars['rr_IF'])*1e-9)
        else:
            txt = '$A_{qb}$ = %.1f mV\nReadout Attn = %d dBm\n$f_r$ = %.4f GHz'%(self.exp_pars['amp_q']*1e3,self.exp_pars['rr_attenuation'],(self.qb_pars['rr_LO']+self.qb_pars['rr_IF'])*1e-9)
        plt.gcf().text(1, 0.15, txt, fontsize=14)
        # fig.set_title(f'Qubit spectroscopy {self.iteration}')
        plt.tight_layout()
        plt.show()
        
    def rr_spec_plot(self,freq,I,Q,df=0.1e6,find_peaks=False):

        I = I*1e3
        Q = Q*1e3
        mag = np.abs(I+1j*Q)
        
        if find_peaks:
            sigma = np.std(mag)
            print(f'Peak threshold at {np.mean(mag)+3*sigma}')
            peaks,_ = scy.signal.find_peaks(mag,height=np.mean(mag)+3*sigma,distance=200,width=3)
            try:
                for i in peaks:
                    print(f'Peaks at: {round(freq[i],5)} GHz\n')
            except:
                print('Peaks not found or do not exist.')
                
        fc,fwhm = self.fit_res(freq,mag)

        fig = plt.figure(figsize=(5,4))

        # Power data
        ax1 = fig.add_subplot(111)
        ax1.plot(freq,mag,'-o', markersize = 3, c='C0')
        ax1.set_xlabel('Frequency (GHz)')
        ax1.set_ylabel('Magnitude (mV)')
        # phase = np.unwrap(np.angle(I+1j*Q))
        # ax2 = fig.add_subplot(212)
        # ax2.plot(freq,phase,'-o', markersize = 3, c='C0')
        # ax2.set_xlabel('Frequency (GHz)')
        # ax2.set_ylabel('Phase (deg)')

        txt = '$\omega_c$ = %.6f GHz\nFWHM = %d kHz\n$T_1$ = %d ns\n$A_d$ = %.1f mV\n$df$ = %1.f kHz'%(fc,fwhm*1e6,1e6/(np.pi*fwhm*1e6),self.exp_pars["amp_q"]*1e3,df*1e-3)
        plt.gcf().text(1, 0.15, txt, fontsize=14)
        # fig.set_title(f'{element} spectroscopy {iteration}')
        plt.tight_layout()
        plt.show()

    def init_IQ_plot(self,):
        '''initialize axes for continuous plotting on the IQ plane'''
        plot = sns.jointplot()
        plot.set_axis_labels('I [mV]', 'Q [mV]')
        plot.ax_marg_x.grid('off')
        plot.ax_marg_y.grid('off')
        plot.fig.tight_layout()
        ax = plt.gca()
        return plot, ax

    def heatplot(self,xdata, ydata, data, xlabel = "", ylabel = "", normalize=False, cbar_label = 'log mag', **kwargs):
        fig,ax = plt.figure(figsize=(4,3), dpi=300)
        if normalize:
            cbar_label += ' (normalized)'

        df = pd.DataFrame(data, columns = xdata, index = ydata)

        if normalize:
            df = df.apply(lambda x: (x-x.mean())/x.std(), axis = 1)

        cbar_options = {
            'label':    cbar_label,
            'ticks':    np.around(np.linspace(np.amin(data),np.amax(data),5),0),
            'pad':      0.1,
            'values':   np.linspace(np.amin(data),np.amax(data),1000),
            'shrink':   1.2,
            'location': 'right',
        }

        heatmap_opts = {
            'ax':    ax,
            'linewidths':  0,
            # 'xticklabels': np.linspace(min(xdata),max(xdata)+0.5,5),
            # 'yticklabels': np.linspace(min(ydata),max(ydata)+0.5,5),
            'vmin':        np.amin(data),
            'vmax':        np.amax(data)
        }

        hm = sns.heatmap(df, ax=ax,cmap = 'viridis', cbar_kws=cbar_options, **kwargs)
        hm.set_xlabel(xlabel, fontsize=12)
        hm.set_ylabel(ylabel, fontsize=12)
        hm.spines[:].set_visible(True)
        ax.tick_params(direction='out',length=0.01,width=0.5,bottom=True, top=True, left=True, right=True,labeltop=False, labelbottom=True,labelrotation=90,labelsize=8,size=8)
        plt.yticks(rotation=0)
        plt.tight_layout()
        plt.show()

        return hm;

    def plot_single_shot(self,datadict, axes=0):

        datadict = {key: value*1e3 for key,value in datadict.items()} # convert to mV
        datadict = {key: value.tolist() for key,value in datadict.items()} # convert to list

        states = []
        # for key,value in datadict.items():
        #     print(key+':'+str(len(value))+'\n')
        [states.append(r'$|g\rangle$') for i in range(len(datadict['I']))]
        [states.append(r'$|e\rangle$') for i in range(len(datadict['Iexc']))]
        data = {
                'I [mV]':   np.hstack((datadict['I'],datadict['Iexc'])),
                'Q [mV]':   np.hstack((datadict['Q'],datadict['Qexc'])),
                'States':   states
                    }
        dataF = pd.DataFrame(data=data)
        plot = sns.jointplot(data=dataF, x='I [mV]',y='Q [mV]',hue='States',ax=axes,space=0)
        plt.show()

    def plot_mixer_opt(self,par1,par2,power_data,cal='LO',mixer='qubit',fc=5e9):
        
        if cal == 'LO':
            par1 = np.around(par1*1e3,1)
            par2 = np.around(par2*1e3,1)
        else:
            par1 = np.around(par1,3)
            par2 = np.around(par2,3)
        par1 = par1.tolist()
        par2 = par2.tolist()
        df = pd.DataFrame(data=power_data,index=par1,columns=par2)

        hm = sns.heatmap(df,cbar_kws={'label': "Power [dBm]"})

        if cal == 'LO':
            hm.set_ylabel('I [mV]')
            hm.set_xlabel('Q [mV]')
        elif cal == 'SB':
            hm.set_ylabel('Gain Imbalance (%)')
            hm.set_xlabel('Phase Imbalance')

        hm.spines[:].set_visible(True)
        hm.tick_params(direction='out',length=0.01,width=0.5,bottom=True, top=False, left=True, right=True,labeltop=False, labelbottom=True,labelrotation=90,labelsize=10,size=10)
        plt.yticks(rotation=0)
        plt.tight_layout()
        if mixer == 'qubit':
            plt.title(f'Qubit Mixer {cal} Calibration at {round(fc*1e-9,4)} GHz')
        elif mixer == 'resonator':
            plt.title(f'Readout Mixer {cal} Calibration at {round(fc*1e-9,4)} GHz')
        plt.show()


    def fit_res(self,f_data,z_data,res_type='notch'):
        fc = f_data[np.argmax(z_data)]
        if res_type == 'notch':
            z_data = z_data-min(z_data)
            idx = np.argwhere(np.diff(np.sign(z_data - 0.5*max(z_data)))).flatten()
            fwhm = f_data[idx[1]] - f_data[idx[0]]

        return fc,fwhm

    def fit_data(self,x_vector,y_vector,dx=0.01,fitFunc='',verbose=0):

        '''
        fit experimental data
        sequence:          'Rabi','ramsey', 'T1' or 'T2'
        x_vector:           time data
        y_vector:           voltage data
        dx:                 sequence stepsize. Used for extracting the frequency of the data
        '''
        
        x_vector = x_vector*1e6
        y_vector = y_vector*1e3

        amp = (max(y_vector)-min(y_vector))/2
        offset = np.mean(y_vector)

        if self.exp == "rabi":
            fitFunction = self.rabi
            x_vector = x_vector*1e3
            period = 1e3/(self.extract_freq(x_vector, y_vector, dx,plot=0))
            print('Period Initial Guess: %.1f ns'%(period))
            phase = pi/2
            lb = [0.1*amp,0.1*period,0,-2*abs(offset)]
            ub = [10*amp,10*period,2*pi,2*abs(offset)]
            p0 = [amp,period,phase,offset]

        elif self.exp == "p-rabi":
            fitFunction = self.rabi
            x_vector = x_vector*1e-6
            period = 1/(self.extract_freq(x_vector, y_vector, dx,plot=0))
            print('Amplitude Initial Guess: %.3f'%(period))
            phase = 0
            lb = [0.1*amp,0.1*period,0,-2*abs(offset)]
            ub = [10*amp,10*period,2*pi,2*abs(offset)]
            p0 = [amp,period,phase,offset]

        elif self.exp == "ramsey":
            f = self.extract_freq(x_vector, y_vector,dx,plot=0)
            # print('Initial Guess for Freq:%.4f MHz'%(f))
            if x_vector[-1] > 20:
                tau = 30
            else:
                tau = 2
            phi = 0
            amp = abs(amp)
            # try:
            if fitFunc != 'envelope':
                p0 = [amp,f,phi,tau,offset]
                lb = [0.75*amp,0.1*f,-pi,0.01,-2*abs(offset)]
                ub = [2*amp,2*f,pi,100,2*abs(offset)]
                fitFunction = self.ramsey
                # fitted_pars, covar = scy.optimize.curve_fit(fitFunction, x_vector, y_vector,p0=p0,method='trf',bounds=[lb,ub],xtol=1e-12,maxfev=20e3)
            elif fitFunc == 'envelope':
                tau = 1
                env = self.get_envelope(y_vector, dx, distance=100)
                env = env(x_vector) + offset
                # env = get_envelope_LPF(x_vector, y_vector)*1e-3
                p0 = [amp,tau,offset]
                if offset < 0:
                    p0 = [amp,tau,offset]
                    lb = [0.95*amp,0.1,2*offset]
                    ub = [2*amp,15,0.5*offset]
                elif offset >= 0:
                    p0 = [amp,tau,offset]
                    lb = [0.9*amp,0.1,0.9*offset]
                    ub = [1.1*amp,15,1.1*offset]
                fitFunction = self.decay
                y_vector = env
                # fitted_pars, covar = scy.optimize.curve_fit(decay, x_vector, env,p0=p0,method='trf',bounds=[lb,ub],xtol=1e-12,maxfev=20e3)
        elif self.exp == "echo":
            if x_vector[-1] < 10:
                tau = 2
                tau_ub = 20
            else:
                tau = 20
                tau_ub = 300
            amp = y_vector[0] - y_vector[-1]
            p0 = [amp, tau, offset]
            amp_bounds = [0.95 * amp, 1.05 * amp]
            off_bounds = [0.95 * offset, 1.05 * offset]
            lb = [min(amp_bounds), 0.1, min(off_bounds)]
            ub = [max(amp_bounds), tau_ub, max(off_bounds)]
            # if offset < 0:
            #     lb = [0.95*amp,0.1,1.05*offset]
            #     ub = [1.05*amp,tau_ub,0.95*offset]
            # elif offset >= 0:
            #     lb = [0.95*amp,0.1,0.95*offset]
            #     ub = [1.05*amp,tau_ub,1.05*offset]
            fitFunction = self.decay
            # fitted_pars, covar = scy.optimize.curve_fit(decay, x_vector, y_vector,p0=p0,method='trf',bounds=[lb,ub],xtol=1e-12,maxfev=6000)
        elif self.exp == "T1":
            tau = 2
            amp = y_vector[0] - y_vector[-1]
            if amp < 0:
                p0 = [amp,tau,offset]
                lb = [10*amp,0.1,-2*abs(offset)]
                ub = [0.5*amp,300,2*abs(offset)]
            elif amp >= 0:
                p0 = [amp,tau,offset]
                lb = [0.5*amp,0.1,-2*abs(offset)]
                ub = [10*amp,300,2*abs(offset)]
            fitFunction = self.decay
            # fitted_pars, covar = scy.optimize.curve_fit(, x_vector, y_vector,p0=p0,method='trf',bounds=[lb,ub],xtol=1e-12,maxfev=6000)

        fitted_pars, covar = scy.optimize.curve_fit(fitFunction, x_vector, y_vector,p0=p0,method='trf',bounds=[lb,ub],xtol=1e-12,maxfev=40e3)
        error = np.sqrt(abs(np.diag(covar)))

        if verbose == 1:
            print('-'*100)
            print('Lower Bounds:',np.around(lb,1))
            print('Initial Guess:',np.around(p0,1))
            print('Upper Bounds:',np.around(ub,1))
            print('Best Fit Pars:',np.around(fitted_pars,1))
            print('Error:',np.around(error,1))
            print('-'*100)
        else:
            pass

        return fitted_pars,error


    def plot_data(self,x_vector,y_vector,fitted_pars=np.zeros(7),fitFunc='',savefig=True):

        x_vector = x_vector*1e6
        y_vector = y_vector*1e3
        qb_drive_freq = (self.qb_pars["qb_LO"] +self.qb_pars["qb_IF"])*1e-9
        fig, ax = plt.subplots()
        
        if self.exp == "p-rabi":
            # ax = self.make_rabi_plot(ax,x_vector,y_vector,fitted_pars)
            x_vector = x_vector*1e-6
            ax.plot(x_vector, y_vector, '-o', markersize = 3, c='C0')
            ax.set_ylabel('Digitizer Voltage (mV)')
            ax.set_xlabel('Pulse Amplitude')
            ax.plot(x_vector,self.rabi(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2],fitted_pars[3]),'r')
            ax.set_title(f'Rabi Measurement {self.iteration:03d}')
            textstr = f'$\omega_d$ = {qb_drive_freq:.4f} GHz\n$N$ = {self.exp_pars["n_avg"]}\n$A_\pi$ = {fitted_pars[1]/2:.3f}\n$T_\pi$ = {self.qb_pars["gauss_len"]/2.4:.1f} ns'

        if self.exp == "rabi":
            x_vector = x_vector*1e3
            ax.plot(x_vector, y_vector, '-o', markersize = 3, c='C0')
            ax.set_ylabel('Digitizer Voltage (mV)')
            ax.set_xlabel('Pulse Duration (ns)')
            ax.plot(x_vector,self.rabi(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2],fitted_pars[3]),'r')
            ax.set_title(f'Rabi Measurement {self.iteration:03d}')
            textstr = f'$\omega_d$ =  {qb_drive_freq:.4f} GHz\n$A_q$ = {self.exp_pars["amp_q"]:.2f} V\n$N$ = {self.exp_pars["n_avg"]}'

        elif self.exp == "ramsey":

            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            fontSize = 16
            tickSize = 11
            markersize= 4
            linewidth = 2
            ax1.plot(x_vector, y_vector, 'o', markersize = markersize, c='C0')
            ax1.set_ylabel('Digitizer Voltage (mV)',fontsize=fontSize)
            ax1.set_xlabel('Pulse Separation ($\mu$s)')
            for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
                label.set_fontsize(tickSize)
            if fitFunc == 'envelope':
                ax1.plot(x_vector,self.decay(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2]),'r',linewidth=linewidth)
            else:
                ax1.plot(x_vector,self.ramsey(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2],fitted_pars[3],fitted_pars[4]),'r',linewidth=linewidth)
            ax1.set_title(f'Ramsey Measurement {self.iteration:03d}')
            textstr = f'$T_pi/2$={self.qb_pars["pi_half_len"]:.1f} ns\n$\omega_d$ = {qb_drive_freq:.4f} GHz\n$\Delta$ = {fitted_pars[1]:.2f} MHz\n$T_2^R$ = {fitted_pars[3]:.1f} $\mu$s\n$N$ = {self.exp_pars["n_avg"]}'

        elif self.exp == "echo":

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x_vector, y_vector, '-o', markersize = 3, c='C0')
            ax.set_ylabel('Digitizer Voltage (mV)')
            ax.set_xlabel('Pulse Separation ($\mu$s)')
            ax.plot(x_vector,self.decay(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2]),'r')
            textstr = f'$T_pi2$ = {self.qb_pars["pi_len"]:.1f} ns\n$\omega_d$ =  {qb_drive_freq:.4f} GHz\n$T_2^E$ = {fitted_pars[1]:.1f} $\mu$s\n$N$ = {self.exp_pars["n_avg"]}'
            ax.set_title(f'Echo Measurement {self.iteration:03d}')


        elif self.exp == "T1":
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x_vector, y_vector, '-o', markersize = 3, c='C0')
            ax.set_ylabel('Digitizer Voltage (mV)')
            ax.set_xlabel('Delay ($\mu$s)')
            ax.plot(x_vector,self.decay(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2]),'r')
            textstr = f'$T_\pi$= {self.qb_pars["pi_len"]/2.4:.1f} ns\n$A_\pi$ = {self.qb_pars["pi_amp"]*1e3:.1f} mV\n$\omega_d$ = {qb_drive_freq:.4f} GHz\n$T_1$ = {fitted_pars[1]:.1f} $\mu$s\n$N$ = {self.exp_pars["n_avg"]}'
            ax.set_title(f'T1 Measurement {self.iteration:03d}')

        plt.gcf().text(0.95, 0.15, textstr, fontsize=14)

        plt.tick_params(axis='both',direction='in',bottom=True, top=True, left=True, right=True,size=8)

        if savefig:
            plt.savefig(f'D:\\{project}\\{device_name}\\{self.exp}-data\\fig_{self.iteration:03d}.png',dpi='figure')

        # plt.show()

        return fig

    def extract_data(self,sweep,B0,nu,tauk,meas_device='CandleQubit_6',nMeasurements=100,nBackMeasurements=100,fileformat='new',sequence='ramsey',iteration=1):

        # get data
        if fileformat == 'old':
            filename = 'B0_%d_uV_nu_%d_kHz_tau_%d_ns' %(round(B0*1e6),round(nu*1e3),round(tauk*1e3))
            datafile = "E:\\generalized-markovian-noise\\%s\\sweep_data\\ramsey\\%s\\data\\data_%s.csv"%(meas_device,sweep,filename)
            tdata_background = pd.read_csv(datafile,on_bad_lines='skip',skiprows=3,header=None,nrows=1).dropna(axis='columns').to_numpy(np.float64)[0]
            tdata = pd.read_csv(datafile,on_bad_lines='skip',skiprows=5,header=None,nrows=1).dropna(axis='columns').to_numpy(np.float64)[0]
            ydata_background = pd.read_csv(datafile,on_bad_lines='skip',skiprows=7,header=None,nrows=nBackMeasurements).dropna(axis='columns').to_numpy(np.float64)
            ydata = pd.read_csv(datafile,on_bad_lines='skip',skiprows=7+2*(nBackMeasurements+1),header=None,nrows=nMeasurements).dropna(axis='columns').to_numpy(np.float64)
            # get exp parameters
            pars = pd.read_csv(datafile,on_bad_lines='skip',header=[0],nrows=1)
            keys = pars.keys()
            values = pars.values
            dictionary = {"mu":loads(pars.loc[0].at['AC_pars'])[0],
                          "sigma":loads(pars.loc[0].at['AC_pars'])[1],
                          "B0":loads(pars.loc[0].at['RT_pars'])[0],
                          # "nu": loads(pars.loc[0].at['RT_pars'])[1],
                          "tauk":loads(pars.loc[0].at['RT_pars'])[1]}
        elif fileformat == 'new':
            # cols_background = [i for i in range(0,94)]
            # cols = [i for i in range(0,250)]
            filename = 'B0_%d_uV_nu_%d_Hz_tau_%d_ns' %(round(B0*1e6),round(nu*1e3),round(tauk*1e3))
            # filename = 'B0_%d_uV_nu_%d_Hz_tau_%d_ns_%d' %(round(B0*1e6),round(nu*1e3),round(tauk*1e3),iteration)
            datafile = "E:\\generalized-markovian-noise\\%s\\sweep_data\\%s\\%s\\data\\data_%s.csv"%(meas_device,sequence,sweep,filename)
            tdata_background = pd.read_csv(datafile,on_bad_lines='skip',skiprows=3,header=None,nrows=1).to_numpy(np.float64)[0]
            tdata = pd.read_csv(datafile,on_bad_lines='skip',skiprows=5,header=None,nrows=1).to_numpy(np.float64)[0]
            ydata_background = pd.read_csv(datafile,on_bad_lines='skip',skiprows=7,header=None,nrows=nBackMeasurements).to_numpy(np.float64)
            ydata = pd.read_csv(datafile,on_bad_lines='skip',skiprows=7+nBackMeasurements+1,header=None,nrows=nMeasurements).to_numpy(np.float64)
            # get exp parameters
            pars = pd.read_csv(datafile,on_bad_lines='skip',header=[0],nrows=1)
            keys = pars.keys()
            values = pars.values
            dictionary = dict(zip(keys,values[0]))
        #print(dictionary)

        return tdata_background, ydata_background, tdata, ydata, dictionary


    def rabi(self,x, amp,period,phase,offset):
        return amp*np.cos(2*pi*x/period+phase)+offset
    
    def make_plot(self,ax,x_vector,y_vector,fitted_pars):
        ax.plot(x_vector, y_vector, '-o', markersize = 3, c='C0')
        ax.set_ylabel('Digitizer Voltage (mV)')
        ax.set_xlabel('Pulse Duration (ns)')
        ax.plot(x_vector*1e3,self.rabi(x_vector*1e3, fitted_pars[0], fitted_pars[1], fitted_pars[2],fitted_pars[3]),'r')
        ax.set_title(f'Rabi Measurement {self.iteration:03d}')
        
        return ax

    def ramsey(self,x,amp,f,phase,tau,offset):
        return amp*np.cos(2*pi*f*x+phase)*np.exp(-x/tau)+offset
    

    def beats(self,x,amp,f1,f2,phase1,phase2,tau,offset):
        return amp*np.cos(pi*(f1+f2)*x+phase1)*np.cos(pi*(f2-f1)*x+phase2)*np.exp(-x/tau)+offset

    def decay(self,x,amp,tau,offset):
        return amp*np.exp(-x/tau)+offset

    def mod_cos(self,x,amp,B0,nu,phi1,phi2,tau,offset):
        return amp*np.cos(B0/nu*np.sin(2*np.pi*nu*x+phi1)+phi2)*np.exp(-x/tau)+offset
        # return amp*np.cos(np.cos(nu*x)*f*x)*np.exp(-x/tau)+offset

    def mod_dec(self,x,amp1,f1,phi1,tau1,amp2,phi2,tau2,offset):
        return amp1*np.cos(2*np.pi*f1*x+phi1)*np.exp(-x/tau1)+ amp2*np.sin(2*np.pi*f1*x+phi2)*np.exp(-x/tau2)+offset
        # return amp*np.cos(np.cos(nu*x)*f*x)*np.exp(-x/tau)+offset

    def extract_freq(self,t_vector,y_vector,dt,plot=0):
        N = len(t_vector)
        yf = scy.fft.fft(y_vector-np.mean(y_vector))
        xf = scy.fft.fftfreq(N,dt)[:round(N/2)]
        # print(len(xf))
        psd = 2.0/N * np.abs(yf[:round(N/2)])
        # print(len(psd))
        # print(psd)
        index_max = np.argmax(psd)
        if plot == 1:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(xf,psd)
            ax.set_xlabel('Frequency (MHz)')
            ax.set_ylabel('Power')
        # print(index_max)

        return xf[index_max]

    def get_envelope(self,sig,dt, distance):
        # split signal into negative and positive parts
        sig = sig - np.mean(sig)
        u_x = np.where(sig > 0)[0]
        l_x = np.where(sig < 0)[0]
        u_y = sig.copy()
        u_y[l_x] = 0
        l_y = -sig.copy()
        l_y[u_x] = 0

        # find upper and lower peaks
        u_peaks, _ = scipy.signal.find_peaks(u_y, distance=distance)
        l_peaks, _ = scipy.signal.find_peaks(l_y, distance=distance)
        # use peaks and peak values to make envelope
        u_x = u_peaks
        u_y = sig[u_peaks]
        l_x = l_peaks
        l_y = sig[l_peaks]

        # add start and end of signal to allow proper indexing
        end = len(sig)
        u_x = np.concatenate(([0],u_x, [end]))*dt
        u_y = np.concatenate(([sig[0]],u_y, [sig[-1]]))
        l_x = np.concatenate(([0],l_x, [end]))*dt
        l_y = np.concatenate(([min(sig)],l_y, [np.mean(sig)]))
        # create envelope functions
        u = scipy.interpolate.interp1d(u_x, u_y,kind='cubic',fill_value="extrapolate")
        # l = scipy.interpolate.interp1d(l_x, l_y,kind='cubic')
        return u

    def get_envelope_LPF(self,x,sig):

        N = len(sig)
        Tmax = x[-1]
        cutoff = 100e6
        fs = N/Tmax

        env = butter_lowpass_filter(sig, cutoff, fs)
        plt.plot(x,env)
        plt.show()

        return env

    def butter_lowpass(self,cutoff, fs, order=5):
        nyq = 0.5 * fs
        normal_cutoff = cutoff / nyq
        sos = butter(order, normal_cutoff, btype='low', output = 'sos', analog=False)
        return sos

    def butter_lowpass_filter(self,data,cutoff,fs):
        sos = butter_lowpass(cutoff, fs, order=5)
        y = sp.signal.sosfilt(sos, data)
        return y

    def Volt2dBm(self,data):

        return 10*np.log10(1e3*data**2/50)

    def Watt2dBm(self,x):
        '''
        converts from units of Watts to dBm
        '''
        return 10.*np.log10(x*1000.)