# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 15:53:55 2023

Description: A script to initialize and setup instruments other than the zurich instruments
@author: Evangelos Vlachos (evlachos@usc.edu)

"""

import sys
from VISAdrivers.sa_api import *
sys.path.append("D:\Program Files\Keysight\Labber\Script")
sys.path.append(r"C:\Users\lfl\measurements_test")
import Labber
import json

client = Labber.connectToServer()

def set_LO(inst,freq,sweep=False):
    
    if not sweep:
        print(f'Setting {inst} LO to {round(freq*1e-9,5)} GHz')
    # initialize LO
    LO = client.connectToInstrument('BNC 845 Signal Generator', dict(name=inst, startup = 'Do nothing'))
    LO.startInstrument()
    LO.setValue('Frequency', freq)
    LO.setValue('Output',True)

def get_LO(inst):
    # initialize qubit LO
    LO = client.connectToInstrument('BNC 845 Signal Generator', dict(name=inst, startup = 'Do nothing'))
    LO.startInstrument()
    return LO.getValue('Frequency')

# def set_rr_LO(freq,sweep=False):
#     if not sweep:
#         print(f'Setting readout LO to {round(freq*1e-9,5)} GHz')
#     # initialize readout LO
#     rrLO = client.connectToInstrument('BNC 845 Signal Generator', dict(name='readout', startup = 'Get config'))
#     rrLO.startInstrument()
#     rrLO.setValue('Frequency', freq)
#     rrLO.setValue('Output',True)

# def get_rr_LO():
#     # initialize qubit LO
#     rrLO = client.connectToInstrument('BNC 845 Signal Generator', dict(name='readout', startup = 'Get config'))
#     rrLO.startInstrument()
#     return rrLO.getValue('Frequency')

def set_output(inst,output=True):
    # initialize qubit LO
    LO = client.connectToInstrument('BNC 845 Signal Generator', dict(name=inst, startup = 'Get config'))
    LO.startInstrument()
    return LO.setValue('Output',output)

def set_attenuator(attenuation):
    # initialize digital attenuator
    attn = client.connectToInstrument('Vaunix Lab Brick Digital Attenuator',dict(name='Readout',address='26777'))
    attn.startInstrument()
    attn.setValue('Attenuation',attenuation)

def get_attenuation():
    # initialize digital attenuator
    attn = client.connectToInstrument('Vaunix Lab Brick Digital Attenuator',dict(name='Readout',address='26777'))
    attn.startInstrument()

    return attn.getValue('Attenuation')

def init_sa():
    sa = client.connectToInstrument('SignalHound SpectrumAnalyzer', dict(name='SA',startup='Get Config'))
    sa.startInstrument()
    sa.setValue('Span',0.5e6)
    sa.setValue('Bandwidth',1e3)
    sa.setValue('Threshold',-20)
    # sa = sa_open_device()["handle"]
    # sa_config_acquisition(device = sa, detector = SA_AVERAGE, scale = SA_LOG_SCALE)
    # sa_config_gain_atten(sa, SA_AUTO_ATTEN, SA_AUTO_GAIN, True)
    # sa_config_sweep_coupling(device = sa, rbw = 1e2, vbw = 1e2, reject=0)

    return sa

# initialize Keithley
# SC = client.connectToInstrument('Keithley 2400 SourceMeter',dict(interface='GPIB',address='24'))
# SC.startInstrument()

if __name__ == '__main__':
    # initialize spectrum analyzer
    pass
