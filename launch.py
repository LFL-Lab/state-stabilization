# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 14:57:30 2023

File run at startup of Spyder.Needs to be different for every project
@author: Evangelos

"""


from tqdm import trange
import gc
import time
import importlib.util
import json
import sys, os
import VISAdrivers.LO845 as LO
from VISAdrivers.LabBrick_LMS_Wrapper import LabBrick_Synthesizer
from VISAdrivers.vaunix_attenuator_wrapper import VaunixAttenuator
from VISAdrivers.sa_api import *
import numpy as np
from zhinst.utils import create_api_session
import experiment_funcs as expf
import matplotlib.pyplot as plt
import csv
import glob
import scipy as scy
import plot_functions as pf
import pandas as pd
import zhinst.toolkit as zt
import instrument_init as inst
# matplotlib inline
# autoreload 2

'''Instruments and connection'''
qa_id = 'dev2528'
awg_id = 'dev8233'

# Instrument Addresses
qubitLO_IP = "USB0::0x03EB::0xAFFF::621-03A100000-0538::0::INSTR"
readoutLO_IP = "USB0::0x03EB::0xAFFF::621-03A100000-0519::0::INSTR"
# readout_attn = 26551

# initialize instruments as python objects
qubitLO = LO.LO(address=qubitLO_IP,reset=False)
readoutLO = LO.LO(readoutLO_IP,reset=False)
attn = VaunixAttenuator()
attn.getNumDevices()
attn.initDevice(26920)
sa = sa_open_device()["handle"]

qubitLO.RF_ON()
readoutLO.RF_ON()

# qubitLO.set_freq(4)
# readoutLO.set_freq(7.2578)
attn.setAttenuationStep(devID=1, step=0.1)
attenuation = 22
attn.setAttenuation(devID=1, atten=attenuation)

'''Initialize connection with Zurich Instruments'''
daq, device_qa = create_api_sessions_uhf('dev2528', use_discovery= 1, ip='127.0.0.1')
awg, device_awg = create_api_sessions_hd('dev8233', use_discovery= 1, ip = '127.0.0.1')

