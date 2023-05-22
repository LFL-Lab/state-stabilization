# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 13:51:23 2022

@author: lfl
"""
from mixer_opt import min_leak
# from non_markovian_noise import options_rabi


#qubit
fq = 3.1e9
lo_isol(sa,awg,fc=qubFreq,mixer='qubit',threshold=-70,amp=0.24,calib=0,plot=1)
min_leak(sa,awg,device='dev8233',mode='coarse',mixer='qubit',threshold=-40,f_LO=fq,amp=0.2,plot=True,measON=True)
min_leak(sa,awg,device='dev8233',mode='fine',mixer='qubit',threshold=-70,f_LO=fq,amp=0.2,plot=True,measON=True)
# mixer_calib(sa,awg,'fine',mixer='qubit',threshold=-75,fc=options_ramsey['qubitDriveFreq'],amp=A_d,sweep=0,config=config)

#ac stark
lo_isol(sa,awg,fc=options_rabi['AC_freq'],mixer='ac',threshold=-10,amp=0.3,calib=0,plot=1)
mixer_calib(sa,awg,'coarse',mixer='ac',threshold=-20,fc=options_rabi['AC_freq'],amp=options_rabi['mu'],sweep=0)
mixer_calib(sa,awg,'fine',mixer='ac',threshold=-75,fc=options_rabi['AC_freq'],amp=options_rabi['mu'],sweep=0)

#readout
rrFreq = 7.2578e9
attenuation = 35
lo_isol(sa,daq,fc=rrFreq,mixer='readout',threshold=-70,amp=0.7,calib=0,plot=1)

attn.setAttenuation(devID = 1, atten = 0)
min_leak(sa,daq,device='dev2528',mode='coarse',mixer='readout',threshold=-40,f_LO=rrFreq,amp=0.7,plot=True,measON=True)
attn.setAttenuation(devID = 1, atten = attenuation)

attn.setAttenuation(devID = 1, atten = 0)
min_leak(sa,daq,device='dev2528',mode='fine',mixer='readout',threshold=-70,f_LO=rrFreq,amp=0.7,plot=True,measON=True)
attn.setAttenuation(devID = 1, atten = attenuation)
