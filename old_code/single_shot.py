# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 15:29:47 2022

@author: lfl
"""

options_single_shot = {
    'nAverages':        2**10,
    'setup':            0,
    'pi2Width':         1/2*pi_pulse*1e-9,
    'measPeriod':       600e-6,
    'qubit_drive_amp':     A_d,
    'cav_resp_time':        0.25e-6,
    'integration_length':   2.3e-6,
    'AC_pars':              [options_rabi['AC_pars'][0],0],
    'rr_IF':            30e6
    }

data_OFF, data_pi = expf.single_shot(daq,awg,**options_single_shot)

#make 2D histogram
pf.plot_single_shot(data_OFF, data_pi)