# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 12:32:09 2021

@author: lfl
"""

'''---------------------------------------------------------Do Ramsey and Sweep pi2 width---------------------------------------------------------'''
iteration_ramsey_sweep_pi2_width = 1
exp = 'ramsey sweep pi2 width'
options_ramsey_sweep_pi2_width = {
    'nAverages':        128,
    'Tmax':             5e-6,
    'stepSize':         10e-9,
    'amplitude_hd':     1.0,
    'sequence':         'ramsey',
    'channel':          0,
    'measPeriod':       200e-6,
    'qubitDriveFreq':   3.8675e9,
    'AC_pars':          [0.4,0],
    'RT_pars':          [0,0],
    'pi2Width':         50
}


pi2Width = 1e-9*np.linspace(20,60,21)

data = np.zeros((375,len(pi2Width)),dtype=complex)

# exp_pars = [options_ramsey_sweep_pi2_width['amplitude_hd'],options_ramsey_sweep_pi2_width['qubitDriveFreq'],options_ramsey_sweep_pi2_width['AC_pars']]
# with open("E:\generalized-markovian-noise\%s_data_%03d.csv"%('ramsey',iteration_ramsey_sweep_pi2_width),"w",newline="") as datafile:
#     writer = csv.writer(datafile)
#     writer.writerow(exp)
#     writer.writerow(exp_pars)

    
for i in range(len(pi2Width)):  
    options_ramsey_sweep_pi2_width['pi2Width'] = pi2Width[i]
    t,ch1Data,ch2Data,nPoints = expf.pulse(daq,awg,qubitLO,setup=[0,0,0],**options_ramsey_sweep_pi2_width)    
    data[:,i] = ch1Data
    # writer.writerow(data[:,i])
    detuning,T_phi,error = pf.pulse_plot1d(sequence='ramsey',dt=options_ramsey_sweep_pi2_width['Tmax']*1e6/nPoints,qubitDriveFreq=options_ramsey_sweep_pi2_width['qubitDriveFreq'],amplitude_hd=options_ramsey_sweep_pi2_width['amplitude_hd'],x_vector=t, y_vector=ch1Data,fitting=1,AC_pars=options_ramsey_sweep_pi2_width['AC_pars'],RT_pars=options_ramsey_sweep_pi2_width['RT_pars'],pi2Width=options_ramsey_sweep_pi2_width['pi2Width'],iteration=iteration_ramsey_sweep_pi2_width)
    
# save data
exp_pars = [options_ramsey1['amplitude_hd'],options_ramsey2['qubitDriveFreq'],detun_arr,T_phi,options_ramsey2['AC_pars']]
# DataPath = 'E:\generalized-markovian-noise\\{:04}\\{:02}\\Data_{:02}{:02}\\'.format(now.year,now.month,now.month,now.day)
with open("E:\generalized-markovian-noise\%s_data_%03d.csv"%('ramsey',iteration_ramsey),"w",newline="") as datafile:
    writer = csv.writer(datafile)
    writer.writerow(exp_pars)

 
iteration_ramsey_sweep = 1