# -*- coding: utf-8 -*-

'''---------------------------------------------------------Ramsey---------------------------------------------------------'''
from configuration import *

try:
    directory = 'E:\\generalized-markovian-noise\\%s\\sweep_data\\ramsey\\' %(meas_device)
    latest_file = max(glob.glob(os.path.join(directory, '*')), key=os.path.getmtime)
    iteration_ramsey = int(latest_file[-3:].lstrip('0')) + 1
except:
    iteration_ramsey = 1

detun = 0e6

options_ramsey = {
    'sampling_rate':    1.2e9,
    'nAverages':        256,
    'Tmax':             50e-6,
    'stepSize':         200e-9,
    'prePulseLength':   1500e-9,
    'postPulseLength':  200e-9,
    'integration_length':   2.3e-6,
    'cav_resp_time':    options_rabi['cav_resp_time'],
    'amplitude_hd':     A_d,
    'active_reset':     False,
    'threshold':        threshold,
    'sequence':         'ramsey',
    'measPeriod':       500e-6,
    'qubitDriveFreq':   options_rabi['qubitDriveFreq']+detun,
    'sweep':            0,
    'pi2Width':         1/2*pi_pulse*1e-9,
    'AC_pars':          [0.0,0],
    'AC_freq':          options_rabi['AC_freq'],
    'RT_pars':          [0,0,0],
    'rr_IF':            30e6
    }

qubitLO.set_freq(options_ramsey['qubitDriveFreq']/1e9)

t,I,Q,nPoints = expf.pulse(daq,awg,setup=[0,0,0],**options_ramsey)

# plot data
data = I
fitted_pars,error = pf.fit_data(x_vector=t,y_vector=data,dt=t[-1]/nPoints,**options_ramsey)
pf.plot_data(awg,x_vector=t,y_vector=data,fitted_pars=fitted_pars,**options_ramsey,iteration=iteration_ramsey,plot_mode=0)

# save data
with open("E:\\generalized-markovian-noise\\%s\\ramsey\\%s_data_%03d.csv"%(meas_device,'ramsey',iteration_ramsey),"w",newline="") as datafile:
    writer = csv.writer(datafile)
    writer.writerow(options_ramsey.keys())
    writer.writerow(options_ramsey.values())
    writer.writerow(t)
    writer.writerow(I)
    writer.writerow(Q)

iteration_ramsey += 1