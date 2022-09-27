# -*- coding: utf-8 -*-

from configuration import *
'''---------------------------------------------------------T1---------------------------------------------------------'''

list_of_files = glob.glob('E:\generalized-markovian-noise\%s\T1\*.csv'%(meas_device))
latest_file = max(list_of_files, key=os.path.getctime)
iteration_T1 = int(latest_file[-7:-4].lstrip('0')) + 1

options_T1 = {
    'sampling_rate':    1.2e9,
    'nAverages':        128,
    'Tmax':             100e-6,
    'stepSize':         1000e-9,
    'integration_length': 2.3e-6,
    'cav_resp_time':    0.5e-6,
    'amplitude_hd':     A_d,
    'sequence':         'T1',
    'measPeriod':       400e-6,
    'qubitDriveFreq':   options_rabi['qubitDriveFreq']+detun,
    'sweep':            0,
    'pi2Width':         1/2*pi_pulse*1e-9,
    'AC_pars':          options_rabi['AC_pars'],
    'RT_pars':          [0,0,0],
    'rr_IF':            30e6
}

qubitLO.set_freq(options_T1['qubitDriveFreq']/1e9)

t,I,Q,nPoints = expf.pulse(daq,awg,setup=[0,0,0],**options_T1)
# plot data
data = I
fitted_pars,error = pf.fit_data(x_vector=t,y_vector=data,dt=t[-1]/nPoints,**options_T1)
pf.plot_data(awg,x_vector=t,y_vector=data,fitted_pars=fitted_pars,**options_T1,iteration=iteration_T1,plot_mode=0)


# save data
with open("E:\\generalized-markovian-noise\\%s\\T1\\%s_data_%03d.csv"%(meas_device,'T1',iteration_T1),"w",newline="") as datafile:
    writer = csv.writer(datafile)
    writer.writerow(options_T1.keys())
    writer.writerow(options_T1.values())
    writer.writerow(t)
    writer.writerow(I)
    writer.writerow(Q)

iteration_T1 += 1