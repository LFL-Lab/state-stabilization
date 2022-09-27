# -*- coding: utf-8 -*-

from configuration import *

'''---------------------------------------------------------Echo---------------------------------------------------------'''

list_of_files = glob.glob('E:\generalized-markovian-noise\%s\Echo\*.csv'%(meas_device))
latest_file = max(list_of_files, key=os.path.getctime)
iteration_echo = int(latest_file[-7:-4].lstrip('0')) + 1

options_echo = {
    'sampling_rate':    1.2e9,
    'nAverages':        256,
    'Tmax':             30e-6,
    'stepSize':         100e-9,
    'integration_length': 2.3e-6,
    'prePulseLength':   1500e-9,
    'postPulseLength':  200e-9,
    'cav_resp_time':    options_rabi['cav_resp_time'],
    'amplitude_hd':     A_d,
    'sequence':         'echo',
    # 'nSteps':           57,
    'measPeriod':       300e-6,
    'qubitDriveFreq':   options_rabi['qubitDriveFreq']+detun,
    'sweep':            0,
    'pi2Width':         1/2*pi_pulse*1e-9,
    # 'piWidth_Y':        44e-9,
    # 'pipulse_position': 104e-9,
    'AC_pars':          [options_rabi['AC_pars'][0],0],
    'AC_freq':          options_rabi['AC_freq'],
    'RT_pars':          [0,0,0],
    'rr_IF':            30e6
}

qubitLO.set_freq(options_echo['qubitDriveFreq']/1e9)

t,I,Q,nPoints = expf.pulse(daq,awg,setup=[0,0,0],**options_echo)
# plot data
data = I
fitted_pars,error = pf.fit_data(x_vector=t,y_vector=data,dt=t[-1]/nPoints,**options_echo)
pf.plot_data(awg,x_vector=t,y_vector=data,fitted_pars=fitted_pars,**options_echo,iteration=iteration_echo,plot_mode=0)

# save data
with open("E:\\generalized-markovian-noise\\%s\\echo\\%s_data_%03d.csv"%(meas_device,'echo',iteration_echo),"w",newline="") as datafile:
    writer = csv.writer(datafile)
    writer.writerow(options_echo.keys())
    writer.writerow(options_echo.values())
    writer.writerow(t)
    writer.writerow(I)
    writer.writerow(Q)

iteration_echo += 1