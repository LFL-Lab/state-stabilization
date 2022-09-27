# -*- coding: utf-8 -*-
'''-----------------------------------------------------spectroscopy------------------------------------------------------'''

from configuration import *



try:
    directory = 'E:\\generalized-markovian-noise\\%s\\sweep_data\\spectroscopy\\' %(meas_device)
    latest_file = max(glob.glob(os.path.join(directory, '*')), key=os.path.getmtime)
    iteration_spec = int(latest_file[-3:].lstrip('0')) + 1
except:
    iteration_spec = 1

options_spec = {
    'frequencies':      np.arange(start=3.0,stop=3.4,step=10e-5), # frequencies are in GHz
    'nAverages':        1024,
    'setup':            0,
    'qubit_drive_amp':     200e-3,
    'readout_drive_amp':     0.7,
    'cav_resp_time':        0.25e-6,
    'integration_length':   2.3e-6,
    'AC_pars':              [0.0,0]
    }

p_data,I,Q = expf.spectroscopy(daq,awg,qubitLO=qubitLO,**options_spec)
# get current QA AWG offset for later reference
# offset_readout_ch1 = daq.get('/dev2528/sigouts/0/offset')['dev2528']['sigouts']['0']['offset']['value']
# # # # # set output of QA AWG to 1 V (what we use for readout)
# daq.setDouble('/dev2528/sigouts/0/offset',offset_readout_ch1+0.7)
# # # # # measure readout ON power
# readout_freq,readout_power = expf.get_power(freq=5.8017e9,plot=1)
# # # # # reset offset value back to original for optimal ON/OFF power ratio
# daq.setDouble('/dev2528/sigouts/0/offset',offset_readout_ch1)
# plot
pf.spec_plot(freq=options_spec['frequencies'],I=I,Q=Q,qubit_drive_amp=options_spec['qubit_drive_amp'])

exp_pars = options_spec
with open("E:\\generalized-markovian-noise\\%s\\spectroscopy\\%s_data_%03d.csv"%(meas_device,'spectroscopy',iteration_spec),"w",newline="") as datafile:
    writer = csv.writer(datafile)
    writer.writerow(options_spec.keys())
    writer.writerow(options_spec.values())
    writer.writerow(I)
    writer.writerow(Q)

iteration_spec += 1


