#!/usr/bin/env python
# coding: utf-8

# ## Demo introduction and preparition for single qubits characterization

# This demo file is for basic qubit charaterization. The including exerimental sequences are
# - Resonator spectroscopy
# - Qubit spectroscopy
# - Rabi flopping
# - T1 and T2
# - Mixer calibration
#
# Required instruments
# - HDAWG and UHFQA
# - LO MW sources for both qubit driving and qubit readout
# - Frequency up and down conversion units
# - Control PC, USB or LAN connection
#
# Recommanded connections for HDAWG and UHFQA
# - Ref clc out of UHFQA connected to Ref clc in of HDAWG
# - Trigger output 1 of HDAWG to Ref/Trigger 1 of UHFQA
# - Enable DIO connection if it is necessary
#

# ## Import Modules

# In[4]:


#  ---------- keyboard shortcuts in 'Help' tap ------------
# restart kernel when somthing changed in subfunctions

# from VISAdrivers import sa_api as sa
import time
import importlib.util
import json
import sys, os
import LO845m as LO
import numpy as np
import UHFQA as qa
import HDAWG as hd
import experiment_funcs as expf
import matplotlib.pyplot as plt
import csv
import glob
import scipy as scy
import plot_functions as pf
import h5py
from VISAdrivers.LabBrick_LMS_Wrapper import LabBrick_Synthesizer
pi = np.pi
# from VISAdrivers.vaunix_attenuator_wrapper import VaunixAttenuator

'''Instruments and connection'''
qa_id = 'dev2528'
awg_id = 'dev8233'
meas_device = "CandleQubit_5"

# Instrument Addresses
qubitLO_IP = "USB0::0x03EB::0xAFFF::621-03A100000-0538::0::INSTR"
readoutLO_IP = "USB0::0x03EB::0xAFFF::621-03A100000-0519::0::INSTR"
acStarkLO = 21841
# readout_attn = 26551

# initialize instruments as python objects
qubitLO = LO.LO(address=qubitLO_IP,reset=False)
readoutLO = LO.LO(readoutLO_IP,reset=False)
acStarkLO = LabBrick_Synthesizer()
acStarkLO.initDevice(21841)
# readout_attn = LabBrick_attn()
# readout_attn.initDevice(26551)

qubitLO.RF_ON()
readoutLO.RF_ON()
acStarkLO.setRFOn(bRFOn=True)

qubitLO.set_freq(3.21)
readoutLO.set_freq(7.2275)
acStarkLO.setPowerLevel(9.0)
acStarkLO.getUseInternalRef()
acStarkLO.setFrequency(7.3586e9)

'''Initialize connection with Zurich Instruments'''
daq, device_qa = qa.create_api_sessions_uhf('dev2528', use_discovery= 1, ip='127.0.0.1')
awg, device_awg = hd.create_api_sessions_hd('dev8233', use_discovery= 1, ip = '127.0.0.1')

# set clocks to 10 MHz reference
daq.setInt('/dev2528/system/extclk', 1)
awg.setInt('/dev8233/system/clocks/referenceclock/source', 1)

'''Channel offsets'''
# read the current channel offset values and store them for future reference in case the QA needs power cycling
# qubit mixer
offset_qubit_ch1 = awg.get('/dev8233/sigouts/0/offset')['dev8233']['sigouts']['0']['offset']['value']
offset_qubit_ch2 = awg.get('/dev8233/sigouts/2/offset')['dev8233']['sigouts']['2']['offset']['value']

# AC stark mixer
offset_ac_stark_ch1 = awg.get('/dev8233/sigouts/1/offset')['dev8233']['sigouts']['1']['offset']['value']
offset_ac_stark_ch2 = awg.get('/dev8233/sigouts/3/offset')['dev8233']['sigouts']['3']['offset']['value']

# readout mixer
offset_readout_ch1 = daq.get('/dev2528/sigouts/0/offset')['dev2528']['sigouts']['0']['offset']['value']
offset_readout_ch2 =  daq.get('/dev2528/sigouts/1/offset')['dev2528']['sigouts']['1']['offset']['value']

awg.setDouble('/dev8233/sigouts/0/offset',offset_qubit_ch1)
awg.setDouble('/dev8233/sigouts/2/offset',offset_qubit_ch2)
awg.setDouble('/dev8233/sigouts/1/offset',offset_ac_stark_ch1)
awg.setDouble('/dev8233/sigouts/3/offset',offset_ac_stark_ch2)

daq.setDouble('/dev2528/sigouts/0/offset',offset_readout_ch1)
daq.setDouble('/dev2528/sigouts/1/offset',offset_readout_ch2)










'''-----------------------------------------------------Sweep AC Stark Amplitude------------------------------------------------------'''
iteration_ramsey_sweep = 1

options_ramsey_sweep = {
    'nAverages':        256,
    'Tmax':             20e-6,
    'stepSize':         20e-9,
    'integration_length': 2e-6,
    'amplitude_hd':     1.0,
    'sequence':         'ramsey',
    'channel':          0,
    'measPeriod':       200e-6,
    'qubitDriveFreq':   3.877e9,
    'sweep':            0,
    'pi2Width':         1/2*pi_pulse*1e-9,
    'AC_pars':          [0.4,0],
    'RT_pars':          [0,0]
}


amp = np.linspace(0.3,0.6,20)
detun_arr = np.zeros(len(amp))
exp_pars = [options_ramsey_sweep['amplitude_hd'],options_ramsey_sweep['qubitDriveFreq'],options_ramsey_sweep['AC_pars']]

with open("E:\generalized-markovian-noise\%s_data_%03d.csv"%('ramsey',iteration_ramsey_sweep),"w",newline="") as datafile:
    writer = csv.writer(datafile)
    writer.writerow(exp_pars)


for i in range(len(amp)):
    options_ramsey_sweep['AC_pars'] = [amp[i],0]
    t,ch1Data,ch2Data,nPoints = expf.pulse(daq,awg,qubitLO,**options_ramsey_sweep)
    # plot data
    detun_arr[i],T_phi,error = pf.pulse_plot1d(sequence='ramsey',x_vector=t,y_vector=ch1Data.imag,plot=1,dt=options_ramsey_sweep['Tmax']*1e6/nPoints,qubitDriveFreq=options_ramsey_sweep['qubitDriveFreq'],amplitude_hd=options_ramsey_sweep['amplitude_hd'],fitting=1,AC_pars=options_ramsey_sweep['AC_pars'],RT_pars=options_ramsey_sweep['RT_pars'],pi2Width=options_ramsey_sweep['pi2Width'],iteration=iteration_ramsey_sweep)

    if max(error) > 1:
        print('Low Quality Data or Bad Fit')
        i -= 1

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(amp,detun_arr, '-o', markersize = 3, c='C0')
ax.set_ylabel('$\Delta$ (MHz)')
ax.set_xlabel('AC stark amplitude (Volts)')
ax.set_title('Ramsey AC stark sweep %03d'%(iteration_ramsey_sweep))
# save data
exp_pars = [options_ramsey_sweep['amplitude_hd'],options_ramsey_sweep['qubitDriveFreq'],detun_arr,T_phi,options_ramsey_sweep['AC_pars']]
# DataPath = 'E:\generalized-markovian-noise\\{:04}\\{:02}\\Data_{:02}{:02}\\'.format(now.year,now.month,now.month,now.day)
with open("E:\generalized-markovian-noise\%s_data_%03d.csv"%('ramsey',iteration_ramsey),"w",newline="") as datafile:
    writer = csv.writer(datafile)
    writer.writerow(exp_pars)


iteration_ramsey_sweep = 1



'''------------------------------------------------------Echo Statistics---------------------------------------------------'''
'''DESCRIPTION: Repeat Echo Measurement for a couple of hours to determine timescale of environmental fluctuations'''
iteration_echo_statistics = 11
exp = 'echo statistics'
options_echo_statistics = {
  'sampling_rate':      2.4e9,
    'nAverages':        512,
    'Tmax':             options_echo['Tmax'],
    'stepSize':         options_echo['stepSize'],
    'integration_length': 2.0e-6,
    'cav_resp_time':    options_echo['cav_resp_time'],
    'amplitude_hd':     A_d,
    'sequence':         'echo',
    'nSteps':           nPoints,
    'measPeriod':       200e-6,
    'qubitDriveFreq':   options_rabi['qubitDriveFreq']+detun,
    'sweep':            0,
    'pi2Width':         options_echo['pi2Width'],
    'pi2Width_Y':       28.5e-9,
    'pipulse_position': 117e-9,
    'AC_pars':          [options_rabi['AC_pars'][0],0.05],
    'RT_pars':          [0,0],
}

nReps = 250
rep = np.arange(nReps)
T2_arr = np.zeros(nReps)
error_arr = np.zeros(nReps)
now = np.zeros(nReps)
plot = 1

start = time.time()
for i in range(nReps):
    I,Q,nPoints = expf.pulse(daq,awg,qubitLO,setup=[1,1,1],**options_echo_statistics)
    T2_arr[i],error = pf.pulse_plot1d(sequence='echo',x_vector=t, y_vector=I,plot=plot,dt=options_echo_statistics['Tmax']*1e6/nPoints,qubitDriveFreq=options_echo_statistics['qubitDriveFreq'],amplitude_hd=options_echo_statistics['amplitude_hd'],fitting=1,AC_pars=options_echo_statistics['AC_pars'],RT_pars=options_echo_statistics['RT_pars'],pi2Width=options_echo_statistics['pi2Width'],iteration=iteration_echo_statistics)
    error_arr[i] = max(error)
    now[i] = time.time()  - start

# # discard values with high errors
T2_arr_clean = np.zeros(nReps)
rep_clean = np.zeros(nReps)
now_clean = np.zeros(nReps)

def condition(x): return x > 5

bad_index_arr = [idx for idx, element in enumerate(error_arr) if condition(element)]

T2_arr_clean = np.delete(T2_arr,bad_index_arr)
rep_clean = np.delete(rep,bad_index_arr)
# time_arr = 14.5*np.arange(len(rep_clean))
time_arr = 14.5*np.arange(len(now))

fig, ax= plt.subplots()
fig.suptitle('Echo Statistics %03d'%(iteration_echo_statistics))
ax.plot(time_arr,T2_arr, '-o', markersize = 3, c='C0')
ax.set_xlabel('time (sec)')
ax.set_ylabel('$T_2 (\mu s)$')
ax.set_ylim((0,2))

# save data
exp_pars = [options_echo_statistics['amplitude_hd'],options_echo_statistics['qubitDriveFreq'],options_echo_statistics['AC_pars']]
with open("E:\\generalized-markovian-noise\\echo\\echo_statistics\\%s_data_%03d.csv"%('echo_statistics',iteration_echo_statistics),"w",newline="") as datafile:
    writer = csv.writer(datafile)
    writer.writerow(exp_pars)
    writer.writerow(T2_arr_clean)
    writer.writerow(now_clean)
    writer.writerow(error_arr)

iteration_echo_statistics += 1

'''------------------------------------------------------T1 Statistics---------------------------------------------------'''
'''DESCRIPTION: Repeat Echo Measurement for a couple of hours to determine timescale of environmental fluctuations'''
iteration_T1_statistics = 1
exp = 'T1 statistics'
options_T1_statistics = {
    'nAverages':        128,
    'Tmax':             40e-6,
    'stepSize':         500e-9,
    'nSteps':           81,
    'amplitude_hd':     1.0,
    'sequence':         'T1',
    'channel':          0,
    'measPeriod':       200e-6,
    'qubitDriveFreq':   3.879e9,
    'sweep':            0,
    'pi2Width':         1/2*pi_pulse*1e-9,
    'AC_pars':          [0,0],
    'RT_pars':          [0,0]
}

nReps = 50
rep = np.arange(nReps)
T1_arr = np.zeros(nReps)
error_arr = np.zeros(nReps)
now = np.zeros(nReps)
plot = 1

start = time.time()
for i in range(nReps):
    t,ch1Data,ch2Data,nPoints = expf.pulse(daq,awg,qubitLO,setup=[1,1,1],**options_T1_statistics)
    T1_arr[i],error = pf.pulse_plot1d(sequence='T1',plot=plot,dt=options_T1_statistics['Tmax']*1e6/nPoints,qubitDriveFreq=options_T1_statistics['qubitDriveFreq'],amplitude_hd=options_T1_statistics['amplitude_hd'],x_vector=t, y_vector=ch1Data,fitting=1,AC_pars=options_T1_statistics['AC_pars'],RT_pars=options_T1_statistics['RT_pars'],pi2Width=options_T1_statistics['pi2Width'],iteration=iteration_T1_statistics)
    error_arr[i] = max(error)
    now[i] = time.time()  - start

# # discard values with high errors
T1_arr_clean = np.zeros(nReps)
rep_clean = np.zeros(nReps)
now_clean = np.zeros(nReps)

def condition(x): return x > 5

bad_index_arr = [idx for idx, element in enumerate(error_arr) if condition(element)]

T1_arr_clean = np.delete(T2_arr,bad_index_arr)
rep_clean = np.delete(rep,bad_index_arr)
# time_arr = 14.5*np.arange(len(rep_clean))
time_arr = 14.5*np.arange(len(now))

fig, ax= plt.subplots()
fig.suptitle('T1 Statistics %03d'%(iteration_T1_statistics))
ax.plot(time_arr,T1_arr, '-o', markersize = 3, c='C0')
ax.set_xlabel('time (sec)')
ax.set_ylabel('$T_1 (\mu s)$')

# save data
exp_pars = [options_T1_statistics['amplitude_hd'],options_T1_statistics['qubitDriveFreq'],options_T1_statistics['AC_pars']]
with open("E:\\generalized-markovian-noise\\T1\\T1_statistics\\%s_data_%03d.csv"%('T1_statistics',iteration_T1_statistics),"w",newline="") as datafile:
    writer = csv.writer(datafile)
    writer.writerow(exp_pars)
    writer.writerow(detun_arr)
    writer.writerow(T_phi_arr)
    writer.writerow(now)
    writer.writerow(error_arr)

iteration_T1_statistics += 1




# average data
T2_avg = np.mean(T2_arr,2)
freq_avg = np.mean(ram_freq_arr,2)
Tb_avg = np.mean(T_b,2)

pf.sweep_plot(tau,1e3*B0,T2_avg)

# plot data
fig, (ax1,ax2) = plt.subplots(2,sharex=True)
ax1.plot(rep,detun_arr, '-o', markersize = 3, c='C0')
ax1.set_ylabel('$\Delta$ (MHz)')
fig.suptitle('Ramsey AC stark sweep %03d'%(iteration_ramsey_statistics))
ax2.plot(rep,T_2_arr, '-o', markersize = 3, c='C0')
ax2.set_xlabel('iteration')
ax2.set_ylabel('$T_{\phi} (\mu s)$')


# save data
exp_pars_ramsey_sweep = [options_ramsey['amplitude_hd'],options_ramsey['qubitDriveFreq'],options_ramsey['AC_pars']]

file = h5py.File("E:\generalized-markovian-noise\%s_data_%03d.h5", 'w')
pars = file.create_group("Experimental Parameters")
data = file.create_group("data")

pars.create_dataset("Experimental Parameters",data=optionsRamsey_par_sweep)
sweep_pars.create_dataset("Tau",data=tau)
sweep_pars.create_dataset("Tau",data=B0)
data.create_dataset("Background Data",data=b_data1)
data.create_dataset("Data",data=data1)
file.close()


'''---------------------------------------Echo Parameter Sweep for Random Telegraph Noise-------------------------------------------'''

'''DESCRIPTION:
    1. Generate instances of telegraph noise for all parameter points and store them in separate folders
    2. Repeat echo measurement with different noise realizations
    3. Every 10 or so instances, do echo without RTN to get statistics for the bare T2
    4. Repeat for every parameter point
'''

sweep_count = 30
numRealizations = 100
b_measurements = 4
numIterations = numRealizations + b_measurements
interval = int(numRealizations/b_measurements) # how often to take background data without generalized markovian noise
B0 = np.linspace(0.01,0.2,10) #the relationship between B0 and frequency of oscillations is Omega_R = 25 MHz * A_q
# tau = np.concatenate((np.linspace(0.01,2,8),np.linspace(3,100,2)))
tau = np.linspace(0.01,0.1,3)
# tau = [0.1,0.7,3]
# B0 = [0.07,0.2]

T2_b = np.zeros((len(B0),len(tau),int(numRealizations/interval)),dtype=float)
b_data1 = np.zeros((len(B0),len(tau),int(numRealizations/interval),112),dtype=float)
b_data2  = np.zeros((len(B0),len(tau),int(numRealizations/interval),112),dtype=float)
error_b = np.zeros((len(B0),len(tau),int(numRealizations/interval)),dtype=float)

T2_arr = np.zeros((len(B0),len(tau),numRealizations),dtype=float)
data1 = np.zeros((len(B0),len(tau),numRealizations,148),dtype=float)
data2 = np.zeros((len(B0),len(tau),numRealizations,148),dtype=float)
error_arr = np.zeros((len(B0),len(tau),numRealizations),dtype=float)
# generate noise instances
Tmax = 4e-6
nPoints = expf.roundToBase(2.4e9*Tmax,16)
str_gen_noise = time.time()
expf.gen_noise_realizations(noiseType='RTN',par1_arr=tau,par2_arr=[0],numRealizations=numRealizations,nPoints=nPoints,T_max=Tmax,sweep_count=sweep_count)
end_gen_noise = time.time()
print('Generating noise realizations took: %.1f s' %(end_gen_noise-str_gen_noise))
sweep_name = 'sweep_%03d'%(sweep_count)
parent_dir = 'E:\\generalized-markovian-noise\\sweep_data\\'
path = os.path.join(parent_dir,sweep_name)
os.mkdir(path)
plt.close('all')

# sweep echo
optionsEcho_par_sweep = {
    'sampling_rate':    2.4e9,
    'nAverages':        512,
    'Tmax':             4e-6,
    'stepSize':         30e-9,
    'integration_length': 2.0e-6,
    'cav_resp_time':    0.5e-6,
    'amplitude_hd':     A_d,
    'sequence':         'echo',
    # 'nSteps':           57,
    'measPeriod':       200e-6,
    'qubitDriveFreq':   options_rabi['qubitDriveFreq']+detun,
    'sweep':            1,
    'pi2Width':         options_echo['pi2Width'],
    'piWidth_Y':       options_echo['piWidth_Y'],
    'pipulse_position': options_echo['pipulse_position'],
    'AC_pars':          options_echo['AC_pars'],
    'RT_pars':          [0,0],
    }

plot = 0
header = ['Amplitude (V)','Drive Freq (Hz)','AC Stark Amplitude (V)','AC stark Noise Amplitude (V)','Pi2 Pulse Width (ns)','Y Pi Pulse Width (ns)','Pi Pulse Position (ns)']
exp_pars = [optionsEcho_par_sweep['amplitude_hd'],optionsEcho_par_sweep['qubitDriveFreq'],optionsEcho_par_sweep['AC_pars'][0],optionsEcho_par_sweep['AC_pars'][1],optionsEcho_par_sweep['pi2Width'],optionsEcho_par_sweep['piWidth_Y'],optionsEcho_par_sweep['pipulse_position']]
start_sweep = time.time()
# B0 = [0.179]
# generate data
for i in range(len(B0)):
    optionsEcho_par_sweep['RT_pars'][0] = B0[i]
    for j in range(len(tau)):
        optionsEcho_par_sweep['RT_pars'][1] = tau[j]
        a = 0 # keeps track of background measurements
        b = 0 # keeps track of noisy measurements
        k = 0
        print('Next parameter point: B_0 = %.2f V and tau = %.2f microseconds' %(B0[i],tau[j]))
        while k < numIterations:
            if k % (interval + 1) == 0 and a != b_measurements:
                # get background T2* every 10 or so measurements
                optionsEcho_par_sweep['RT_pars'] = [0,0]
                optionsEcho_par_sweep['nAverages'] = 1024
                optionsEcho_par_sweep['Tmax'] = 3e-6
                optionsEcho_par_sweep['stepSize'] = 30e-9
                # optionsEcho_par_sweep['nAverages'] = 1 #pars used for testing
                # optionsEcho_par_sweep['Tmax'] = 3e-6
                # optionsEcho_par_sweep['stepSize'] = 1000e-9
                print('----------------------------------\nExecuting background Echo measurement')
                t2,I,Q,nPoints = expf.pulse(daq,awg,qubitLO,setup=[0,1,0],sweep_name=sweep_name,**optionsEcho_par_sweep)
                b_data1[i,j,a,:] = I
                b_data2[i,j,a,:] = Q
                T2_b[i,j,a],error = pf.pulse_plot1d(sequence='echo',x_vector=t2, y_vector=I,plot=plot,dt=optionsEcho_par_sweep['Tmax']*1e6/nPoints,qubitDriveFreq=optionsEcho_par_sweep['qubitDriveFreq'],amplitude_hd=optionsEcho_par_sweep['amplitude_hd'],fitting=1,AC_pars=optionsEcho_par_sweep['AC_pars'],RT_pars=optionsEcho_par_sweep['RT_pars'],pi2Width=optionsEcho_par_sweep['pi2Width'])
                error_b[i,j,a] = max(error)
                a += 1
                print('End measurement\n----------------------------------' )
            else:
                if k % (interval + 1) == 1:
                    setup = [0,1,0]
                else:
                    setup = [2,1,1]

                # if k % 10 == 0:
                #     plot = 1
                # else:
                #     plot = 0

                optionsEcho_par_sweep['RT_pars'] = [B0[i],tau[j]]
                optionsEcho_par_sweep['Tmax'] = 4e-6
                optionsEcho_par_sweep['nAverages'] = 512
                optionsEcho_par_sweep['stepSize'] = 30e-9
                # optionsEcho_par_sweep['Tmax'] = 4e-6 #pars used for testing
                # optionsEcho_par_sweep['nAverages'] = 1
                # optionsEcho_par_sweep['stepSize'] = 30e-9
                print('----------------------------------\nStart %s measurement' %("Echo"))
                print('Implementing noise realization %d' %(b+1))
                t1,I,Q,nPoints = expf.pulse(daq,awg,qubitLO,setup=setup,sweep_name=sweep_name,instance=b,**optionsEcho_par_sweep)
                data1[i,j,b,0:nPoints] = I
                data2[i,j,b,0:nPoints] = Q
                T2_arr[i,j,b],error = pf.pulse_plot1d(sequence='echo',x_vector=t1, y_vector=I,plot=plot,dt=optionsEcho_par_sweep['Tmax']*1e6/nPoints,qubitDriveFreq=optionsEcho_par_sweep['qubitDriveFreq'],amplitude_hd=optionsEcho_par_sweep['amplitude_hd'],fitting=1,AC_pars=optionsEcho_par_sweep['AC_pars'],RT_pars=optionsEcho_par_sweep['RT_pars'],pi2Width=optionsEcho_par_sweep['pi2Width'])
                error_arr[i,j,b] = max(error)
                b += 1
                print('End measurement\n----------------------------------' )
            k += 1

        # save data after each parameter sweep point
        filename = 'RTN_B0_%d_mV_tau_%d_ns' %(round(B0[i]*1e3),round(tau[j]*1e3))
        with open("E:\generalized-markovian-noise\sweep_data\%s\data_%s.csv"%(sweep_name,filename),"w",newline="") as datafile:
        # with open("E:\\generalized-markovian-noise\\sweep_data\\%s\\redo\\data_%s.csv"%(sweep_name,filename),"w",newline="") as datafile:
            writer = csv.writer(datafile)
            writer.writerow(['Background Time Data'])
            writer.writerow(t2)
            writer.writerow(['Time Data'])
            writer.writerow(t1)
            writer.writerow(['Background Data: Channel 1'])
            writer.writerows(b_data1[i,j,:,:])
            writer.writerow(['Background Data: Channel 2'])
            writer.writerows(b_data2[i,j,:,:])
            writer.writerow(['Data: Channel 1'])
            writer.writerows(data1[i,j,:,:])
            writer.writerow(['Data: Channel 2'])
            writer.writerows(data2[i,j,:,:])

end_sweep = time.time()
print('Total Sweep Duration: %.1f s or %.1f hours, or %.1f days' %(end_sweep-start_sweep,(end_sweep-start_sweep)/3600,(end_sweep-start_sweep)/(3600*24)))
sweep_count += 1