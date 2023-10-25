# -*- coding: utf-8 -*-
"""
Created on Thu May 18 11:46:18 2023

@author: lfl
"""

atten = np.arange(10,40,2)

qb.exp_pars = {
    'num_samples':          512,
    'n_avg':                1,
    'fsAWG':                2.4e9,
    'qubit_reset_time':     300e-6,
    }


for i in atten:
    qb.update_qb_value('rr_atten',i)
    data = qb.single_shot()
    #make 2D histogram
    qb.plot_single_shot(data)
    
    
pi_amp = np.arange(0.21,0.25,1e-3)

qb.exp_pars = {
    'num_samples':          512,
    'n_avg':                1,
    'fsAWG':                2.4e9,
    'qubit_reset_time':     300e-6,
    }


for i in pi_amp:
    qb.update_qb_value('pi_amp',i)
    data = qb.single_shot()
    #make 2D histogram
    qb.plot_single_shot(data)


#%% Punchout Measurement
'''------------------------------------------------------Punchout----------------------------------------------'''


freqs = np.arange(start=6.703,stop=6.706,step=50e-6) # frequencies are in GHz
atten = np.arange(0,35,2.5)
#print(atten)
p_data = np.zeros((len(atten),len(freqs)))

qb.exp_pars = {
    'exp':                  'spectroscopy',
    'n_avg':                64,
    'element':              'rr',
    'rr_reset_time':        30e-6,
    'satur_dur':            10e-6,
    'rr_atten':             0,
    'on_off':               True,
    }

for i,a in enumerate(atten):
    qb.exp_pars['rr_atten'] = a
    p_data[i,:],I,Q = qb.rr_spectroscopy(freqs=freqs,device_name = device_name,save_data=False)
    pt.rr_spec_plot(freq=freqs,I=I,Q=Q,mag=p_data[i,:],exp_pars=qb.exp_pars,qb_pars=qb.qb_pars,df=1e9*(freqs[1]-freqs[0]),find_pks=True)
    
df = pt.heatplot(xdata=np.around(freqs,6),ydata=atten,z_data = p_data*1e3,xlabel='Frequency (GHz)',
            ylabel='Attenuation (dB)', normalize = True,cbar_label='Magnitude',title='Punchout Measurement')

#%% rabi amplitude sweep
'''------------------------------------------------------Sweep Rabi Amplitude----------------------------------------------'''
detuning = 0
qb.exp_pars = {
    'initial-state':        '0',
    'exp':                  't-rabi',
    'n_avg':                1024,
    'x0':                   13e-9,
    'xmax':                 2e-6,
    'dx':                   6e-9,
    'fsAWG':                2.4e9,
    'amp_q':                0.3,
    'active_reset':         False,
    'qubit_drive_freq':     qb.qb_pars['qb_freq']+detuning,
    'tomographic-axis':     'Z',
    }



amp = np.linspace(0.01,0.501,21)
Omega_Rabi = []
data = np.zeros((len(amp),299))

for i in range(len(amp)):
    qb.exp_pars['amp_q'] = amp[i]
    # generate data
    print(' Qubit Drive Amplitude = %.3f V'%(amp[i]))
    t,data[i,:],nSteps = qb.pulsed_exp(qb=qb_name,verbose=1,device_name=device_name,check_mixers=False,save_data=False)
    # plot data
    fitted_pars,error = pt.fit_data(x_vector=t,y_vector=data[i,:],exp='t-rabi',dx=t[-1]/nSteps*1e6,verbose=0)
    # pt.plot_t_rabi_data(t,data,fitted_pars,qb=qb_name,exp_pars=qb.exp_pars,qb_pars=qb.qb_pars,device_name=device_name,project=project)
    Omega_Rabi.append(1e3/fitted_pars[1])
    # pf.plot_data(awg,x_vector=t,y_vector=I,fitted_pars=fitted_pars,**options_rabi_sweep,plot_mode=0)

# load data from previous measurement
# data = pd.read_csv("E:\\generalized-markovian-noise\\%s\\Rabi\\RabiPowerSweep\\%s_data_%03d.csv"%(meas_device,'RabiPowerSweep',iteration_rabi_power_sweep),nrows=2,on_bad_lines='skip',skiprows=2,header=None).to_numpy(np.float64)
# amp = data[0,:]*1e3
# Omega_Rabi = data[1,:]
start = 1
cutoff = 10
best_vals, covar = scy.optimize.curve_fit(line, amp[start:cutoff],Omega_Rabi[start:cutoff],xtol=1e-6,maxfev=3000)

fig, ax1 = plt.subplots(dpi=300)
# plt.xticks(np.arange(-0.1e3,1.1e3,step=0.2))
# plt.yticks(np.arange(0,12,step=2))
left,bottom,width,height = [0.5, 0.25, 0.3, 0.4]
ax2 = fig.add_axes([left,bottom,width,height])
ax1.plot(amp[start:]*1e3,Omega_Rabi[start:], '-o', markersize = 3, c='C0')
ax1.set_xlabel('Qubit Drive Amplitude (mV)')
ax1.set_ylabel('$\Omega_R/2\pi$ (MHz)')
ax2.plot(amp[start:cutoff]*1e3,Omega_Rabi[start:cutoff],amp[start:cutoff]*1e3,line(amp[start:cutoff],best_vals[0],best_vals[1]))
ax2.set_xlabel('Qubit Drive Amplitude (mV)',fontsize=9)
ax2.set_ylabel('$\Omega_R/2\pi$ (MHz) ',fontsize=10)
# plt.xticks(np.arange(0,amp[35],step=1e-4),fontsize=9)
# plt.yticks(np.arange(0,np.max(Omega_Rabi),step=2),fontsize=9)
inset_box_txt = '$f_R=$'+f"{best_vals[0]:.2e} MHz/V"+'$\\cdot A_d +$' +f"{best_vals[1]:.2e} MHz/V"
plt.gcf().text(0.45, 0.675, inset_box_txt, fontsize=12)
plt.show()



#%%  Sweep sigma (sigma calibration)

try:
    list_of_files = glob.glob('E:\generalized-markovian-noise\%s\Ramsey\AC_stark_calibrations\sigma_calibration_data_*.csv'%(meas_device))
    latest_file = max(list_of_files, key=os.path.getctime)
    iteration_ramsey_sigma_calibration = int(latest_file[-7:-4].lstrip('0')) + 1
except:
    iteration_ramsey_sigma_calibration = 1

options_ramsey_sigma_calibration = {
    'nAverages':        128,
    'Tmax':             60e-6,
    'stepSize':         100e-9,
    'prePulseLength':   options_ramsey['prePulseLength'],
    'postPulseLength':  options_ramsey['postPulseLength'],
    'integration_length': 2.3e-6,
    'cav_resp_time':    options_ramsey['cav_resp_time'],
    'amp_q':     A_d,
    'active_reset':     True,
    'threshold':        threshold,
    'measPeriod':       600e-6,
    'qubitDriveFreq':   options_rabi['qubitDriveFreq']+detun,
    'pi2Width':         1/2*pi_pulse*1e-9,
    'mu':               0,
    'sigma':            0,
    'AC_freq':          options_ramsey['AC_freq'],
    'B0':               0,
    'nu':               0,
    'tauk':             0,
    'source':           2,
}

numRealizations = 30
# sigma = np.linspace(0.0005,0.15,31)
sigma = np.arange(0.0005,0.125,step=5e-3)
T_phi_arr = np.zeros(len(sigma))
error_arr = np.zeros(len(sigma))


for i in range(len(sigma)):
    expf.mixer_calib(sa,awg,'fine',mixer='qubit',threshold=-75,fc=options_ramsey['qubitDriveFreq'],amp=A_d,sweep=1,config=config)
    options_ramsey_sigma_calibration['sigma'] = sigma[i]
    if sigma[i] < 30e-3:
        options_ramsey_sigma_calibration['Tmax'] = 100e-6
        options_ramsey_sigma_calibration['stepSize'] = 400e-9
    elif sigma[i] >= 30e-3 and sigma[i] < 40e-3:
        options_ramsey_sigma_calibration['Tmax'] = 40e-6
        options_ramsey_sigma_calibration['stepSize'] = 100e-9
    elif sigma[i] >= 40e-3 and sigma[i] < 80e-3:
        options_ramsey_sigma_calibration['Tmax'] = 10e-6
        options_ramsey_sigma_calibration['stepSize'] = 60e-9
        numRealizations = 50
    elif sigma[i] >= 80e-3:
        options_ramsey_sigma_calibration['Tmax'] = 5e-6
        options_ramsey_sigma_calibration['stepSize'] = 60e-9
        numRealizations = 100
    nPoints,nSteps, pulse_length_increment, pulse_length_start = expf.calc_nSteps(sequence='ramsey',fsAWG=1.2e9,stepSize=options_ramsey_sigma_calibration['stepSize'],Tmax=options_ramsey_sigma_calibration['Tmax'])
    data_arr = np.zeros(nSteps-1)
    for k in range(numRealizations):
        if k == 0:
            setup = [0,1,0]
        else:
            setup = [2,1,1]
            options_ramsey_sigma_calibration['nSteps'] = nSteps
        try:
            t,data,nSteps = expf.pulse(daq,awg,setup=setup,sequence='ramsey',**options_ramsey_sigma_calibration)
        # pf.plot_data(awg,x_vector=t,y_vector=I,**options_ramsey_sigma_calibration,plot_mode=1)
            data_arr += data
        except Exception as exc:
            print(exc)
            k -= 1
    data_avg = data_arr/numRealizations
    fitted_pars,error = pf.fit_data(x_vector=t,y_vector=data_avg,sequence='echo',dt=t[-1]/nSteps,**options_ramsey_sigma_calibration)
    # pf.plot_data(awg,x_vector=t,y_vector=data_avg,fitted_pars=fitted_pars,**options_ramsey_sigma_calibration,plot_mode=1)
    T_phi_arr[i] = fitted_pars[1]
    error_arr[i] = error[1]

# save data
with open("E:\\generalized-markovian-noise\\%s\\Ramsey\\AC_stark_calibrations\\%s_data_%03d.csv"%(meas_device,'sigma_calibration',iteration_ramsey_sigma_calibration),"w",newline="") as datafile:
    writer = csv.writer(datafile)
    writer.writerow(options_ramsey_sigma_calibration.keys())
    writer.writerow(options_ramsey_sigma_calibration.values())
    writer.writerow(sigma)
    writer.writerow(T_phi_arr)

iteration_ramsey_sigma_calibration += 1

# load data from previous measurement
data = pd.read_csv("E:\\generalized-markovian-noise\\%s\\Ramsey\\AC_stark_calibrations\\%s_data_%03d.csv"%(meas_device,'sigma_calibration',iteration_ramsey_sigma_calibration),nrows=2,on_bad_lines='skip',skiprows=2,header=None).to_numpy(np.float64)
sigma = data[0,:]
T_phi_arr = data[1,:]

fig = plt.figure(figsize=(10,7.5),dpi=600)
ax = fig.add_subplot(111)
ax.errorbar(sigma*1e3,T_phi_arr,yerr=error_arr/np.sqrt(numRealizations),fmt='-o', linewidth=3,markersize = 10, c='r',ecolor='r')
ax.set_ylabel('$T_2^R$ ($\mu$s)',fontsize=30)
ax.set_xlabel('Noise Amplitude (mV)',fontsize=30)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
# ax.set_title('Ramsey AC Stark Sigma Calibration')
ax.tick_params(axis='both',direction='in',bottom=True, top=True, left=True, right=True,size=16)
plt.savefig("G:\\Shared drives\\LFL\Projects\\Generalized Markovian noise\\paper-figures\\sigma_calibration\\"+"sigma_calibration", bbox_inches="tight",
            pad_inches=0.3, transparent=True,format='pdf')



#%% Statistics Measurements

'''------------------------------------------------------Pipulse Statistics---------------------------------------------------'''
'''DESCRIPTION: Repeat Rabi Measurement for a couple of hours to determine timescale of environmental fluctuations'''


options_rabi_statistics = {
    'sampling_rate':    options_rabi['sampling_rate'],
    'nAverages':        options_rabi['nAverages'],
    'Tmax':             options_rabi['Tmax'],
    'stepSize':         options_rabi['stepSize'],
    'integration_length': options_rabi['integration_length'],
    'cav_resp_time':    options_rabi['cav_resp_time'],
    'amp_q':     A_d,
    'nSteps':           nSteps,
    'measPeriod':       options_rabi['measPeriod'],
    'qubitDriveFreq':   options_rabi['qubitDriveFreq'],
    'mu':               0.305,
    'sigma':            0,
    'B0':               0,
    'nu':               0,
    'tauk':             0,
    'source':           2
}

nReps = 200
rep = np.arange(nReps)
rabi_freq_arr = np.zeros(nReps)
error_arr = np.zeros(nReps)
time_arr = np.zeros(nReps)
plot = 1

start = time.time()
for i in range(nReps):
    t,data,nSteps = expf.pulse(daq,awg,setup=[1,1,1],sequence='rabi',**options_rabi_statistics)
    fitted_pars,error = pf.fit_data(x_vector=t,y_vector=data,sequence='rabi',dt=t[-1]/nSteps,**options_rabi_statistics)
    pf.plot_data(awg,x_vector=t,y_vector=I,fitted_pars=fitted_pars,**options_rabi_statistics,iteration=iteration_rabi_statistics,plot_mode=0)
    rabi_freq_arr[i] = 1/(np.round(fitted_pars[1]))*1e9
    error_arr[i] = max(error)
    time_arr[i] = time.time()  - start




bad_index_arr = [idx for idx, element in enumerate(error_arr) if condition(element)]

# discard values with high errors

fig, ax= plt.subplots()
fig.suptitle('Rabi Statistics %03d'%(iteration_rabi_statistics))
ax.plot(time_arr,rabi_freq_arr, '-o', markersize = 3, c='C0')
ax.set_xlabel('time (sec)')
ax.set_ylabel('$\Omega_R$ (MHz)')

# save data
with open("E:\\generalized-markovian-noise\\%s\\rabi\\Rabi_statistics\\%s_data_%03d.csv"%(meas_device,'rabi_statistics',iteration_rabi_statistics),"w",newline="") as datafile:
    writer = csv.writer(datafile)
    writer.writerow(options_rabi_statistics.keys())
    writer.writerow(options_rabi_statistics.values())
    writer.writerow(rabi_freq_arr_clean)
    writer.writerow(now)
    writer.writerow(error_arr)

iteration_rabi_statistics += 1


#%% ramsey statistics
'''------------------------------------------------------Ramsey Statistics---------------------------------------------------'''
'''DESCRIPTION: Repeat Ramsey Measurement for a couple of hours to determine timescale of environmental fluctuations'''


qb.exp_pars = {
    'exp':                  'ramsey',
    'n_avg':                512,
    'x0':                   100e-9,
    'xmax':                 300e-6,
    'dx':                   2500e-9,
    'fsAWG':                0.6e9,
    'amp_q':                0.1,
    'active_reset':         False,
    'qubit_drive_freq':     qb.qb_pars['qb_freq']+detuning,
    'tomographic-axis':     'Z',
}

nReps = 1000
rep = np.arange(nReps)
detun_arr = np.zeros(nReps)
T_phi_arr = np.zeros(nReps)
data_arr = np.zeros((nReps,119))
# error_arr = np.zeros(nReps)
# time_arr = np.zeros(nReps)

# start = time.time()

for i in range(nReps):
    t,data,nSteps = qb.pulsed_exp(qb=qb_name,device_name=device_name,verbose=1,check_mixers=False,save_data=True)
    data_arr[i,:] = data
    fitted_pars,error = pt.fit_data(x_vector=t,y_vector=data,exp='ramsey',dx=t[-1]/nSteps*1e6,verbose=0)
    detun_arr[i] = fitted_pars[1]
    T_phi_arr[i] = fitted_pars[3]
    # time_arr[i] = time.time()  - start

# # discard values with high errors
# detun_arr_clean = np.zeros(nReps)
# T_phi_arr_clean = np.zeros(nReps)
# rep_clean = np.zeros(nReps)

# bad_index_arr = [idx for idx, element in enumerate(error_arr) if condition(element)]

# detun_arr_clean = np.delete(detun_arr,bad_index_arr)
# T_phi_arr_clean = np.delete(T_phi_arr,bad_index_arr)
# rep_clean = np.delete(rep,bad_index_arr)

# load data from previous measurement
# data = pd.read_csv("E:\\generalized-markovian-noise\\%s\\ramsey\\ramsey_statistics\\%s_data_%03d.csv"%(meas_device,'ramsey_statistics',iteration_ramsey_statistics),nrows=2,on_bad_lines='skip',skiprows=2,header=None).to_numpy(np.float64)
# re = data[0,:]
# detun_arr = data[1,:]

fig, (ax1,ax2) = plt.subplots(2,sharex=True)
ax1.plot(np.arange(0,nReps),detun_arr, '-o', markersize = 3, c='C0')
# ax1.plot(rep_clean,np.abs(detun_arr_clean), '-o', markersize = 3, c='C0')
ax1.set_ylim((0,1))
ax1.set_ylabel('$\Delta$ (MHz)')
ax2.plot(np.arange(0,nReps),T_phi_arr, '-o', markersize = 3, c='C0')

# ax2.plot(rep_clean,T_phi_arr_clean, '-o', markersize = 3, c='C0')
ax2.set_xlabel('time (sec)')
ax2.set_ylabel('$T_{\phi} (\mu s)$')
# ax2.set_ylim((0,10))

plt.gcf().text(1, 0.25, textstr, fontsize=14)


#%%
'''------------------------------------------------------Test of Markovianity---------------------------------------------------'''
'''DESCRIPTION: Repeat Ramsey Measurement for a couple of hours to determine timescale of environmental fluctuations'''

try:
    list_of_files = glob.glob('E:\\generalized-markovian-noise\\%s\\Ramsey\\markovianity_check\\*.csv'%(meas_device))
    latest_file = max(list_of_files, key=os.path.getctime)
    iteration_markovianity_check = int(latest_file[-7:-4].lstrip('0')) + 1
except:
    iteration_markovianity_check = 1

options_markovianity_check = {
    'nAverages':        options_ramsey['nAverages'],
    'Tmax':             options_ramsey['Tmax'],
    'stepSize':         options_ramsey['stepSize'],
    'integration_length': options_ramsey['integration_length'],
    'cav_resp_time':    options_ramsey['cav_resp_time'],
    'amp_q':     A_d,
    'measPeriod':       options_ramsey['measPeriod'],
    'qubitDriveFreq':   options_ramsey['qubitDriveFreq'],
    'active_reset':     True,
    'threshold':        threshold,
    'pi2Width':         options_ramsey['pi2Width'],
    'AC_freq':          options_ramsey['AC_freq'],
    'mu':               0.315,
    'sigma':            0.015,
    'B0':               0,
    'nu':               0,
    'tauk':             0,
    'source':           2,
    'noise_rate':       1,
}

nReps = 500

for j in range(1,13):
    data_arr = np.zeros((nReps,nSteps))
    #calibrate mixers
    expf.mixer_calib(sa,awg,'fine',mixer='qubit',threshold=-75,fc=options_markovianity_check['qubitDriveFreq'],amp=options_markovianity_check['amp_q'])
    expf.mixer_calib(sa,awg,'fine',mixer='ac',threshold=-75,fc=options_markovianity_check['AC_freq'],amp=options_markovianity_check['mu'])
    attn.setAttenuation(devID = 1, atten = 0)
    expf.mixer_calib(sa,daq,'fine',mixer='readout',threshold=-75,fc=7.2581e9,amp=0.7)
    attn.setAttenuation(devID = 1, atten = attenuation)

    for i in range(nReps):
        if i == 0:
            setup = [0,1,0]
        else:
            setup = [2,1,1]
            options_markovianity_check['nSteps'] = nSteps

        options_markovianity_check['noise_rate'] = j
        t,data,nSteps = expf.pulse(daq,awg,setup=setup,sequence='ramsey',**options_markovianity_check)
        data_arr[i,:] = data

    # save data
    with open("E:\\generalized-markovian-noise\\%s\\ramsey\\markovianity_check\\%s_data_awg_rate_%03d_kHz.csv"%(meas_device,'markovianity_check',int(1.2e6*2**(-(j-1)))),"w",newline="") as datafile:
        writer = csv.writer(datafile)
        writer.writerow(options_markovianity_check.keys())
        writer.writerow(options_markovianity_check.values())
        writer.writerow(t)
        writer.writerows(data_arr)

# load data from previous measurement
data = pd.read_csv("E:\\generalized-markovian-noise\\%s\\ramsey\\markovianity_check\\%s_data_%03d.csv"%(meas_device,'markovianity_check',iteration_markovianity_check),nrows=2,on_bad_lines='skip',skiprows=2,header=None).to_numpy(np.float64)
re = data[0,:]

fig, (ax1,ax2) = plt.subplots(2,sharex=True)
ax1.plot(time_arr,np.abs(detun_arr), '-o', markersize = 3, c='C0')
# ax1.plot(rep_clean,np.abs(detun_arr_clean), '-o', markersize = 3, c='C0')
ax1.set_ylim((0,1))
ax1.set_ylabel('$\Delta$ (MHz)')
fig.suptitle('Ramsey Statistics %03d'%(iteration_ramsey_statistics))
ax2.plot(time_arr,T_phi_arr, '-o', markersize = 3, c='C0')

# ax2.plot(rep_clean,T_phi_arr_clean, '-o', markersize = 3, c='C0')
ax2.set_xlabel('time (sec)')
ax2.set_ylabel('$T_{\phi} (\mu s)$')
# ax2.set_ylim((0,10))
textstr = "$\omega_d = 2\pi\\times%.4f$ GHz" %(options_markovianity_check['qubitDriveFreq']*1e-9)
plt.gcf().text(1, 0.25, textstr, fontsize=14)



iteration_ramsey_statistics += 1

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
    'amp_q':     A_d,
    'nSteps':           nSteps,
    'measPeriod':       200e-6,
    'qubitDriveFreq':   options_rabi['qubitDriveFreq']+detun,
    'sweep':            0,
    'pi2Width':         options_echo['pi2Width'],
    'pi2Width_Y':       28.5e-9,
    'pipulse_position': 117e-9,
    'mu':               0.305,
    'sigma':            0,
    'B0':               0,
    'nu':               0,
    'tauk':             0,
}

nReps = 250
rep = np.arange(nReps)
T2_arr = np.zeros(nReps)
error_arr = np.zeros(nReps)
now = np.zeros(nReps)
plot = 1

start = time.time()
for i in range(nReps):
    t,data,nSteps = expf.pulse(daq,awg,qubitLO,setup=[1,1,1],sequence='echo',**options_echo_statistics)
    T2_arr[i],error = pf.pulse_plot1d(sequence='echo',x_vector=t, y_vector=data,plot=plot,dt=options_echo_statistics['Tmax']*1e6/nSteps,qubitDriveFreq=options_echo_statistics['qubitDriveFreq'],amp_q=options_echo_statistics['amp_q'],fitting=1,AC_pars=options_echo_statistics['AC_pars'],RT_pars=options_echo_statistics['RT_pars'],pi2Width=options_echo_statistics['pi2Width'],iteration=iteration_echo_statistics)
    error_arr[i] = max(error)
    now[i] = time.time()  - start

# # discard values with high errors
T2_arr_clean = np.zeros(nReps)
rep_clean = np.zeros(nReps)
now_clean = np.zeros(nReps)


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
exp_pars = [options_echo_statistics['amp_q'],options_echo_statistics['qubitDriveFreq'],options_echo_statistics['AC_pars']]
with open("E:\\generalized-markovian-noise\\echo\\echo_statistics\\%s_data_%03d.csv"%('echo_statistics',iteration_echo_statistics),"w",newline="") as datafile:
    writer = csv.writer(datafile)
    writer.writerow(exp_pars)
    writer.writerow(T2_arr_clean)
    writer.writerow(now_clean)
    writer.writerow(error_arr)

iteration_echo_statistics += 1

#%% T1 (noise)
'''------------------------------------------------------T1 Statistics---------------------------------------------------'''
'''DESCRIPTION: Repeat Echo Measurement for a couple of hours to determine timescale of environmental fluctuations'''
detuning = 0
qb.wfm_pars = {
    't0':                   0.1e-6,
    'tmax':                 16e-6,
    'dt':                   0.2e-6,
    'fsAWG':                2.4e9,
    'mu':                   0,
    'sigma':                150e-3,
    }

qb.exp_pars = {
    'exp':                  'T1',
    'n_avg':                512,
    'x0':                   qb.wfm_pars['t0'],
    'xmax':                 qb.wfm_pars['tmax'],
    'dx':                   qb.wfm_pars['dt'],
    'fsAWG':                qb.wfm_pars['fsAWG'],
    'active_reset':         True,
    'qubit_drive_freq':     qb.qb_pars['qb_freq']+detuning,
    'tomographic-axis':     'Z',
}

timeSteps = round((qb.wfm_pars['tmax']-qb.wfm_pars['t0'])/qb.wfm_pars['dt'])-1
nReps =2000
data = np.zeros((nReps,timeSteps))
for i in range(nReps):
    print('repetition number',i)
    t,data[i,:],nSteps = qb.pulsed_exp(qb=qb_name,device_name=device_name, verbose=1,check_mixers=False,save_data=True)
    
    # Update Plots every 10 % of the measurement
    if (np.mod(i+1,400) == 0 ):
        fitted_pars,error = pt.fit_data(x_vector=t[1:],y_vector=np.mean(data[:,1:],0),exp='T1',dx=t[-1]/nSteps,verbose=0)
        pt.plot_T1_data(t[1:],np.mean(data[:,1:],0),fitted_pars,qb=qb_name,exp_pars=qb.exp_pars,qb_pars=qb.qb_pars,device_name=device_name,project=project)
        
pt.plot_T1_data(t[1:],np.mean(data[:,1:],0),fitted_pars,qb=qb_name,exp_pars=qb.exp_pars,qb_pars=qb.qb_pars,device_name=device_name,project=project)

# start = time.time()
# for i in range(nReps):
#     t,data,nSteps = expf.pulse(daq,awg,setup=[1,1,1],sequence='T1',**options_T1_statistics)
#     fitted_pars,error = pf.fit_data(x_vector=t,y_vector=data,sequence='T1',dt=t[-1]/nSteps,**options_T1_statistics)
#     T1_arr[i] = fitted_pars[1]
#     pf.plot_data(awg,x_vector=t,y_vector=I,fitted_pars=fitted_pars,sequence='T1',**options_T1_statistics,iteration=iteration_T1_statistics,plot_mode=0)
#     error_arr[i] = max(error)
#     now[i] = time.time()  - start

# # # discard values with high errors
# T1_arr_clean = np.zeros(nReps)
# rep_clean = np.zeros(nReps)
# now_clean = np.zeros(nReps)


# bad_index_arr = [idx for idx, element in enumerate(error_arr) if condition(element)]

# T1_arr_clean = np.delete(T1_arr,bad_index_arr)
# rep_clean = np.delete(rep,bad_index_arr)
# # time_arr = 14.5*np.arange(len(rep_clean))
# time_arr = 14.5*np.arange(len(now))

# fig, ax= plt.subplots()
# fig.suptitle('T1 Statistics %03d'%(iteration_T1_statistics))
# ax.plot(time_arr,T1_arr, '-o', markersize = 3, c='C0')
# ax.set_xlabel('time (sec)')
# ax.set_ylabel('$T_1 (\mu s)$')

#%% T2 (noise)
'''------------------------------------------------------T2 Statistics---------------------------------------------------'''
'''DESCRIPTION: Repeat T2 Measurement with different noise realizations '''

qb.wfm_pars = {
    't0':                   0.1e-6,
    'tmax':                 30e-6,
    'dt':                   0.3e-6,
    'fsAWG':                1.2e9,
    'mu':                   0,
    'sigma':                171e-3,
    }

qb.exp_pars = {
    'exp':                  'ramsey',
    'n_avg':                128,
   'x0':                   qb.wfm_pars['t0'],
   'xmax':                 qb.wfm_pars['tmax'],
   'dx':                   qb.wfm_pars['dt'],
   'fsAWG':                qb.wfm_pars['fsAWG'],
    'amp_q':                0.4,
    'active_reset':         False,
    'qubit_drive_freq':     qb.qb_pars['qb_freq'],
    'tomographic-axis':     'Z',
}

nReps = 200 
data = np.zeros((nReps,101))
for i in range(nReps):
    t,data[i,:],nSteps = qb.pulsed_exp(qb=qb_name,device_name=device_name, verbose=1,check_mixers=False,save_data=False)

fitted_pars,error = pt.fit_data(x_vector=t[1:],y_vector=np.mean(data,0),exp='ramsey',dx=t[-1]/nSteps,verbose=0)
pt.plot_ramsey_data(t,np.mean(data,0),fitted_pars,qb=qb_name,exp_pars=qb.exp_pars,qb_pars=qb.qb_pars,device_name=device_name,project=project)


#%% something else
'''------------------------------------------------------Measuring the Qubit's Bandwidth via Echo'---------------------------------------------------'''
'''DESCRIPTION: The goal of this experiment is to determine the lower bound of the qubit's bandwidth by doing echo with increasing noise amplitude'''

options_qubitBW_meas = {
  'sampling_rate':      options_echo['sampling_rate'],
    'nAverages':        128,
    'Tmax':             options_echo['Tmax'],
    'stepSize':         options_echo['stepSize'],
    'integration_length': options_echo['integration_length'],
    'cav_resp_time':    options_echo['cav_resp_time'],
    'amp_q':     A_d,
    'measPeriod':       600e-6,
    'qubitDriveFreq':   options_rabi['qubitDriveFreq']+detun,
    'sweep':            0,
    'active_reset':     False,
    'prePulseLength':   options_echo['prePulseLength'],
    'postPulseLength':   options_echo['postPulseLength'],
    'pi2Width':         options_echo['pi2Width'],
    'mu':               options_echo['mu'],
    'AC_freq':          options_echo['AC_freq'],
    'sigma':            0,
    'B0':               0,
    'nu':               0,
    'tauk':             0,

}

nReps = 600
sigma_arr = np.linspace(0.02,0.06,5)
data_arr = np.zeros((nReps,60))

for j in range(4,5):
    data_arr = np.zeros((nReps,nSteps))
    #calibrate mixers

    for i in trange(nReps,desc='Iteration'):
        if i % 200 == 0:
            expf.mixer_calib(sa,awg,'fine',mixer='qubit',threshold=-75,fc=options_qubitBW_meas['qubitDriveFreq'],amp=options_qubitBW_meas['amp_q'])
            expf.mixer_calib(sa,awg,'fine',mixer='ac',threshold=-75,fc=options_qubitBW_meas['AC_freq'],amp=options_qubitBW_meas['mu'])
            attn.setAttenuation(devID = 1, atten = 0)
            expf.mixer_calib(sa,daq,'fine',mixer='readout',threshold=-75,fc=7.2581e9,amp=0.7)
            attn.setAttenuation(devID = 1, atten = attenuation)
        else:
            pass
        setup = [0,1,1]
        options_qubitBW_meas['sigma'] = 0.06
        try:
            t,data,nSteps = expf.pulse(daq,awg,setup=setup,sequence='echo_v2',verbose=0,**options_qubitBW_meas)
            data_arr[i,:] = data
        except Exception as exc:
            print(exc)
            i -= 1

    # save data
    with open("E:\\generalized-markovian-noise\\%s\\echo\\qubitBW_meas\\f_AC_%d_MHz_sigma_%d_mV.csv"%(meas_device,int(options_qubitBW_meas['AC_freq']*1e-6),int(options_qubitBW_meas['sigma']*1e3)),"w",newline="") as datafile:
        writer = csv.writer(datafile)
        writer.writerow(options_qubitBW_meas.keys())
        writer.writerow(options_qubitBW_meas.values())
        writer.writerow(t)
        writer.writerows(data_arr)


#%% Ramsey Sweep
'''---------------------------------------Ramsey Parameter Sweep for Random Telegraph Noise-------------------------------------------'''

'''DESCRIPTION:
    1. Generate ~200 instances of telegraph noise for all parameter points and store them in separate folders
    2. Repeat ramsey measurement with different noise realizations
    3. Every x instances, do ramsey without RTN to get statistics for the bare T2*
    4. Repeat for every parameter point
'''

try:
    directory = 'E:\\generalized-markovian-noise\\%s\\sweep_data\\ramsey\\' %(meas_device)
    latest_sweep = max(glob.glob(os.path.join(directory, '*')), key=os.path.getmtime)
    sweep_count = int(latest_sweep[-3:].lstrip('0')) + 1
except:
    sweep_count = 1

numRealizations = 100
b_measurements = 100
numIterations = numRealizations + b_measurements
interval = int(numRealizations/b_measurements) # how often to take background data without generalized markovian noise
# B0 = np.linspace(0.0005,0.05,10)
# B0 = [0.0005 ,0.05]
B0_kHz = [25,500,750,1000,1250,2000,2150,2250,2500,2750,3000]
# B0_kHz = np.arange(100,5000,step=500)
# B0_kHz = np.arange(25,4000,step=100)
# B0_kHz = [25,500,1000,2000,3000]
B0 = np.around([expf.convertB0_kHz_to_mV(i) for i in B0_kHz],5)*1e-3
B0[i] = expf.convertB0_kHz_to_mV(np.sqrt(2/9*1e6*(1-1/tau[k])**2+2*(2*np.pi*nu[j])**2)/(2*np.pi))*1e-3
# nu = np.concatenate(([i/2 for i in B0_kHz],[3000]))
# nu = np.concatenate((np.linspace(0.1,2,5),np.linspace(3,20,4),np.linspace(25,100,5),[1000]))
# nu = np.concatenate((np.linspace(0.1,2,4),np.linspace(3,100,6),[1000]))
# B0 =[ expf.convertB0_kHz_to_mV(1500)*1e-3]
# nu = [10,500,1e3,1.25e3,1.5e3,1.75e3,2e3,2.25e3,2.5e3]

# B0 = [expf.convertB0_kHz_to_mV(np.sqrt(2/9*1e12*(1-1/tau[0])**2+2*(2*np.pi*i*1e3)**2)*1e-3) for i in nu]
# nu = [0.01,10]
# nu = [1500]

B0_kHz = [750,1500,2122,3000]
# B0_kHz = [2*i for i in nu]
# B0_2_kHz = [np.sqrt(2/9*1e6+2*(2*np.pi*i)**2)/(2*np.pi) for i in nu]
# B0_kHz = np.concatenate((B0_1_kHz,B0_2_kHz))
# B0_kHz = np.sort(B0_kHz)
B0 = np.around([expf.convertB0_kHz_to_mV(i) for i in B0_kHz],5)*1e-3

# nu = np.concatenate(([i/2 for i in B0_kHz],[8000]))
# tau = np.concatenate((np.linspace(0.1,2,3),np.linspace(3,8,3),[20]))
# tau = [0.1,0.5,1,2,5,10,20,100]
# tau = [0.1,10]
tau = [0.5]
# tau = [0.1,2,10,100]

# initialize arrays for data
stepSize_background = 200e-9
Tmax_background = 6e-6
stepSize = 13e-9
Tmax = 6e-6
nPoints, nStepsBackground, pulse_length_increment, pulse_length_start = expf.calc_nSteps(sequence='ramsey',fsAWG=1.2e9,stepSize=stepSize_background,Tmax=Tmax_background)
nPoints,nSteps, pulse_length_increment, pulse_length_start = expf.calc_nSteps(sequence='ramsey',fsAWG=1.2e9,stepSize=stepSize,Tmax=Tmax)
# initialize data arrays
backData, measData = expf.init_arrays(numRealizations=numRealizations,interval=interval,nPointsBackground=nStepsBackground-1,nPoints=nSteps-1)

# make appropriate directories
sweep_name = 'sweep_%03d'%(sweep_count)
parent_dir = 'E:\\generalized-markovian-noise\\%s\\sweep_data\\ramsey\\'%(meas_device)
main_path = os.path.join(parent_dir,sweep_name)
data_path = os.path.join(main_path,'data')
noise_path = os.path.join(main_path,'noise_instances')
plot_path = os.path.join(main_path,'plot_images')
try:
    os.mkdir(main_path)
    os.mkdir(data_path)
    os.mkdir(noise_path)
    os.mkdir(plot_path)
except:
    print('Directories already exist')

# generate noise instances
str_gen_noise = time.time()
use_wk = 0
expf.gen_noise_realizations(par1_arr=tau,par2_arr=nu,numRealizations=numRealizations,nPoints=nPoints,T_max=Tmax,sweep_count=sweep_count,
                            sequence='ramsey',meas_device=meas_device,wk=use_wk)
end_gen_noise = time.time()
print('Generating noise realizations took: %.1f s' %(end_gen_noise-str_gen_noise))

plt.close('all')

# sweep ramsey
options_ramsey_par_sweep = {
    'sampling_rate':    1.2e9,
    'active_reset':     True,
    'threshold':        options_ramsey['threshold'],
    'prePulseLength':   options_ramsey['prePulseLength'],
    'postPulseLength':  options_ramsey['postPulseLength'],
    'nAverages':        128,
    'Tmax':             3e-6,
    'stepSize':         100e-9,
    'pi2Width':         1/2*pi_pulse*1e-9,
    'amp_q':     A_d,
    'cav_resp_time':    options_ramsey['cav_resp_time'],
    'measPeriod':       options_ramsey['measPeriod'],
    'qubitDriveFreq':   options_ramsey['qubitDriveFreq'],
    'mu':               0,
    'sigma':            0.1,
    'B0':               0,
    'nu':               0,
    'tauk':             0,
    'AC_freq':          0,
    'source':           2,
    'sweep':            1,
    'wk':               use_wk,
    }


start_sweep = time.time()
errors = 0

# calibrate UHFQA before running sweep
daq.setInt('/dev2528/system/calib/calibrate',1)

# a = 1

for i in trange(len(B0),desc='B0 Loop'):
    for j in trange(len(nu),desc='nu Loop'):
        for k in trange(len(tau),desc='tauk Loop'):
            start_point_t = time.time()
            # load noise instances
            # B0[i] = np.around(expf.convertB0_kHz_to_mV(B0_kHz[i]),5)*1e-3
            # B0[i] = expf.convertB0_kHz_to_mV(np.sqrt(2/9*1e6*(1-1/tau[k])**2+2*(2*np.pi*nu[j])**2)/(2*np.pi))*1e-3 # amplitude is in linear units. The factor of 1e-3 at the end is to convert mV to Volts.
            noise_realizations = B0[i]*expf.pull_wfm(sweep_name,nu[j],tau[k],sequence='ramsey')
            # filename = 'B0_%d_uV_nu_%d_Hz_tau_%d_ns_%d' %(round(B0[i]*1e6),round(nu[j]*1e3),round(tau[k]*1e3),int(q))
            filename = 'B0_%d_uV_nu_%d_Hz_tau_%d_ns' %(round(B0[i]*1e6),round(nu[j]*1e3),round(tau[k]*1e3))
            b = 0 # index that keeps track of background measurements
            n = 0 # index that keeps track of noisy measurements
            m = 0 # index that keeps track of overall measurements
            # optimize mixers
            expf.mixer_calib(sa,awg,'fine',mixer='qubit',threshold=-75,fc=options_ramsey_par_sweep['qubitDriveFreq'],amp=options_ramsey_par_sweep['amp_q'],config='XY')
            if options_ramsey_par_sweep['mu'] != 0:
                expf.mixer_calib(sa,awg,'fine',mixer='ac',threshold=-75,fc=options_ramsey_par_sweep['AC_freq'],amp=options_ramsey_par_sweep['mu'])
            else:
                pass
            attn.setAttenuation(devID = 1, atten = 0)
            expf.mixer_calib(sa,daq,'fine',mixer='readout',threshold=-75,fc=7.2578e9,amp=0.7)
            attn.setAttenuation(devID = 1, atten = attenuation)

            #calibrate pi_pulse
            t,data,nSteps = expf.pulse(daq,awg,setup=[0,1,0],sequence='rabi',**options_rabi)
            fitted_pars,error = pf.fit_data(x_vector=t,y_vector=data,sequence='rabi',dt=t[-1]/nSteps,**options_rabi)
            options_ramsey_par_sweep['pi2Width'] = np.round(1/4*fitted_pars[1])*1e-9

            for m in trange(numIterations,desc='Measurements'):
                if m % (interval + 1) == 0 and b <= b_measurements: # executes a background measurement every n=interval noisy measurements
                # get background T2* every 10 or so measurements
                    options_ramsey_par_sweep['nAverages'] = 128
                    options_ramsey_par_sweep['B0'] = 0
                    options_ramsey_par_sweep['Tmax'] = Tmax_background
                    options_ramsey_par_sweep['stepSize'] = stepSize_background
                    white_noise_instance = np.random.normal(loc=options_ramsey_par_sweep['mu'], scale=options_ramsey_par_sweep['sigma'], size=nPoints)
                    print('----------------------------------\nExecuting background Ramsey measurement')
                    try:
                        t2,data,nSteps = expf.pulse(daq,awg,setup=[0,1,0],sequence='ramsey',verbose=0,white_noise_instance=white_noise_instance,
                                                     **options_ramsey_par_sweep)
                        backData[b,:] = data
                        b += 1
                        m += 1
                    except Exception as exc:
                        print(exc)
                        errors += 1
                else:
                    if m % (interval + 1) == 1:
                        setup = [0,1,0] # upload a new seqc file, and re-configure QA unit
                    else:
                        setup = [2,1,0] # replace noisy waveforms only; everything else stays the same
                    options_ramsey_par_sweep['B0'] = B0[i]
                    options_ramsey_par_sweep['nu'] = nu[j]
                    options_ramsey_par_sweep['tauk'] = tau[k]
                    options_ramsey_par_sweep['Tmax'] = Tmax
                    options_ramsey_par_sweep['nAverages'] = 128
                    options_ramsey_par_sweep['stepSize'] = stepSize
                    print('----------------------------------\nImplementing noise realization %d' %(n+1))
                    noise_instance = noise_realizations[n,:]
                    try:
                        t1,data,nSteps = expf.pulse(daq,awg,setup=setup,sequence='ramsey',noise_instance=noise_instance,white_noise_instance=white_noise_instance,
                                                     verbose=0,**options_ramsey_par_sweep)
                        measData[n,:] = data
                        n += 1
                        m += 1
                    except Exception as exc:
                        print(exc)
                        errors += 1

            # save data after each parameter sweep point
            with open("E:\\generalized-markovian-noise\\%s\\sweep_data\\ramsey\\%s\\data\\data_%s.csv"%(meas_device,sweep_name,filename),"w",newline="") as datafile:
                writer = csv.writer(datafile)
                writer.writerow(options_ramsey_par_sweep.keys())
                writer.writerow(options_ramsey_par_sweep.values())
                writer.writerow(['Background Time Data'])
                writer.writerow(t2)
                writer.writerow(['Time Data'])
                writer.writerow(t1)
                writer.writerow(['Background Data'])
                writer.writerows(backData)
                writer.writerow(['Data'])
                writer.writerows(measData)

            # a += 1
            end_point_t = time.time()
            print('Parameter point took %s\n%d errors'%(end_point_t-start_point_t,errors))






end_sweep = time.time()
print('Total Sweep Duration: %.1f s = %.1f hours = %.1f days' %(end_sweep-start_sweep,(end_sweep-start_sweep)/3600,(end_sweep-start_sweep)/(3600*24)))
sweep_count += 1


'''---------------------------------------Echo Parameter Sweep for Random Telegraph Noise-------------------------------------------'''

'''DESCRIPTION:
    1. Generate ~200 instances of telegraph noise for all parameter points and store them in separate folders
    2. Repeat ramsey measurement with different noise realizations
    3. Every x instances, do ramsey without RTN to get statistics for the bare T2*
    4. Repeat for every parameter point
'''

try:
    directory = 'E:\\generalized-markovian-noise\\%s\\sweep_data\\echo\\' %(meas_device)
    latest_sweep = max(glob.glob(os.path.join(directory, '*')), key=os.path.getmtime)
    sweep_count = int(latest_sweep[-3:].lstrip('0')) + 1
except:
    sweep_count = 1

numRealizations = 128
b_measurements = 32
numIterations = numRealizations + b_measurements
interval = int(numRealizations/b_measurements) # how often to take background data without generalized markovian noise
B0_kHz = [25,100,500,1000,1500,2000,2500,3000,3500,4000]
B0 = np.around([expf.convertB0_kHz_to_mV(i) for i in B0_kHz],4) #the relationship between B0 and frequency of oscillations is Omega_R = 25 MHz * A_q
# B0 = [0.0005 ,0.05]
# B0 = [0.039]
# nu = np.concatenate((np.linspace(0.1,2,5),np.linspace(3,20,4),np.linspace(25,100,5),[1000]))
# nu = np.concatenate((np.linspace(0.1,2,4),np.linspace(3,100,6),[1000]))
nu = [i/2 for i in B0_kHz]
# nu = [0.01,10]
# nu = [1e3]
# tau = np.concatenate((np.linspace(0.1,2,3),np.linspace(3,20,2),[1000]))
# tau = [0.1,0.5,1.5,3,5,10,20,100]
# tau = [0.1,10]
tau = [1e12]

# initialize arrays for data
Tmax_background = 40e-6
Tmax = 150e-6
stepSize_background = 1000e-9
stepSize = 1500e-9

nPoints, nStepsBackground, pulse_length_increment, pulse_length_start = expf.calc_nSteps(sequence='echo',fsAWG=1.2e9,stepSize=stepSize_background,Tmax=Tmax_background)
nPoints,nSteps, pulse_length_increment, pulse_length_start = expf.calc_nSteps(sequence='echo',fsAWG=1.2e9,stepSize=stepSize,Tmax=Tmax)
backData, measData = expf.init_arrays(numRealizations=numRealizations,interval=interval,nPointsBackground=nStepsBackground,nPoints=nSteps)

# make appropriate directories
sweep_name = 'sweep_%03d'%(sweep_count)
parent_dir = 'E:\\generalized-markovian-noise\\%s\\sweep_data\\echo\\'%(meas_device)
main_path = os.path.join(parent_dir,sweep_name)
data_path = os.path.join(main_path,'data')
noise_path = os.path.join(main_path,'noise_instances')
plot_path = os.path.join(main_path,'plot_images')
try:
    os.mkdir(main_path)
    os.mkdir(data_path)
    os.mkdir(noise_path)
    os.mkdir(plot_path)
except:
    print('Directories already exist')

# generate noise instances
str_gen_noise = time.time()
expf.gen_noise_realizations(par1_arr=tau,par2_arr=nu,numRealizations=numRealizations,nPoints=nPoints,T_max=Tmax,sweep_count=sweep_count,sequence='echo',meas_device=meas_device)
end_gen_noise = time.time()
print('Generating noise realizations took: %.1f s' %(end_gen_noise-str_gen_noise))

plt.close('all')

# sweep ramsey
options_echo_par_sweep = {
    'sampling_rate':    1.2e9,
    'active_reset':     False,
    'prePulseLength':   options_echo['prePulseLength'],
    'postPulseLength':  options_echo['postPulseLength'],
    'nAverages':        128,
    'Tmax':             3e-6,
    'stepSize':         100e-9,
    'pi2Width':         1/2*pi_pulse*1e-9,
    'amp_q':     A_d,
    'cav_resp_time':    options_echo['cav_resp_time'],
    'measPeriod':       options_echo['measPeriod'],
    'qubitDriveFreq':   options_echo['qubitDriveFreq'],
    'mu':               0.3,
    'sigma':            0.2,
    'B0':               0,
    'nu':               0,
    'tauk':             0,
    'AC_freq':          acStarkLO.getFrequency(),
    'source':           2,
    'sweep':            1,
    'phi':              1,
    }

start_sweep = time.time()

errors = 0

for i in trange(2,len(B0),desc='B0 Loop'):
    for j in trange(len(nu),desc='nu Loop'):
        for k in trange(len(tau),desc='tauk Loop'):

            # load noise instances
            noise_realizations = B0[i]*expf.pull_wfm(sweep_name,nu[j],tau[k],sequence='echo')
            start_point_t = time.time()
            filename = 'B0_%d_uV_nu_%d_Hz_tau_%d_ns' %(round(B0[i]*1e3),round(nu[j]*1e3),round(tau[k]*1e3))
            b = 0 # index that keeps track of background measurements
            n = 0 # index that keeps track of noisy measurements
            m = 0 # index that keeps track of overall measurements

            # optimize mixers
            expf.mixer_calib(sa,awg,'fine',mixer='qubit',threshold=-75,fc=options_echo_par_sweep['qubitDriveFreq'],amp=options_rabi['amp_q'],config='echo')
            attn.setAttenuation(devID = 1, atten = 0)
            expf.mixer_calib(sa,daq,'fine',mixer='readout',threshold=-75,fc=7.2581e9,amp=0.7)
            attn.setAttenuation(devID = 1, atten = attenuation)

            #calibrate pi_pulse
            t,data,nPoints = expf.pulse(daq,awg,setup=[0,1,0],sequence='rabi',**options_rabi)
            fitted_pars,error = pf.fit_data(x_vector=t,y_vector=data,sequence='rabi',dt=t[-1]/nPoints,**options_rabi)
            options_echo_par_sweep['pi2Width'] = np.round(1/4*fitted_pars[1])*1e-9

            for m in trange(numIterations,desc='Measurements'):
                if m % (interval + 1) == 0 and b <= b_measurements: # executes a background measurement every n=interval noisy measurements
                    # get background T2* every 10 or so measurements
                    options_echo_par_sweep['nAverages'] = 128
                    options_echo_par_sweep['B0'] = 0
                    options_echo_par_sweep['Tmax'] = Tmax_background
                    options_echo_par_sweep['stepSize'] = stepSize_background
                    print('----------------------------------\nExecuting background Ramsey measurement')
                    # try:
                    t2,data,nPoints = expf.pulse(daq,awg,setup=[0,1,0],sequence='echo',**options_echo_par_sweep)
                    backData[b,:] = data
                    b += 1
                    # m += 1
                    # except:
                        # print('Measurement Failed')
                        # errors += 1
                else:
                    if m % (interval + 1) == 1:
                        setup = [0,1,0] # upload a new seqc file, and re-configure QA unit
                    else:
                        setup = [2,1,1] # replace noisy waveforms only; everything else stays the same
                    options_echo_par_sweep['B0'] = B0[i]
                    options_echo_par_sweep['nu'] = nu[j]
                    options_echo_par_sweep['tauk'] = tau[k]
                    options_echo_par_sweep['Tmax'] = Tmax
                    options_echo_par_sweep['nAverages'] = 128
                    options_echo_par_sweep['stepSize'] = stepSize
                    print('----------------------------------\nImplementing noise realization %d' %(n+1))
                    noise_instance = noise_realizations[n,:]
                    # try:
                    t1,data,nPoints = expf.pulse(daq,awg,setup=setup,sequence='echo',noise_instance=noise_instance,**options_echo_par_sweep)
                    measData[n,:] = data
                    n += 1
                    # m += 1
                    # except:
                    #     print('Measurement Failed')
                    #     errors += 1

            # save data after each parameter sweep point
            with open("E:\\generalized-markovian-noise\\%s\\sweep_data\\echo\\%s\\data\\data_%s.csv"%(meas_device,sweep_name,filename),"w",newline="") as datafile:
                writer = csv.writer(datafile)
                writer.writerow(options_echo_par_sweep.keys())
                writer.writerow(options_echo_par_sweep.values())
                writer.writerow(['Background Time Data'])
                writer.writerow(t2)
                writer.writerow(['Time Data'])
                writer.writerow(t1)
                writer.writerow(['Background Data'])
                writer.writerows(backData)
                writer.writerow(['Data'])
                writer.writerows(measData)

            end_point_t = time.time()
            print('Parameter point took %s\n%d errors'%(end_point_t-start_point_t,errors))

end_sweep = time.time()
print('Total Sweep Duration: %.1f s = %.1f hours = %.1f days' %(end_sweep-start_sweep,(end_sweep-start_sweep)/3600,(end_sweep-start_sweep)/(3600*24)))
sweep_count += 1
