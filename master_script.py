"""
Created on Thu Apr 14 11:35:33 2022

@author: Evangelos Vlachos <evlachos@usc.edu>
@edits: Malida Hecht <mohecht@usc.edu>
"""

#%% Initialization
from qubit import qubit
import numpy as np
import plotting as pt
import matplotlib.pyplot as plt
from utilities import compute_bloch,compute_rho,compute_coherence,compute_purity

device_name = "WM1"
project = 'coherence-stabilization'
qb_name = 'qb6'
qb = qubit(qb_name)

#%% Spectroscopy (Resonator)
'''-----------------------------------------------------Resonator Spectroscopy------------------------------------------------------'''

freqs = np.arange(start=6.7035,stop=6.7055,step=25e-6) # frequencies are in GHz

qb.exp_pars = {
    'exp':                  'spectroscopy',
    'n_avg':                128,
    'element':              'rr',
    'rr_reset_time':        30e-6,
    'satur_dur':            2e-6,
    'rr_atten':             30,
    'on_off':               True,
    'tomographic_axis':     'Z',
    }

p_data,I,Q = qb.rr_spectroscopy(freqs,device_name)
fc = pt.rr_spec_plot(freq=freqs,I=I,Q=Q,mag=p_data,exp_pars=qb.exp_pars,qb_pars=qb.qb_pars,df=1e9*(freqs[1]-freqs[0]),find_pks=True)
qb.update_qb_value('rr_LO',fc*1e9)

#%% Spectroscopy (Qubit)

'''-----------------------------------------------------Qubit Spectroscopy------------------------------------------------------'''

freqs = np.arange(start=4.33,stop=4.5,step=100e-6) # frequencies are in GHz

qb.exp_pars = {
    'exp':                  'spectroscopy',
    'n_avg':                256,
    'element':              'qubit',
    'qubit_reset_time':     200e-6,
    'amp_q':                0.03,
    'satur_dur':            40e-6,
    'on_off':               False,
    'tomographic-axis':     'Z',
    'fsAWG':                2.4e9,
    }

p_data,I,Q = qb.qb_spectroscopy(freqs,device_name)
pt.qb_spec_plot(freq=freqs,I=I,Q=Q,mag=p_data*1e3,exp_pars=qb.exp_pars,qb_pars=qb.qb_pars,find_pks=True)

## Need to add title to this plot ^^


#%% Time Rabi
'''-----------------------------------------------------Time Rabi------------------------------------------------------'''
detuning = 0

qb.exp_pars = {
    'initial-state':        '0',
    'exp':                  't-rabi',
    'n_avg':                128,
    'x0':                   13e-9,
    'xmax':                 2e-6,
    'dx':                   6e-9,
    'fsAWG':                2.4e9,
    'amp_q':                0.3,
    'active_reset':         False,
    'qubit_drive_freq':     qb.qb_pars['qb_freq']+detuning,
    'tomographic-axis':     'Z',
    }

t,data,nSteps = qb.pulsed_exp(qb=qb_name,verbose=1,device_name=device_name,check_mixers=False,save_data=False)
# plot data
fitted_pars,error = pt.fit_data(x_vector=t,y_vector=data,exp='t-rabi',dx=t[-1]/nSteps*1e6,verbose=0)
pt.plot_t_rabi_data(t,data,fitted_pars,qb=qb_name,exp_pars=qb.exp_pars,qb_pars=qb.qb_pars,device_name=device_name,project=project)

#%% Power Rabi
'''-----------------------------------------------------Power Rabi------------------------------------------------------'''

qb.exp_pars = {
    'exp':                  'p-rabi',
    'n_avg':                512,
    'x0':                   0,
    'xmax':                 0.6,
    'dx':                   10e-3,
    'fsAWG':                2.4e9,
    'active_reset':         False,
    'qubit_drive_freq':     qb.qb_pars['qb_freq']+detuning,
    'tomographic-axis':     'Z',
    }

amp,data,nSteps = qb.pulsed_exp(qb=qb_name,device_name=device_name,verbose=1,check_mixers=False)
# plot data
fitted_pars,error = pt.fit_data(x_vector=amp,y_vector=data,exp='p-rabi',dx=qb.dx,verbose=0)
pt.plot_p_rabi_data(amp,data,fitted_pars,qb=qb_name,exp_pars=qb.exp_pars,qb_pars=qb.qb_pars,device_name=device_name,project=project)

qb.update_pi(pi_amp=fitted_pars[1]/2)


#%% Single Shot Experiment

qb.exp_pars = {
    'exp':                  'single-shot',
    'theta':                0,
    'num_samples':          512,
    'n_avg':                1,
    'fsAWG':                2.4e9,
    }

data = qb.single_shot(device_name,qb_name)

#make 2D histogram
pt.plot_single_shot(data, qb.exp_pars,qb.qb_pars,qb.iteration)

#%% T1 Measurement

qb.wfm_pars = {
    't0':                   0.1e-6,
    'tmax':                 100e-6,
    'dt':                   1e-6,
    'fsAWG':                0.6e9,
    'mu':                   0,
    'sigma':                0e-3,
    }

qb.exp_pars = {
    'exp':                  'T1',
    'n_avg':                512,
    'x0':                   qb.wfm_pars['t0'],
    'xmax':                 qb.wfm_pars['tmax'],
    'dx':                   qb.wfm_pars['dt'],
    'fsAWG':                qb.wfm_pars['fsAWG'],
    'active_reset':         False,
    'qubit_drive_freq':     qb.qb_pars['qb_freq']+detuning,
    'tomographic-axis':     'Z',
}


t,data,nSteps = qb.pulsed_exp(qb=qb_name,device_name=device_name, verbose=1,check_mixers=False,save_data=False)
# plot data
fitted_pars,error = pt.fit_data(x_vector=t,y_vector=data,exp='T1',dx=t[-1]/nSteps,verbose=0)
pt.plot_T1_data(t,data,fitted_pars,qb=qb_name,exp_pars=qb.exp_pars,qb_pars=qb.qb_pars,device_name=device_name,project=project)


#%% Ramsey Experiment

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

t,data,nSteps = qb.pulsed_exp(qb=qb_name,device_name=device_name,verbose=1,check_mixers=False,save_data=True)

# fit & plot data
fitted_pars,error = pt.fit_data(x_vector=t,y_vector=data,exp='ramsey',dx=t[-1]/nSteps*1e6,verbose=0)
pt.plot_ramsey_data(t,data,fitted_pars,qb=qb_name,exp_pars=qb.exp_pars,qb_pars=qb.qb_pars,device_name=device_name,project=project)
# qb.update_qb_value('qb_freq', qb.exp_pars['qubit_drive_freq']-fitted_pars[1]*1e6)

#%% tomography Experiment

qb.exp_pars = {
    'exp':                  'tomography',
    'n_avg':                2,
    'x0':                   100e-9,
    'xmax':                 120e-6,
    'dx':                   1e-6,
    'fsAWG':                1.2e9,
    'amp_q':                0.1,
    'active_reset':         False,
    'qubit_reset_time':     2,
    'qubit_drive_freq':     qb.qb_pars['qb_freq']+detuning,
    'tomographic-axis':     'Z',
}

t,data,nSteps = qb.pulsed_exp(qb=qb_name,device_name=device_name,verbose=1,check_mixers=False,save_data=True)

# fit & plot data
fitted_pars,error = qb.fit_data(x_vector=t,y_vector=data,dx=t[-1]/nSteps*1e6,verbose=0)
qb.plot_ramsey_data(t,data,fitted_pars,qb=qb_name,exp_pars=qb.exp_pars,qb_pars=qb.qb_pars,device_name=device_name,project=project)

#%% Coordinate System Calibration
'''Execute time-based rabi experiment, and measure along 3 different axis at the end'''

qb.exp_pars = {
    'exp':                  't-rabi',
    'n_avg':                1024,
    'x0':                   13e-9,
    'xmax':                 0.4e-6,
    'dx':                   6e-9,
    'fsAWG':                2.4e9,
    'amp_q':                0.4,
    'active_reset':         False,
    'qubit_drive_freq':     qb.qb_pars['qb_freq']+detuning,
}

t,data_cal,nSteps = qb.tomography_calibration(qb=qb_name,device_name=device_name,verbose=1,check_mixers=False,save_data=False)
# fit & plot data
calib_states = qb.cal_coord(data_cal)
pt.tom_calib_plot(x_data=t, y_data=data_cal, coords=calib_states)

#%% coherence-stabilization Experiment

qb.wfm_pars = {
    'x0':                   0.1e-6,
    'xmax':                 200e-6,
    'dx':                   4e-6,
    'fsAWG':                0.6e9,
    'mu':                   0,
    'sigma':                100e-3,
    'T1':                   2.5e-6,         
    }

qb.exp_pars = {
    'exp':                  'coherence-stabilization',
    'initial-state':        '7',
    'n_avg':                2048,
    'n_realizations':       1,
    'x0':                   qb.wfm_pars['x0'],
    'xmax':                 qb.wfm_pars['xmax'],
    'dx':                   qb.wfm_pars['dx'],
    'fsAWG':                qb.wfm_pars['fsAWG'],
    'amp_q':                0.1,
    'active_reset':         False,
    'qubit_drive_freq':     qb.qb_pars['qb_freq']+detuning,
}

wfms,data,nSteps = qb.coherence_stabilization(qb=qb_name,device_name=device_name,verbose=1,save_data=False)
t = np.linspace(qb.wfm_pars['x0'],qb.wfm_pars['xmax'],data.shape[1])
# fit & plot data
v0 = compute_bloch(data[:,0], calib_states)
vf = compute_bloch(data[:,-1], calib_states)
rho_0 = compute_rho(v0)
rho_f = compute_rho(vf)

pt.tom_calib_plot(x_data=t, y_data=data, coords=calib_states,data='data')
pt.plot_tomography(rho_0,qb.exp_pars['initial-state'], tmax=t[0]*1e6,cal_states=calib_states)
pt.plot_tomography(rho_f,qb.exp_pars['initial-state'], tmax=t[-1]*1e6,cal_states=calib_states)

pt.plot_coherence(t,data,wfms,qb.exp_pars,qb.qb_pars,qb.wfm_pars,calib_states)

#%% Virtual Z-Gate Experiment

qb.exp_pars = {
    'n_avg':                8192,
    'x0':                   0,
    'xmax':                 720,
    'dx':                   2,
    'fsAWG':                0.6e9,
    'active_reset':         False,
    'qubit_drive_freq':     qb.qb_pars['qb_freq']+detuning,
    'tomographic-axis':     'Z',
}


phase,data,nSteps = qb.pulsed_exp(qb=qb_name,device_name=device_name,verbose=1,check_mixers=False,save_data=True)
# fit & plot data
fitted_pars,error = qb.fit_data(x_vector=phase,y_vector=data,dx=phase[-1]/nSteps,verbose=0)
qb.plot_data(phase,data,fitted_pars,qb=qb_name)

#%% Echo Measurement

qb.exp_pars = {
    'n_avg':                2048,
    'x0':                   50e-9,
    'xmax':                 4e-6,
    'dx':                   0.05e-6,
    'fsAWG':                0.6e9,
    'amp_q':                0.1,
    'active_reset':         False,
    'qubit_reset_time':     300e-6,
    'qubit_drive_freq':     qb.qb_pars['qb_freq']+detuning,
    'tomographic-axis':     'Z',
}


t,data,nSteps = qb.pulsed_exp(qb=qb_name,device_name=device_name,verbose=1,check_mixers=False,save_data=True)
fitted_pars,error = qb.fit_data(x_vector=t,y_vector=data,dx=t[-1]/nSteps*1e6,verbose=0)
qb.plot_echo_data(x_vector=t,y_vector=data,fitted_pars=fitted_pars)


#%% Mixer Optimization
qb.min_leak(inst=qb.awg,f_LO=qb.qb_pars['qb_LO'],mixer='qubit',cal='lo',plot=True,span=0.1e6)
qb.min_leak(inst=qb.qa,f_LO=qb.qb_pars['rr_LO'],mixer='rr',cal='lo',plot=True,span=0.1e6)


# qb.min_leak(inst=qb.awg,f_LO=qb.qb_pars['qb_LO'],f_IF=qb.qb_pars['qb_IF']+detuning,cal='ssb',amp=0.4,threshold=-60,span=0.1e6)

qb.suppr_image(inst=qb.awg,f_LO=qb.qb_pars['qb_LO'],f_IF=qb.qb_pars['qb_IF'],amp=0.4)







