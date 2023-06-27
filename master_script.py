"""
Created on Thu Apr 14 11:35:33 2022

@author: Evangelos Vlachos <evlachos@usc.edu>
@edits: Malida Hecht <mohecht@usc.edu>
"""

#%% Initialization
from qubit import qubit
import numpy as np
import plotting as pt

device_name = "WM1"
project =  'dynamical-decoupling'
qb_name = 'qb6'
qb = qubit(qb_name)

#%% Spectroscopy (Resonator)
'''-----------------------------------------------------Resonator Spectroscopy------------------------------------------------------'''

freqs = np.arange(start=6.7035,stop=6.706,step=50e-6) # frequencies are in GHz

qb.exp_pars = {
    'n_avg':                512,
    'element':              'rr',
    'rr_reset_time':        20e-6,
    'satur_dur':            2e-6,
    'rr_atten':             30,
    'on_off':               True,
    'tomographic_axis':     'Z',
    }

p_data,I,Q = qb.rr_spectroscopy(freqs)
fc = pt.rr_spec_plot(freq=freqs,I=I,Q=Q,mag=p_data,exp_pars=qb.exp_pars,qb_pars=qb.qb_pars,df=1e9*(freqs[1]-freqs[0]),find_pks=True)
qb.update_qb_value('rr_LO',fc*1e9)

#%% Spectroscopy (Qubit)

'''-----------------------------------------------------Qubit Spectroscopy------------------------------------------------------'''

freqs = np.arange(start=3.8,stop=4,step=100e-6) # frequencies are in GHz

qb.exp_pars = {
    'n_avg':                512,
    'element':              'qubit',
    'qubit_reset_time':     200e-6,
    'amp_q':                0.1,
    'satur_dur':            40e-6,
    'on_off':               True,
    'tomographic-axis':     'Z',
    'fsAWG':                2.4e9,
    }

p_data,I,Q = qb.qb_spectroscopy(freqs)
pt.qb_spec_plot(freq=freqs,I=I,Q=Q,mag=p_data*1e3,exp_pars=qb.exp_pars,qb_pars=qb.qb_pars,find_pks=True)

## Need to add title to this plot ^^


#%% Time Rabi
'''-----------------------------------------------------Time Rabi------------------------------------------------------'''
detuning = -1.67e6

qb.exp_pars = {
    'n_avg':                256,
    'x0':                   13e-9,
    'xmax':                 1e-6,
    'dx':                   6e-9,
    'fsAWG':                2.4e9,
    'amp_q':                0.3,
    'active_reset':         False,
    'qubit_reset_time':     300e-6,
    'qubit_drive_freq':     qb.qb_pars['qb_freq']+detuning,
    'tomographic-axis':     'Z',
    }

t,data,nSteps = qb.pulsed_exp(exp='t-rabi',verbose=1,check_mixers=False)
# plot data
fitted_pars,error = pt.fit_data(x_vector=t,y_vector=data,exp='t-rabi',dx=t[-1]/nSteps*1e6,verbose=0)
pt.plot_t_rabi_data(x_vector=t,y_vector=data,exp_pars=qb.exp_pars,qb_pars=qb.qb_pars,fitted_pars=fitted_pars,device_name=device_name,project=project)

#%% Power Rabi
'''-----------------------------------------------------Power Rabi------------------------------------------------------'''

qb.exp_pars = {
    'n_avg':                512,
    'x0':                   0.01,
    'xmax':                 0.5,
    'dx':                   10e-3,
    'fsAWG':                2.4e9,
    'active_reset':         False,
    'qubit_reset_time':     300e-6,
    'qubit_drive_freq':     qb.qb_pars['qb_freq']+detuning,
    'tomographic-axis':     'Z',
    }

amp,data,nSteps = qb.pulsed_exp(exp='p-rabi',verbose=1,check_mixers=False)
# plot data
fitted_pars,error = pt.fit_data(x_vector=amp,y_vector=data,exp='p-rabi',dx=qb.dx,verbose=0)
pt.plot_p_rabi_data(x_vector=amp,y_vector=data,exp_pars=qb.exp_pars,qb_pars=qb.qb_pars,fitted_pars=fitted_pars,device_name=device_name,project=project)

qb.update_pi(pi_amp=fitted_pars[1]/2)


#%% Single Shot Experiment

qb.exp_pars = {
    'num_samples':          512,
    'n_avg':                1,
    'fsAWG':                2.4e9,
    'qubit_reset_time':     600e-6,
    }

data = qb.single_shot()

#make 2D histogram
pt.plot_single_shot(data, qb.exp_pars,qb.qb_pars,qb.iteration)

#%% T1 Measurement
qb.exp_pars = {
    'n_avg':                1024,
    'x0':                   100e-9,
    'xmax':                 250e-6,
    'dx':                   2e-6,
    'fsAWG':                0.3e9,
    'active_reset':         False,
    'qubit_drive_freq':     qb.qb_pars['qb_freq']+detuning,
    'tomographic-axis':     'Z',
}


t,data,nSteps = qb.pulsed_exp(exp='T1',verbose=1,check_mixers=False,save_data=True)
# plot data
fitted_pars,error = pt.fit_data(x_vector=t,y_vector=data,exp='T1',dx=t[-1]/nSteps,verbose=0)
pt.plot_T1_data(x_vector=t,y_vector=data,exp_pars=qb.exp_pars,qb_pars=qb.qb_pars,fitted_pars=fitted_pars,device_name=device_name,project=project)


#%% Ramsey Experiment

qb.exp_pars = {
    'n_avg':                512,
    'x0':                   60e-9,
    'xmax':                 120e-6,
    'dx':                   1000e-9,
    'fsAWG':                0.6e9,
    'amp_q':                0.1,
    'active_reset':         False,
    'qubit_drive_freq':     qb.qb_pars['qb_freq']+detuning,
    'tomographic-axis':     'Z',
}

t,data,nSteps = qb.pulsed_exp(exp='ramsey',verbose=1,check_mixers=False,save_data=True)

# fit & plot data
fitted_pars,error = pt.fit_data(x_vector=t,y_vector=data,exp='ramsey',dx=t[-1]/nSteps*1e6,verbose=0)
pt.plot_ramsey_data(x_vector=t,y_vector=data,exp_pars=qb.exp_pars,qb_pars=qb.qb_pars,fitted_pars=fitted_pars,device_name=device_name,project=project)

#%% tomography Experiment

qb.exp_pars = {
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

t,data2,nSteps = qb.pulsed_exp(exp='tomography',verbose=1,check_mixers=False,save_data=True)

# fit & plot data
fitted_pars,error = qb.fit_data(x_vector=t,y_vector=data,dx=t[-1]/nSteps*1e6,verbose=0)
qb.plot_ramsey_data(x_vector=t,y_vector=data,fitted_pars=fitted_pars)

#%% Virtual Z-Gate Experiment

qb.exp_pars = {
    'n_avg':                8192,
    'x0':                   0,
    'xmax':                 720,
    'dx':                   2,
    'fsAWG':                0.6e9,
    'active_reset':         False,
    'qubit_reset_time':     300e-6,
    'qubit_drive_freq':     qb.qb_pars['qb_freq']+detuning,
    'tomographic-axis':     'Z',
}


phase,data,nSteps = qb.pulsed_exp(exp='z-gate',verbose=1,check_mixers=False,save_data=True)
# fit & plot data
fitted_pars,error = qb.fit_data(x_vector=phase,y_vector=data,dx=phase[-1]/nSteps,verbose=0)
qb.plot_data(x_vector=phase,y_vector=data,fitted_pars=fitted_pars)

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


t,data,nSteps = qb.pulsed_exp(exp='echo',verbose=1,check_mixers=False,save_data=True)
fitted_pars,error = qb.fit_data(x_vector=t,y_vector=data,dx=t[-1]/nSteps*1e6,verbose=0)
qb.plot_echo_data(x_vector=t,y_vector=data,fitted_pars=fitted_pars)


#%% Mixer Optimization
qb.min_leak(inst=qb.awg,f_LO=qb.qb_pars['qb_LO'],mixer='qubit',mode='coarse',plot=True)
qb.min_leak(inst=qb.qa,f_LO=qb.qb_pars['rr_LO'],mixer='rr',mode='coarse',plot=True)

#f
qb.suppr_image(inst=qb.awg,f_LO=qb.qb_pars['qb_LO'],f_IF=qb.qb_pars['qb_IF']+detuning,mode='coarse',amp=0.2,threshold=-30)
qb.suppr_image(inst=qb.awg,f_LO=qb.qb_pars['qb_LO'],f_IF=qb.qb_pars['qb_IF']+detuning,mode='fine',amp=0.2,threshold=-60)

qb.suppr_image(inst=qb.qa,f_LO=qb.qb_pars['rr_LO'],f_IF=qb.qb_pars['rr_IF'],mode='coarse',amp=1)







