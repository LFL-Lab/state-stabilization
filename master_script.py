"""
Created on Thu Apr 14 11:35:33 2022

@author: Evangelos Vlachos <evlachos@usc.edu>
@edits: Malida Hecht <mohecht@usc.edu>
"""

#%% Initialization
from qubit import qubit
import numpy as np

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
    'rr_atten':             25,
    'on_off':               True,
    }

p_data,I,Q = qb.rr_spectroscopy(freqs)
fc = qb.rr_spec_plot(freq=freqs,I=I,Q=Q,mag=p_data,df=1e9*(freqs[1]-freqs[0]),find_pks=True)
qb.update_qb_value('rr_LO',fc*1e9)

#%% Spectroscopy (Qubit)

'''-----------------------------------------------------Qubit Spectroscopy------------------------------------------------------'''

freqs = np.arange(start=3.9,stop=4.2,step=200e-6) # frequencies are in GHz

qb.exp_pars = {
    'n_avg':                1024,
    'element':              'qubit',
    'qubit_reset_time':     200e-6,
    'amp_q':                0.25,
    'satur_dur':            40e-6,
    'on_off':               True,
    }

p_data,I,Q = qb.qb_spectroscopy(freqs)
qb.qb_spec_plot(freq=freqs,I=I,Q=Q,mag=p_data*1e3,find_pks=True)

## Need to add title to this plot ^^


#%% Time Rabi
'''-----------------------------------------------------Time Rabi------------------------------------------------------'''

qb.exp_pars = {
    'n_avg':                512,
    'x0':                   13e-9,
    'xmax':                 1e-6,
    'dx':                   6e-9,
    'fsAWG':                2.4e9,
    'amp_q':                0.3,
    'active_reset':         False,
    'qubit_reset_time':     300e-6,
    'qubit_drive_freq':     3.879e9+0.61e6,
    }

t,data,nSteps = qb.pulsed_exp(exp='rabi',verbose=1,check_mixers=False)
# plot data
fitted_pars,error = qb.fit_data(x_vector=t,y_vector=data,dx=t[-1]/nSteps*1e6,verbose=0)
qb.plot_data(x_vector=t,y_vector=data,fitted_pars=fitted_pars)

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
    'qubit_drive_freq':     3.879e9+0.61e6,
    }

amp,data,nSteps = qb.pulsed_exp(exp='p-rabi',verbose=1,check_mixers=False)
# plot data
fitted_pars,error = qb.fit_data(x_vector=amp,y_vector=data,dx=qb.dx,verbose=0)
qb.plot_data(x_vector=amp,y_vector=data,fitted_pars=fitted_pars)

qb.update_pi(pi_amp=fitted_pars[1]/2)


#%% Single Shot Experiment

qb.exp_pars = {
    'num_samples':          512,
    'n_avg':                1,
    'fsAWG':                2.4e9,
    'qubit_reset_time':     300e-6,
    }

data = qb.single_shot()

#make 2D histogram
qb.plot_single_shot(data)

#%% T1 Measurement
qb.exp_pars = {
    'n_avg':                512,
    'x0':                   100e-9,
    'xmax':                 250e-6,
    'dx':                   2e-6,
    'fsAWG':                0.3e9,
    'active_reset':         False,
    'qubit_reset_time':     300e-6,
    'qubit_drive_freq':     3.879e9+0.61e6,
}


t,data,nSteps = qb.pulsed_exp(exp='T1',verbose=1,check_mixers=False,save_data=True)
# plot data
fitted_pars,error = qb.fit_data(x_vector=t,y_vector=data,dx=t[-1]/nSteps,verbose=0)
qb.plot_data(x_vector=t,y_vector=data,fitted_pars=fitted_pars)


#%% Ramsey Experiment

qb.exp_pars = {
    'n_avg':                512,
    'x0':                   50e-9,
    'xmax':                 200e-6,
    'dx':                   2e-6,
    'fsAWG':                0.6e9,
    'amp_q':                0.1,
    'active_reset':         False,
    'qubit_reset_time':     300e-6,
    'qubit_drive_freq':     3.879e9+0.61e6,
}


t,data,nSteps = qb.pulsed_exp(exp='ramsey',verbose=1,check_mixers=False,save_data=True)

# fit & plot data
fitted_pars,error = qb.fit_data(x_vector=t,y_vector=data,dx=t[-1]/nSteps*1e6,verbose=0)
qb.plot_data(x_vector=t,y_vector=data,fitted_pars=fitted_pars)


#%% Virtual Z-Gate Experiment

qb.exp_pars = {
    'n_avg':                8192,
    'x0':                   0,
    'xmax':                 720,
    'dx':                   2,
    'fsAWG':                0.6e9,
    'active_reset':         False,
    'qubit_reset_time':     300e-6,
    'qubit_drive_freq':     3.879e9+0.61e6,
}


phase,data,nSteps = qb.pulsed_exp(exp='z-gate',verbose=1,check_mixers=False,save_data=True)
# fit & plot data
fitted_pars,error = qb.fit_data(x_vector=phase,y_vector=data,dx=phase[-1]/nSteps,verbose=0)
qb.plot_data(x_vector=phase,y_vector=data,fitted_pars=fitted_pars)

#%% Echo Measurement

qb.exp_pars = {
    'n_avg':                512,
    'x0':                   100e-9,
    'xmax':                 100e-6,
    'dx':                   2e-6,
    'fsAWG':                0.6e9,
    'amp_q':                0.1,
    'active_reset':         False,
    'qubit_reset_time':     200e-6,
    'qubit_drive_freq':     3.879e9+0.61e6,
}


t,data,nSteps = qb.pulsed_exp(exp='echo',verbose=1,check_mixers=False,save_data=True)
fitted_pars,error = qb.fit_data(x_vector=t,y_vector=data,dx=t[-1]/nSteps*1e6,verbose=0)
qb.plot_data(x_vector=t,y_vector=data,fitted_pars=fitted_pars)


#%% Mixer Optimization
qb.min_leak(inst=qb.awg,f_LO=qb.qb_pars['qb_LO'],mixer='qubit',mode='coarse',plot=True)
qb.min_leak(inst=qb.qa,f_LO=qb.qb_pars['rr_LO'],mixer='rr',mode='coarse',plot=True)

#f
qb.suppr_image(inst=qb.awg,f_LO=qb.qb_pars['qb_LO'],f_IF=qb.qb_pars['qb_IF'],mode='coarse')
qb.suppr_image(inst=qb.qa,f_LO=qb.qb_pars['rr_LO'],f_IF=qb.qb_pars['rr_IF'],mode='coarse',amp=1)







