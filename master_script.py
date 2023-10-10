"""
Created on Thu Apr 14 11:35:33 2022

@author: Evangelos Vlachos <evlachos@usc.edu>
@edits: Malida Hecht <mohecht@usc.edu>
"""

#%% Initialization
#import qubit
from qubit import qubit
import numpy as np
import plotting as pt
#import matplotlib.pyplot as plt
from tqdm import tqdm
#from utilities import compute_bloch,compute_rho,compute_coherence,compute_purity, calc_steps

device_name = "WM1"
project = 'coherence-stabilization'
qb_name = 'qb5'
qb = qubit(qb_name)

#%% Spectroscopy (Resonator)
'''-----------------------------------------------------Resonator Spectroscopy------------------------------------------------------'''

freqs = np.arange(start=6.1155,stop=6.1165,step=5e-6) # frequencies are in GHz

qb.exp_pars = {
    'exp':                  'spectroscopy',
    'n_avg':                512,
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

freqs = np.arange(start=3.8,stop=4.5,step=1000e-6) # frequencies are in GHz

qb.exp_pars = {
    'exp':                  'spectroscopy',
    'n_avg':                2048,
    'element':              'qubit',
    'qubit_reset_time':     350e-6,
    'amp_q':                0.25,
    'satur_dur':            20e-6,
    'on_off':               True,
    'tomographic-axis':     'Z',
    'fsAWG':                2.4e9,
    'active_reset':         False,
    }

p_data,I,Q = qb.qb_spectroscopy(freqs,device_name)
pt.qb_spec_plot(freq=freqs,I=I,Q=Q,mag=p_data*1e3,exp_pars=qb.exp_pars,qb_pars=qb.qb_pars,find_pks=True,iteration = 1)



#%% Time Rabi
'''-----------------------------------------------------Time Rabi------------------------------------------------------'''
detuning = 0

qb.exp_pars = {
    'initial-state':        '0',
    'exp':                  't-rabi',
    'n_avg':                128,
    'x0':                   13e-9,
    'xmax':                 1e-6,
    'dx':                   6e-9,
    'fsAWG':                2.4e9,
    'amp_q':                0.3,
    'active_reset':         False,
    'qubit_drive_freq':     qb.qb_pars['qb_freq'],
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
    'xmax':                 0.5,
    'dx':                   5e-3,
    'fsAWG':                2.4e9,
    'active_reset':         False,
    'qubit_drive_freq':     qb.qb_pars['qb_freq'],
    'tomographic-axis':     'Z',
    }

amp,data,nSteps = qb.pulsed_exp(qb=qb_name,device_name=device_name,verbose=1,check_mixers=False)
# get fitted data

fitted_pars,error = pt.fit_data(x_vector=amp,y_vector=data,exp='p-rabi',dx=qb.dx,verbose=0)
#update threshold:
#if(qb.exp_pars['active_reset']):
#    pass
#else:
#    qb.update_qb_value('threshold',fitted_pars[2]*4096*1e-3)

# plot pulse
pt.plot_p_rabi_data(amp,data,fitted_pars,
                    qb=qb_name,
                    exp_pars=qb.exp_pars,
                    qb_pars=qb.qb_pars,
                    device_name=device_name,
                    project=project,
                    iteration =1)

qb.update_pi(pi_amp=fitted_pars[1]/2)

#%% Calibrate Rabi

qb.exp_pars = {
    'exp':                  'calibrate-rabi',
    'n_avg':                2048,
    'x0':                   0,
    'xmax':                 128,
    'dx':                   1,
    'fsAWG':                2.4e9,
    'active_reset':         False,
    'qubit_drive_freq':     qb.qb_pars['qb_freq'],
    'tomographic-axis':     'Z',
    }

amp,data,nSteps = qb.pulsed_exp(qb=qb_name,device_name=device_name,verbose=1,check_mixers=False)
# plot data



fitted_pars,error = pt.fit_data(x_vector=amp,y_vector=data,exp='p-rabi',dx=qb.dx,verbose=0)
#update threshold:
#if(qb.exp_pars['active_reset']):
#    pass
#else:
#    qb.update_qb_value('threshold',fitted_pars[2]*4096*1e-3)

pt.plot_p_rabi_data(amp,data,fitted_pars,qb=qb_name,exp_pars=qb.exp_pars,qb_pars=qb.qb_pars,device_name=device_name,project=project)

#qb.update_pi(pi_amp=fitted_pars[1]/2)

#%% Single Shot Experiment

qb.exp_pars = {
    'exp':                  'single-shot',
    'num_samples':          512,
    'n_avg':                1,
    'fsAWG':                2.4e9,
    'active_reset':         True,
    }

theta =0

data,offset_off,offset_pi = qb.single_shot(theta,device_name,qb_name)

#make 2D histogram
pt.plot_single_shot(data, qb.exp_pars,qb.qb_pars,qb.iteration)

#%% T1 Measurement

qb.wfm_pars = {
    't0':                   0.1e-6,
    'tmax':                 300e-6,
    'dt':                   4e-6,
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
    'qubit_drive_freq':     qb.qb_pars['qb_freq'],
    'tomographic-axis':     'Z',
}


t,data,nSteps = qb.pulsed_exp(qb=qb_name,device_name=device_name, verbose=1,check_mixers=False,save_data=False)
# plot data
fitted_pars,error = pt.fit_data(x_vector=t,y_vector=data,exp='T1',dx=t[-1]/nSteps,verbose=0)
pt.plot_T1_data(t,data,fitted_pars,qb=qb_name,exp_pars=qb.exp_pars,qb_pars=qb.qb_pars,device_name=device_name,project=project)


#%% Ramsey Experiment
qb.wfm_pars = {
    't0':                   0.1e-6,
    'tmax':                 100e-6,
    'dt':                   .025e-6,
    'fsAWG':                0.6e9,
    'mu':                   0,
    'sigma':                0e-3,
    }
qb.exp_pars = {
    'exp':                  'ramsey',
    'n_avg':                512,
    'x0':                   100e-9,
    'xmax':                 300e-6,
    'dx':                   1500e-9,
    'fsAWG':                0.6e9,
    'amp_q':                0.4,
    'active_reset':         False,
    'qubit_drive_freq':     qb.qb_pars['qb_freq'],
    'tomographic-axis':     'Z',
}

t,data,nSteps = qb.pulsed_exp(qb=qb_name,device_name=device_name,verbose=1,check_mixers=False,save_data=True)

# fit & plot data
fitted_pars,error = pt.fit_data(x_vector=t,y_vector=data,exp='ramsey',dx=t[-1]/nSteps*1e6,verbose=0)
pt.plot_ramsey_data(t,data,
                    fitted_pars,
                    qb=qb_name,
                    exp_pars=qb.exp_pars,
                    qb_pars=qb.qb_pars,
                    device_name=device_name,
                    project=project,
                    iteration =1)
# qb.update_qb_value('qb_freq', qb.exp_pars['qubit_drive_freq']-fitted_pars[1]*1e6)

## FOR RAMSEY TO UPDATE QB VALUE:

#qb.update_qb_value('qb_freq',qb.qb_pars['qb_freq']+0.002e6)
#%% Tomography Experiment
''' this one isnt working right'''
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
qb.plot_ramsey_data(t,data,fitted_pars,
                    qb=qb_name,
                    exp_pars=qb.exp_pars,
                    qb_pars=qb.qb_pars,
                    device_name=device_name,
                    project=project,
                    iteration = qb.iteration)



#%% Coordinate System Calibration
'''Execute time-based rabi experiment, and measure along 3 different axis at the end'''
detuning = 0

qb.exp_pars = {
    'exp':                  'tomography-calibration',
    'n_avg':                2048,
    'x0':                   13e-9,
    'xmax':                 0.8e-6,
    'dx':                   6e-9,
    'fsAWG':                2.4e9,
    'amp_q':                0.15,
    'active_reset':         False,
    'qubit_drive_freq':     qb.qb_pars['qb_freq']+detuning,
}

t,data_cal,nSteps = qb.tomography_calibration(qb=qb_name,device_name=device_name,verbose=1,check_mixers=False,save_data=False)
# fit & plot data
fitted_pars,error = pt.fit_data(x_vector=t,y_vector=data_cal[2,:],exp='t-rabi',dx=t[-1]/nSteps*1e6,verbose=0)
calib_states = qb.cal_coord(fitted_pars[-1]*1e-3,fitted_pars[0]*1e-3)
pt.tom_calib_plot(x_data=t, y_data=data_cal, coords=calib_states)

#%% coherence stabilization Experiment
qb.wfm_pars = {
    'x0':                   50e-9,
    'xmax':                 .5e-6,
    'dx':                   0.05e-6,
    'fsAWG':                1.2e9,
    'mu':                   0,
    'sigma':                211e-3,
    'T2':                   6.6e-6,         
    }


qb.exp_pars = {
    'exp':                  'coherence-stabilization',
    'initial-state':        (1/2,1/5),  # in units of pi (azimuthal,polar)
    'n_avg':                1024,
    'n_realizations':       1,
    'x0':                   qb.wfm_pars['x0'],
    'dx':                   qb.wfm_pars['dx'],
    'fsAWG':                qb.wfm_pars['fsAWG'],
    'amp_q':                1,
    'active_reset':         True,
    'qubit_drive_freq':     qb.qb_pars['qb_freq'],
}

wfms,data,v_b,calib_states,nSteps = qb.coherence_stabilization(qb=qb_name,device_name=device_name,verbose=1,save_data=True)

t = np.linspace(qb.wfm_pars['x0'],qb.wfm_pars['xmax'],data.shape[1]-2)
# fit & plot data
# v0 = compute_bloch(data[:,0], calib_states)
# vf = compute_bloch(data[:,-1], calib_states)
# rho_0 = compute_rho(v0)
# rho_f = compute_rho(vf)


# pt.tom_calib_plot(x_data=t, y_data=data, coords=calib_states,data='data')
# pt.plot_tomography(rho_0,qb.exp_pars['initial-state'], tmax=t[0]*1e6,cal_states=calib_states)
# pt.plot_tomography(rho_f,qb.exp_pars['initial-state'], tmax=t[-1]*1e6,cal_states=calib_states)

pt.plot_coherence(t,v_b,wfms,qb.exp_pars,qb.qb_pars,qb.wfm_pars,calib_states)

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
fitted_pars,error = pt.fit_data(x_vector=phase,y_vector=data,dx=phase[-1]/nSteps,verbose=0)
pt.plot_data(phase,data,fitted_pars,qb=qb_name)

#%% Echo Measurement
detuning = 0
qb.exp_pars = {
    'exp':                  'echo',
    'n_avg':                1024,
    'x0':                   50e-9,
    'xmax':                 300e-6,
    'dx':                   4000e-9,
    'fsAWG':                0.6e9,
    'amp_q':                0.1,
    'active_reset':         True,
    'qubit_reset_time':     300e-6,
    'qubit_drive_freq':     qb.qb_pars['qb_freq']+detuning,
    'tomographic-axis':     'Z',
}


t,data,nSteps = qb.pulsed_exp(qb=qb_name,device_name=device_name,verbose=1,check_mixers=False,save_data=True)
fitted_pars,error = pt.fit_data(x_vector=t,y_vector=data,exp = 'echo',dx=t[-1]/nSteps*1e6,verbose=0)
pt.plot_echo_data(x_vector=t,y_vector=data,fitted_pars=fitted_pars,exp_pars = qb.exp_pars,qb_pars = qb.qb_pars,iteration =1)

#%% Mixer Optimization
qb.min_leak(inst=qb.awg,f_LO=qb.qb_pars['qb_LO'],mixer='qubit',cal='lo',plot=True,span=0.4e6)
 
qb.min_leak(inst=qb.qa,f_LO=qb.qb_pars['rr_LO'],mixer='rr',cal='lo',plot=True,span=0.1e6)

# qb.min_leak(inst=qb.awg,f_LO=qb.qb_pars['qb_LO'],f_IF=qb.qb_pars['qb_IF']+detuning,cal='ssb',amp=0.4,threshold=-60,span=0.1e6)

qb.suppr_image(inst=qb.awg,f_LO=qb.qb_pars['qb_LO'],f_IF=qb.qb_pars['qb_IF'],amp=0.3)
 

## if it screams, run t1 (if need more do rabi/ramsey) to save parameters that arent in the JSON :]



#%% New Coherence Stabilization:
    
qb.wfm_pars = {
    'x0':                   50e-9,
    'xmax':                 5e-6,
    'dx':                   0.25e-6,
    'fsAWG':                1.2e9,
    'mu':                   0,
    'sigma':                0e-3,
    'T2':                   3e-6,         
    }


qb.exp_pars = {
    'exp':                  'coherence-stabilization',
    'initial-state':        (1/2,1/4),  # in units of pi (azimuthal,polar)
    'n_avg':                1024,
    'n_realizations':       1,
    'x0':                   qb.wfm_pars['x0'],
    'dx':                   qb.wfm_pars['dx'],
    'fsAWG':                qb.wfm_pars['fsAWG'],
    'amp_q':                1,
    'active_reset':         False,
    'qubit_drive_freq':     qb.qb_pars['qb_freq']
}

# initiate data and v_b 

#print('Initizalizing...')
data, v_b, timeout, n_steps, filepath = qb.setup_coherence_stabilization(device_name=device_name, 
                                    qb=qb_name, 
                                    calibrate=False, 
                                    verbose=True)
#wfms = []
t = np.linspace(qb.wfm_pars['x0'],qb.wfm_pars['xmax'],v_b.shape[1])

#print('Initizalization Complete, starting measurements...')
with tqdm(total=qb.exp_pars['n_realizations']) as pbar:
    for ii in range(qb.exp_pars['n_realizations']):
        wfms,data[:,ii],v_b[:,:,ii],calib_states,nSteps = qb.coherence_stabilization_single(timeout=timeout,
                                       n_steps=n_steps,
                                       filepath=filepath,
                                       device_name=device_name, 
                                       qb = qb_name, 
                                       verbose=True)
        
        
        pbar.update(1)
        print('\n')
# compute average for bloch vectors
v_ave = np.mean(v_b, 2)

# compute x-axis


pt.plot_coherence(t,v_ave,wfms,
                  exp_pars=qb.exp_pars,
                  qb_pars = qb.qb_pars,
                  wfm_pars = qb.wfm_pars,
                  calib_states=calib_states, 
                  iteration =qb.iteration)
