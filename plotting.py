# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 16:30:32 2023

@author: Evangelos Vlachos <evlachos@usc.edu>
"""

import numpy as np
from scipy.signal import butter,find_peaks
import seaborn as sns; sns.set() # styling
import pandas as pd
import scipy as scy
import qutip as qt
import matplotlib.pyplot as plt
import lmfit
from utilities import normalize_data,compute_bloch,compute_coherence,compute_purity
pi = np.pi


''' Plotting Parameters '''
# plt.style.use(['science','no-latex'])
# plt.rcParams['figure.dpi'] = 150
# plt.rcParams['axes.facecolor'] = 'white'
# plt.rcParams['axes.edgecolor'] = 'black'
# plt.rcParams['axes.grid'] = False
# plt.rcParams['figure.frameon'] = True
# plt.rcParams.update(plt.rcParamsDefault)
# for a complete set of parameters "print(plt.rcParams)"
sns.set_style('ticks')
plt.rcParams['figure.dpi'] = 150
plt.rcParams['text.usetex'] = False
plt.rcParams['font.family'] =  'Arial'
plt.rcParams['font.size'] = 16
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["xtick.major.top"] = True
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.titlesize'] = 16
plt.rcParams["xtick.major.bottom"] = True
plt.rcParams["xtick.top"] = True
plt.rcParams["xtick.bottom"] = True
plt.rcParams["ytick.left"] = True
plt.rcParams["ytick.right"] = True
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams["ytick.major.right"] = True
plt.rcParams["ytick.labelright"] = False
plt.rcParams["ytick.minor.visible"] = False


#%% spectroscopy_plots

def qb_spec_plot(freq,I,Q,mag,exp_pars={},qb_pars={},attenuation=-30,iteration=1,find_pks=False):

    I = I*1e3
    Q = Q*1e3

    phase = np.unwrap(np.angle(I+1j*Q))
    if find_pks:
        sigma = np.std(mag)
        print(f'Peak threshold at {np.mean(mag)+3*sigma:.1f}')
        peaks,_ = find_peaks(mag,height=np.mean(mag)+3*sigma,distance=200,width=3)
        try:
            for i in peaks:
                print(f'Peaks at: {round(freq[i],5)} GHz\n')
        except:
            print('Peaks not found or do not exist.')

    fig = plt.figure(figsize=(5,4))

    # # I data
    # ax1 = fig.add_subplot(221)
    # ax1.plot(freq,I,'-o', markersize = 3, c='C0')
    # ax1.set_xlabel('Frequency (GHz)')
    # ax1.set_ylabel('I (mV)')
    # # Q data
    # ax1 = fig.add_subplot(222)
    # ax1.plot(freq,Q,'-o', markersize = 3, c='C0')
    # ax1.set_xlabel('Frequency (GHz)')
    # ax1.set_ylabel('Q (mV)')
    # Power data
    ax1 = fig.add_subplot()
    ax1.plot(freq,mag,'-o', markersize = 3, c='C0')
    ax1.set_xlabel('Frequency (GHz)')
    ax1.set_ylabel('Magnitude (mV)')
    # phase data
    # ax1 = fig.add_subplot(212)
    # ax1.plot(freq,phase,'-o', markersize = 3, c='C0')
    # ax1.set_xlabel('Frequency (GHz)')
    # ax1.set_ylabel('Phase (rad)')

    if len(peaks) == 2:
        txt = '$\omega_{01}$ = %.4f GHz\n$\omega_{02}/2$ = %.4f GHz\n$\\alpha$ = %.1f MHz\n$A_{qb}$ = %.2f V\nReadout Attn = -%d dB\n$f_r$ = %.4f GHz'%(freq[peaks[1]],freq[peaks[0]],2*(freq[peaks[0]]-freq[peaks[1]])*1e3,exp_pars['amp_q'],qb_pars['rr_atten'],qb_pars['rr_LO']*1e-9)
    elif len(peaks) == 1:
        txt = '$\omega_{01}$ = %.4f GHz\n$A_{qb}$ = %.1f mV\nReadout Attn = -%d dB\n$f_r$ = %.4f GHz'%(freq[peaks[0]],exp_pars['amp_q']*1e3,qb_pars['rr_atten'],qb_pars['rr_LO']*1e-9)
    else:
        txt = '$A_{qb}$ = %.1f mV\nReadout Attn = -%d dB\n$f_r$ = %.4f GHz'%(exp_pars['amp_q']*1e3,qb_pars['rr_atten'],qb_pars['rr_LO']*1e-9)
    plt.gcf().text(1, 0.15, txt, fontsize=14)
    ax1.set_title(f'{exp_pars["element"]} spectroscopy {iteration}')
    plt.tight_layout()
    plt.show()
    
def rr_spec_plot(freq,I,Q,mag,exp_pars={},qb_pars={},df=0.1e6,iteration=1,find_pks=False):

    I = I*1e3
    Q = Q*1e3
    freq = freq[1:]
    mag = mag[1:]*1e3
    
    # if find_pks:
    #     sigma = np.std(mag)
    #     print(f'Peak threshold at {np.mean(mag)+3*sigma}')
    #     peaks,_ = find_peaks(mag,height=np.mean(mag)+3*sigma,distance=200,width=3)
    #     try:
    #         for i in peaks:
    #             print(f'Peaks at: {round(freq[i],5)} GHz\n')
    #     except:
    #         print('Peaks not found or do not exist.')
            
    fc,fwhm = fit_res(freq,mag)

    fig = plt.figure(figsize=(5,4))

    # Power data
    ax1 = fig.add_subplot(111)
    ax1.plot(freq,mag,'-o', markersize = 3, c='C0')
    ax1.set_xlabel('Frequency (GHz)')
    ax1.set_ylabel('Magnitude (mV)')
    # phase = np.unwrap(np.angle(I+1j*Q))
    # ax2 = fig.add_subplot(212)
    # ax2.plot(freq,phase,'-o', markersize = 3, c='C0')
    # ax2.set_xlabel('Frequency (GHz)')
    # ax2.set_ylabel('Phase (deg)')

    txt = f'$\omega_c$ = {fc:.6f} GHz\nFWHM = {fwhm*1e6:.1f} kHz\n$T_1$ = {1e6/(np.pi*fwhm*1e6):.1f} ns\n$P_r$ = -{exp_pars["rr_atten"]:.1f} dB\n$df$ = {df*1e-3:.1f} kHz'
    plt.gcf().text(1, 0.15, txt, fontsize=14)
    ax1.set_title(f'{exp_pars["element"]} spectroscopy {iteration}')
    plt.tight_layout()
    plt.show()
    
    return fc
    
def heatplot(xdata, ydata, z_data, xlabel = "", ylabel = "", normalize=False, 
             cbar_label = 'log mag',title='', **kwargs):
    
    fig = plt.figure(figsize=(6,5), dpi=300)
    ax1  = fig.add_subplot(6,1,(1,4))
    # if normalize:
    #     cbar_label += ' (normalized)'

    df = pd.DataFrame(z_data, columns = xdata, index = ydata)

    if normalize:
        # df = df.apply(lambda x: (x-x.mean())/x.std(), axis = 1)
        df = df.apply(lambda x: (x/x.max()), axis = 1)
    
    cbar_options = {
        'label':    cbar_label,
        # 'ticks':    np.around(np.linspace(np.amin(z_data),np.amax(z_data),5),1),
        'pad':      0.05,
        # 'values':   np.linspace(np.amin(z_data),np.amax(z_data),1000),
        'shrink':   1.1,
        'location': 'top',
    }

    kwargs = {
        'linewidths':  0,
        # 'xticklabels': np.linspace(min(xdata),max(xdata)+0.5,5),
        # 'yticklabels': np.linspace(min(ydata),max(ydata)+0.5,5),
        'vmin':        np.amin(z_data),
        'vmax':        np.amax(z_data)
    }
    

    hm = sns.heatmap(df, ax=ax1,cmap = 'seismic', cbar_kws=cbar_options)
    # hm.set_xlabel(xlabel, fontsize=12)
    hm.set_ylabel(ylabel, fontsize=14)
    hm.spines[:].set_visible(True)
    # ax1.set_title(title,fontsize=12)
    ax1.tick_params(direction='out',length=0.01,width=0.5,bottom=False, top=False, left=True, right=False,labeltop=False, labelbottom=True,labelrotation=90,labelsize=10,size=8)
    plt.yticks(rotation=0)
    
    ax2 = fig.add_subplot(6,1,(5,6))
    ax2.plot(xdata,z_data[0,:]/max(z_data[0,:]),'o', markersize = 3, c='b',label=f'$P_r$ = -{ydata[0]} dB')
    ax2.plot(xdata,z_data[-1,:]/max(z_data[-1,:]),'o', markersize = 3, c='r',label=f'$P_r$ = -{ydata[-1]} dB')
    ax2.legend()
    ax2.set_xlabel(xlabel, fontsize=12)
    ax2.set_ylabel(cbar_label,fontsize=12)
    fc1,_ = fit_res(xdata, z_data[0,:])
    fc2,_ = fit_res(xdata, z_data[-1,:])
    txt = f'$f_1$ = {fc1:.5f} GHz\n$f_2$ = {fc2:.5f} GHz\n$2\chi/2\pi$ = {(fc2-fc1)*1e3:.1f} MHz'
    plt.gcf().text(0.95, 0.15, txt, fontsize=14)
    
    # plt.tight_layout()
        
    return df
    
#%% time-based_plots
def plot_p_rabi_data(x_vector,y_vector,fitted_pars,qb='',exp_pars={},qb_pars={},iteration=1,device_name='',project='',savefig=True):

    qb_drive_freq = exp_pars['qubit_drive_freq']*1e-9
    fig, ax = plt.subplots()
    y_vector = y_vector*1e3
    
    ax.plot(x_vector, y_vector, '-o', markersize = 3, c='C0')
    ax.set_ylabel('Digitizer Voltage (mV)')
    ax.set_xlabel('Pulse Amplitude')
    ax.plot(x_vector,p_rabi(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2]),'r')
    ax.set_title(f'Rabi Measurement {iteration:03d}')
    textstr = f'$\omega_d$ = {qb_drive_freq:.4f} GHz\n$N$ = {exp_pars["n_avg"]}\n$A_\pi$ = {fitted_pars[1]/2:.3f}\n$T_\pi$ = {qb_pars["gauss_len"]/2.4:.1f} ns'
    
    plt.gcf().text(0.95, 0.15, textstr, fontsize=14)

    plt.tick_params(axis='both',direction='in',bottom=True, top=True, left=True, right=True,size=8)

    if savefig:
        plt.savefig(f'D:\\{project}\\{device_name}\\{qb}\\p-rabi-data\\fig_{iteration:03d}.png',dpi='figure')


def plot_t_rabi_data(x_vector,y_vector,fitted_pars,qb='',exp_pars={},qb_pars={},device_name='',project='',iteration=1,savefig=True):
    
    qb_drive_freq = exp_pars['qubit_drive_freq']*1e-9
    fig, ax = plt.subplots()
    y_vector = y_vector*1e3
    x_vector = x_vector*1e9
    
    ax.plot(x_vector, y_vector, '-o', markersize = 3, c='C0')
    ax.set_ylabel('Digitizer Voltage (mV)')
    ax.set_xlabel('Pulse Duration (ns)')
    ax.plot(x_vector,rabi(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2],fitted_pars[3]),'r')
    ax.set_title(f'Rabi Measurement {iteration:03d}')
    textstr = f'$\omega_d$ =  {qb_drive_freq:.4f} GHz\n$A_q$ = {exp_pars["amp_q"]:.2f} V\n$N$ = {exp_pars["n_avg"]}'
    
    plt.gcf().text(0.95, 0.15, textstr, fontsize=14)

    plt.tick_params(axis='both',direction='in',bottom=True, top=True, left=True, right=True,size=8)

    if savefig:
        plt.savefig(f'D:\\{project}\\{device_name}\\{qb}\\t-rabi-data\\fig_{iteration:03d}.png',dpi='figure')

def plot_ramsey_data(x_vector,y_vector,fitted_pars,qb='',exp_pars={},qb_pars={},fitFunc='',iteration=1,device_name='',project='',savefig=True):
     
    qb_drive_freq = exp_pars['qubit_drive_freq']*1e-9
    fig, ax = plt.subplots()
    x_vector = x_vector*1e6
    y_vector = y_vector*1e3
    
    ax.set_xlabel('Pulse Separation ($\mu$s)')
    ax.plot(x_vector, y_vector, '-o', markersize = 3, c='C0')
    
    if fitFunc == 'envelope':
        ax.plot(x_vector,decay(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2]),'r')
    else:
        ax.plot(x_vector,ramsey(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2],fitted_pars[3],fitted_pars[4]),'r')
    ax.set_title(f'Ramsey Measurement {iteration:03d}')
    textstr = f'$T_\pi$ = {qb_pars["pi_len"]/2.4:.1f} ns\n$A_\pi$ = {qb_pars["pi_amp"]*1e3:.1f} mV\n$\omega_d$ = {qb_drive_freq:.4f} GHz\n$\Delta$ = {fitted_pars[1]:.3f} MHz\n$T_2^R$ = {fitted_pars[3]:.1f} $\mu$s\n$N$ = {exp_pars["n_avg"]}'

    plt.gcf().text(0.95, 0.15, textstr, fontsize=14)

    plt.tick_params(axis='both',direction='in',bottom=True, top=True, left=True, right=True,size=8)

    if savefig:
        plt.savefig(f'D:\\{project}\\{device_name}\\{qb}\\ramsey-data\\fig_{iteration:03d}.png',dpi='figure')


def plot_z_gate_data(x_vector,y_vector,fitted_pars,qb='',exp_pars={},qb_pars={},iteration=1,device_name='',project='',savefig=True):
    
    qb_drive_freq = exp_pars['qubit_drive_freq']*1e-9
    fig, ax = plt.subplots()
    y_vector = y_vector*1e3
    
    ax.set_ylabel('Digitizer Voltage (mV)')
    ax.set_xlabel('$\phi$ (deg)')
    ax.plot(x_vector, y_vector, '-o', markersize = 3, c='C0')
    ax.plot(x_vector,cos(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2],fitted_pars[3]),'r')
    ax.set_title(f'Virtual Z-Gate Measurement {iteration:03d}')
    textstr = f'$T_\pi$ = {qb_pars["pi_len"]/2.4:.1f} ns\n$A_\pi$ = {qb_pars["pi_amp"]*1e3:.1f} mV\n$\omega_d$ = {qb_drive_freq:.4f} GHz\n$N$ = {exp_pars["n_avg"]}'

    plt.gcf().text(0.95, 0.15, textstr, fontsize=14)

    plt.tick_params(axis='both',direction='in',bottom=True, top=True, left=True, right=True,size=8)

    if savefig:
        plt.savefig(f'D:\\{project}\\{device_name}\\{qb}\\z-gate-data\\fig_{iteration:03d}.png',dpi='figure')

def plot_echo_data(x_vector,y_vector,fitted_pars,qb='',exp_pars={},qb_pars={},iteration=1,device_name='',project='',savefig=True):
     
    qb_drive_freq = exp_pars['qubit_drive_freq']*1e-9
    fig, ax = plt.subplots()
    x_vector = x_vector*1e6
    y_vector = y_vector*1e3
    
    ax.plot(x_vector, y_vector, '-o', markersize = 3, c='C0')
    ax.set_ylabel('Digitizer Voltage (mV)')
    ax.set_xlabel('Pulse Separation ($\mu$s)')
    ax.plot(x_vector,decay(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2]),'r')
    textstr = f'$T_\pi$ = {qb_pars["pi_len"]/2.4:.1f} ns\n$A_\pi$ = {qb_pars["pi_amp"]*1e3:.1f} mV\n$\omega_d$ = {qb_drive_freq:.4f} GHz\n$T_2^E$ = {fitted_pars[1]:.1f} $\mu$s\n$N$ = {exp_pars["n_avg"]}'
    ax.set_title(f'Echo Measurement {iteration:03d}')
    
    plt.gcf().text(0.95, 0.15, textstr, fontsize=14)

    plt.tick_params(axis='both',direction='in',bottom=True, top=True, left=True, right=True,size=8)

    if savefig:
        plt.savefig(f'D:\\{project}\\{device_name}\\{qb}\\echo-data\\fig_{iteration:03d}.png',dpi='figure')


def plot_T1_data(x_vector,y_vector,fitted_pars,qb='',exp_pars={},qb_pars={},iteration=1,device_name='',project='',savefig=True):
    
    qb_drive_freq = exp_pars['qubit_drive_freq']*1e-9
    fig, ax = plt.subplots()
    x_vector = x_vector*1e6
    y_vector = y_vector*1e3
    
    ax.plot(x_vector, y_vector, '-o', markersize = 3, c='C0')
    ax.set_ylabel('Digitizer Voltage (mV)')
    ax.set_xlabel('Delay ($\mu$s)')
    ax.plot(x_vector,decay(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2]),'r')
    textstr = f'$T_\pi$ = {qb_pars["pi_len"]/2.4:.1f} ns\n$A_\pi$ = {qb_pars["pi_amp"]*1e3:.1f} mV\n$\omega_d$ = {qb_drive_freq:.4f} GHz\n$T_1$ = {fitted_pars[1]:.1f} $\mu$s\n$N$ = {exp_pars["n_avg"]}'
    ax.set_title(f'T1 Measurement {iteration:03d}')

    plt.gcf().text(0.95, 0.15, textstr, fontsize=14)

    plt.tick_params(axis='both',direction='in',bottom=True, top=True, left=True, right=True,size=8)

    if savefig:
        plt.savefig(f'D:\\{project}\\{device_name}\\{qb}\\T1-data\\fig_{iteration:03d}.png',dpi='figure')
        
        
#%% tomography_plots
def plot_tomography(data,initial_state,tmax,cal_states):
    
    tick_loc = [1,2]
    yticklabels = [r'$|1\rangle$',r'$|0\rangle$']
    xticklabels = [r'$\langle1|$',r'$\langle0|$']
    ztick_loc = [-1,-0.5,0,0.5,1]
    states = [['$0$','$1$']]
    fig,[ax1,ax2] = qt.qpt_plot(data,states,title='Quantum State Tomography\n'+r'$\rho_0$'+f'= {initial_state[0]} | T = {tmax:.2f}' +r'$\mu$s') 
    ax1.set_xticks(ticks=tick_loc,labels=xticklabels)
    ax1.set_yticks(ticks=tick_loc,labels=yticklabels)
    ax1.set_zticks(ticks=ztick_loc)
    ax2.set_xticks(ticks=tick_loc,labels=xticklabels)
    ax2.set_yticks(ticks=tick_loc,labels=yticklabels)
    ax2.set_zticks(ticks=ztick_loc)
    
    # plt.xlim(0,2)
    # plt.ylim(0,2)
    # plt.tight_layout()
    # ax1.set_zlabel(r'$|{\rho}|$')
    return ax1

def tom_calib_plot(x_data,y_data,coords,data='cal'):
    
    x_data = x_data*1e6
    plot_data = np.zeros(y_data.shape)
    # plot_data = y_data
    
    if data == 'cal':
        # y_data = y_data*1e3
        plot_data = normalize_data(y_data)
        labels = [r'$\langle X|\psi\rangle$',r'$\langle Y|\psi\rangle$',r'$\langle Z|\psi\rangle$']
    else:
        for i in range(y_data.shape[1]):
            plot_data[:,i] = compute_bloch(y_data[:,i],calib_pars=coords) 
        labels = ['$v_x$','$v_y$','$v_z$']
    
    fig, ax = plt.subplots(figsize=(6,4),dpi=150)
    
    ax.plot(x_data,plot_data[0,:],'-o',color='b',label=labels[0])
    ax.plot(x_data,plot_data[1,:],'-x',color='r',label=labels[1])
    ax.plot(x_data,plot_data[2,:],'-<',color='k',label=labels[2])
    
    ax.set_xlabel('Drive Duration ($\mu$s)')
    # ax.set_ylabel('Bloch Vector')
    ax.set_ylim(-1,1)
    
    txt = r"$|+X\rangle$"+f" = {coords[0]*1e3:.1f} mV\n"+r"$|+Y\rangle$" +f"= {coords[1]*1e3:.1f} mV\n"+r"$|1\rangle$"+f" = {coords[2]*1e3:.1f} mV\n"
    plt.gcf().text(0.95, 0.15, txt, fontsize=14)
    ax.legend(loc='upper right')

def plot_bloch(state):

    b = qt.Bloch()
    b.add_vectors(state)
    b.show()
     
#%% mixer_opt_plots
def power_plot(freqs,signal,power,fc):
    plt.plot(freqs*1e-6, signal,'-')
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Power [dBm]')
    plt.show()
    
def plot_mixer_opt(par1,par2,power_data,cal='LO',mixer='qubit',fc=5e9):
    
    if cal == 'LO':
        par1 = np.around(par1*1e3,1)
        par2 = np.around(par2*1e3,1)
    else:
        par1 = np.around(par1,3)
        par2 = np.around(par2,3)
    par1 = par1.tolist()
    par2 = par2.tolist()
    df = pd.DataFrame(data=power_data,index=par1,columns=par2)

    hm = sns.heatmap(df,cbar_kws={'label': "Power [dBm]"})

    if cal == 'LO':
        hm.set_ylabel('I [mV]')
        hm.set_xlabel('Q [mV]')
    elif cal == 'SB':
        hm.set_ylabel('Gain Imbalance (%)')
        hm.set_xlabel('Phase Imbalance')

    hm.spines[:].set_visible(True)
    hm.tick_params(direction='out',length=0.01,width=0.5,bottom=True, top=False, left=True, right=True,labeltop=False, labelbottom=True,labelrotation=90,labelsize=10,size=10)
    plt.yticks(rotation=0)
    plt.tight_layout()
    if mixer == 'qubit':
        plt.title(f'Qubit Mixer {cal} Calibration at {round(fc*1e-9,4)} GHz')
    elif mixer == 'rr':
        plt.title(f'Readout Mixer {cal} Calibration at {round(fc*1e-9,4)} GHz')
    plt.show()
    
#%% misc_plots
def tof_plot(adc1,adc2):
    plt.figure()
    plt.title('time-of-flight calibration analysis')
    plt.plot(adc1)
    plt.plot(adc2)
    plt.legend(["adc1", "adc2"])
    plt.show()
    


def init_IQ_plot():
    '''initialize axes for continuous plotting on the IQ plane'''
    plot = sns.jointplot()
    plot.set_axis_labels('I [mV]', 'Q [mV]')
    plot.ax_marg_x.grid('off')
    plot.ax_marg_y.grid('off')
    plot.fig.tight_layout()
    ax = plt.gca()
    return plot, ax

def plot_single_shot(datadict, exp_pars={},qb_pars={},iteration=1, axes=0):

    datadict = {key: value*1e3 for key,value in datadict.items()} # convert to mV
    datadict = {key: value.tolist() for key,value in datadict.items()} # convert to list

    states = []
    # for key,value in datadict.items():
    #     print(key+':'+str(len(value))+'\n')
    [states.append(r'$|g\rangle$') for i in range(len(datadict['I']))]
    [states.append(r'$|e\rangle$') for i in range(len(datadict['Iexc']))]
    data = {
            'I [mV]':   np.hstack((datadict['I'],datadict['Iexc'])),
            'Q [mV]':   np.hstack((datadict['Q'],datadict['Qexc'])),
            'States':   states
                }
    dataF = pd.DataFrame(data=data)
    ax = sns.jointplot(data=dataF, x='I [mV]',y='Q [mV]',hue='States',space=0)
    txt=f'Single Shot {iteration}\n$P_r$ = -{qb_pars["rr_atten"]} dB | $T_\pi$ = {int(qb_pars["pi_len"]/2.4):d} ns | $A_\pi$ = {qb_pars["pi_amp"]:.3f}'
    # ax.ax_joint.text(s=txt, size='medium', horizontalalignment='center', color='black')
    ax.fig.suptitle(txt)
    ax.fig.subplots_adjust(top = 0.85)
    # plt.show()

def plot_coherence(t_data,v_b,wfms,exp_pars={},qb_pars={},wfm_pars={},calib_states={},savefig=False,project='',device_name='',qb=''):
    
    v_b = v_b[:,1:-1]
    qb_drive_freq = exp_pars['qubit_drive_freq']*1e-9
    fig, axs = plt.subplots(2,2,figsize=(12,9),dpi=150)
    # x_vector = x_vector*1e9
    # x_vector = data[0]
    t_data = t_data*1e6
    coherence = compute_coherence(v_b, calib_states)
    purity = compute_purity(v_b, calib_states)
    # wfms = data[4]
    plot_data = v_b
    # for i in range(v_b.shape[1]):
    #     plot_data[:,i] = compute_bloch(data[:,i],calib_pars=calib_states) 
    labels = ['$v_x$','$v_y$','$v_z$']

    print(coherence.shape)
    # coherence
    axs[0,0].plot(t_data, coherence, '-o', markersize = 3, c='C0')
    axs[0,0].set_ylabel('C(t)')
    axs[0,0].set_xlabel('Drive Duration ($\mu$s)')
    # purity
    axs[1,0].plot(t_data, purity, '-o', markersize = 3, c='C0')
    axs[1,0].set_ylabel(r'Tr[$\rho^2$(t)]')
    axs[1,0].set_xlabel('Drive Duration ($\mu$s)')
    axs[1,0].set_ylim([0,1.2])
    # bloch vector components
    axs[0,1].plot(t_data,plot_data[0,:],'-o',color='b',label=labels[0])
    axs[0,1].plot(t_data,plot_data[1,:],'-x',color='r',label=labels[1])
    axs[0,1].plot(t_data,plot_data[2,:],'-<',color='k',label=labels[2])
    axs[0,1].set_xlabel('Drive Duration ($\mu$s)')
    axs[0,1].legend()
    axs[0,1].set_ylim([-1.0,1.0])
    # # sx, sy waveforms
    axs[1,1].plot(wfms[0]*1e6, wfms[2], '-o', markersize = 1, c='r',label='$\sigma_y$',alpha=0.25)
    axs[1,1].plot(wfms[0]*1e6,wfms[1], '-o', markersize = 3, c='b',label='$\sigma_x$')
    
    axs[1,1].set_ylabel('A(t)')
    axs[1,1].set_xlabel('Drive Duration ($\mu$s)')
    axs[1,1].legend()
    # # ax.set_title(f'Rabi Measurement {iteration:03d}')
    textstr = f'Initial State: {exp_pars["initial-state"]}'+f'\n$\omega_d$ =  {qb_drive_freq:.4f} GHz\n'+f'$N$ = {exp_pars["n_avg"]}\n'+r'$\sigma$ = '+f'{wfm_pars["sigma"]*1e3:.1f} mV\n$T_1$ = {wfm_pars["T2"]*1e6:.1f}'+r'$\mu$s'
    # plt.tight_layout()
    plt.gcf().text(0.95, 0.15, textstr, fontsize=14)

    plt.tick_params(axis='both',direction='in',bottom=True, top=True, left=True, right=True,size=8)

    # if savefig:
        # plt.savefig(f'D:\\{project}\\{device_name}\\{qb}\\t-rabi-data\\fig_{iteration:03d}.png',dpi='figure')
        
#%% fit_funcs
def fit_res(f_data,z_data,res_type='notch'):
    fc = f_data[np.argmin(z_data)]
    if res_type == 'notch':
        # z_data = z_data-min(z_data)
        
        idx = np.argwhere(np.diff(np.sign(-(z_data - 0.5*max(z_data))))).flatten()
        try:
            fwhm = f_data[idx[1]] - f_data[idx[0]]
        except:
            print('FWHM not found')
            fwhm = 0
    return fc,fwhm

def fit_data(x_vector,y_vector,exp='t-rabi',dx=0.01,fitFunc='',verbose=0):

    '''
    fit experimental data
    sequence:          'Rabi','ramsey', 'T1' or 'T2'
    x_vector:           time data
    y_vector:           voltage data
    dx:                 sequence stepsize. Used for extracting the frequency of the data
    '''
    
    x_vector = x_vector*1e6
    y_vector = y_vector*1e3

    amp = (max(y_vector)-min(y_vector))/2
    offset = np.mean(y_vector)

    if exp == "t-rabi":
        fitFunction = rabi
        x_vector = x_vector*1e3
        period = 1e3/(extract_freq(x_vector, y_vector, dx,plot=0))
        print('Period Initial Guess: %.1f ns'%(period))
        phase = pi/2
        lb = [0.1*amp,0.1*period,0,-2*abs(offset)]
        ub = [10*amp,10*period,2*pi,2*abs(offset)]
        p0 = [amp,period,phase,offset]

    elif exp == "p-rabi":
        fitFunction = p_rabi
        x_vector = x_vector*1e-6
        period = 1/(extract_freq(x_vector, y_vector, dx,plot=0))
        print('Amplitude Initial Guess: %.3f'%(period))
        lb = [-2*amp,0.1*period,-2*abs(offset)]
        ub = [2*amp,10*period,2*abs(offset)]
        # p0 = [amp,period,phase,offset]
        p0 = [amp,period,offset]

    elif exp == "ramsey" or exp == 'tomography':
        f = extract_freq(x_vector, y_vector,dx,plot=0)
        print('Initial Guess for Freq:%.4f MHz'%(f))
        if x_vector[-1] > 20:
            tau = 30
        else:
            tau = 2
        phi = 0
        amp = abs(amp)
        # try:
        if fitFunc != 'envelope':
            p0 = [amp,f,phi,tau,offset]
            lb = [0.75*amp,0.1*f,-pi,0.01,-2*abs(offset)]
            ub = [2*amp,2*f,pi,500,2*abs(offset)]
            fitFunction = ramsey
            # fitted_pars, covar = scy.optimize.curve_fit(fitFunction, x_vector, y_vector,p0=p0,method='trf',bounds=[lb,ub],xtol=1e-12,maxfev=20e3)
        elif fitFunc == 'envelope':
            tau = 1
            env = get_envelope(y_vector, dx, distance=100)
            env = env(x_vector) + offset
            # env = get_envelope_LPF(x_vector, y_vector)*1e-3
            p0 = [amp,tau,offset]
            if offset < 0:
                p0 = [amp,tau,offset]
                lb = [0.95*amp,0.1,2*offset]
                ub = [2*amp,15,0.5*offset]
            elif offset >= 0:
                p0 = [amp,tau,offset]
                lb = [0.9*amp,0.1,0.9*offset]
                ub = [1.1*amp,15,1.1*offset]
            fitFunction = decay
            y_vector = env
    
    elif exp == 'z-gate':
        f = extract_freq(x_vector*1e-6, y_vector,dx,plot=0)
        print('Initial Guess for Freq:%.4f Hz'%(f))
        phi = pi
        amp = abs(amp)
        p0 = [amp,f,phi,offset]
        lb = [0.75*amp,0.1*f,-pi,-2*abs(offset)]
        ub = [2*amp,2*f,pi,2*abs(offset)]
        fitFunction = cos
        
    elif exp == "echo":
        if x_vector[-1] < 10:
            tau = 2
            tau_ub = 20
        else:
            tau = 20
            tau_ub = 300
        amp = y_vector[0] - y_vector[-1]
        p0 = [amp, tau, offset]
        amp_bounds = [0.95 * amp, 1.05 * amp]
        off_bounds = [0.95 * offset, 1.05 * offset]
        lb = [min(amp_bounds), 0.1, min(off_bounds)]
        ub = [max(amp_bounds), tau_ub, max(off_bounds)]
        # if offset < 0:
        #     lb = [0.95*amp,0.1,1.05*offset]
        #     ub = [1.05*amp,tau_ub,0.95*offset]
        # elif offset >= 0:
        #     lb = [0.95*amp,0.1,0.95*offset]
        #     ub = [1.05*amp,tau_ub,1.05*offset]
        fitFunction = decay
        # fitted_pars, covar = scy.optimize.curve_fit(decay, x_vector, y_vector,p0=p0,method='trf',bounds=[lb,ub],xtol=1e-12,maxfev=6000)
    elif exp == "T1":
        tau = 2
        amp = y_vector[0] - y_vector[-1]
        if amp < 0:
            p0 = [amp,tau,offset]
            lb = [10*amp,0.1,-2*abs(offset)]
            ub = [0.5*amp,300,2*abs(offset)]
        elif amp >= 0:
            p0 = [amp,tau,offset]
            lb = [0.5*amp,0.1,-2*abs(offset)]
            ub = [10*amp,300,2*abs(offset)]
        fitFunction = decay
        # fitted_pars, covar = scy.optimize.curve_fit(, x_vector, y_vector,p0=p0,method='trf',bounds=[lb,ub],xtol=1e-12,maxfev=6000)

    fitted_pars, covar = scy.optimize.curve_fit(fitFunction, x_vector, y_vector,p0=p0,method='trf',bounds=[lb,ub],xtol=1e-12,maxfev=40e3)
    error = np.sqrt(abs(np.diag(covar)))

    if verbose == 1:
        print('-'*100)
        print('Lower Bounds:',np.around(lb,1))
        print('Initial Guess:',np.around(p0,1))
        print('Upper Bounds:',np.around(ub,1))
        print('Best Fit Pars:',np.around(fitted_pars,1))
        print('Error:',np.around(error,1))
        print('-'*100)
    else:
        pass

    return fitted_pars,error

def rabi(x, amp,period,phase,offset):
    return amp*np.cos(2*pi*x/period+phase)+offset

def p_rabi(x, amp,period,offset):
    return amp*np.cos(2*pi*x/period)+offset

def ramsey(x,amp,f,phase,tau,offset):
    return amp*np.cos(2*pi*f*x+phase)*np.exp(-x/tau)+offset

def beats(x,amp,f1,f2,phase1,phase2,tau,offset):
    return amp*np.cos(pi*(f1+f2)*x+phase1)*np.cos(pi*(f2-f1)*x+phase2)*np.exp(-x/tau)+offset

def decay(x,amp,tau,offset):
    return amp*np.exp(-x/tau)+offset

def cos(x,amp,f,phase,offset):
    return amp*np.cos(2*pi*f*x+phase)+offset

def mod_cos(x,amp,B0,nu,phi1,phi2,tau,offset):
    return amp*np.cos(B0/nu*np.sin(2*np.pi*nu*x+phi1)+phi2)*np.exp(-x/tau)+offset
    # return amp*np.cos(np.cos(nu*x)*f*x)*np.exp(-x/tau)+offset

def mod_dec(x,amp1,f1,phi1,tau1,amp2,phi2,tau2,offset):
    return amp1*np.cos(2*np.pi*f1*x+phi1)*np.exp(-x/tau1)+ amp2*np.sin(2*np.pi*f1*x+phi2)*np.exp(-x/tau2)+offset
    # return amp*np.cos(np.cos(nu*x)*f*x)*np.exp(-x/tau)+offset

def extract_freq(t_vector,y_vector,dt,plot=0):
    N = len(t_vector)
    yf = scy.fft.fft(y_vector-np.mean(y_vector))
    xf = scy.fft.fftfreq(N,dt)[:round(N/2)]
    # print(len(xf))
    psd = 2.0/N * np.abs(yf[:round(N/2)])
    # print(len(psd))
    # print(psd)
    index_max = np.argmax(psd)
    if plot == 1:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(xf,psd)
        ax.set_xlabel('Frequency (MHz)')
        ax.set_ylabel('Power')
    # print(index_max)

    return xf[index_max]

def get_envelope(sig,dt, distance):
    # split signal into negative and positive parts
    sig = sig - np.mean(sig)
    u_x = np.where(sig > 0)[0]
    l_x = np.where(sig < 0)[0]
    u_y = sig.copy()
    u_y[l_x] = 0
    l_y = -sig.copy()
    l_y[u_x] = 0

    # find upper and lower peaks
    u_peaks, _ = scy.signal.find_peaks(u_y, distance=distance)
    l_peaks, _ = scy.signal.find_peaks(l_y, distance=distance)
    # use peaks and peak values to make envelope
    u_x = u_peaks
    u_y = sig[u_peaks]
    l_x = l_peaks
    l_y = sig[l_peaks]

    # add start and end of signal to allow proper indexing
    end = len(sig)
    u_x = np.concatenate(([0],u_x, [end]))*dt
    u_y = np.concatenate(([sig[0]],u_y, [sig[-1]]))
    l_x = np.concatenate(([0],l_x, [end]))*dt
    l_y = np.concatenate(([min(sig)],l_y, [np.mean(sig)]))
    # create envelope functions
    u = scy.interpolate.interp1d(u_x, u_y,kind='cubic',fill_value="extrapolate")
    # l = scipy.interpolate.interp1d(l_x, l_y,kind='cubic')
    return u

# def fit_2D_gauss(data_off,data_pi):
#     '''Takes single shot data and extracts mu and sigma for the 2 gaussian distributions'''
#     model = lmfit.models.gaussian2dModel()
#     params_off = model.guess
    
    
    
def get_envelope_LPF(x,sig):

    N = len(sig)
    Tmax = x[-1]
    cutoff = 100e6
    fs = N/Tmax

    env = butter_lowpass_filter(sig, cutoff, fs)
    plt.plot(x,env)
    plt.show()

    return env

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    sos = butter(order, normal_cutoff, btype='low', output = 'sos', analog=False)
    return sos

def butter_lowpass_filter(data,cutoff,fs):
    sos = butter_lowpass(cutoff, fs, order=5)
    y = scy.signal.sosfilt(sos, data)
    return y

def gen_WK_sig(fs,nu,tauk,Tmax):

    '''generate waveform using Wiener-Kinchin method'''

    N = int(fs*Tmax)+1
    dt = 1/fs
    df = 1/Tmax
    t = np.linspace(-Tmax/2,Tmax/2,N)

    # define autocorrelation and compute PSD
    autocorr = np.cos(2*np.pi*nu*t)*np.exp(-np.abs(t)/tauk)
    psd = 2/N*scy.fft.fft(autocorr)
    freqs = scy.fft.fftfreq(N,dt*1e6)[:round(N/2)]
    power = np.sqrt(np.abs(psd*df))
    # freq = np.fft.fftfreq(N,dt)[:round(2*nu/df)]
    # plt.plot(freq,psd[:round(2*nu/df)])
    # plt.show()

    # generate array of random phases
    phi_arr = np.random.rand(int(len(power)/2)) * 2*np.pi
    P_psd = np.zeros(len(phi_arr),dtype=complex)

    for i in range (1,len(phi_arr)):
        P_psd[i] = power[i]*np.exp(1j*phi_arr[i])

    # construct the 2 sided, symmetric fourier spectrum
    rand_psd = np.concatenate(([power[0]],np.conjugate(P_psd),np.flip(P_psd)))
    # inverse fourier transform to get timeseries
    signal = np.fft.ifft(rand_psd)
    signal = signal[int(len(signal)/2)+1:]

    if np.random.rand() < 0.5:
        signal = - signal
    else:
        pass

    return np.real(signal)/max(np.real(signal)),np.abs(psd[:round(N/2)]),freqs,autocorr[round(N/2)+1:]



