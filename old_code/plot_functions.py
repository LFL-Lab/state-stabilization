# -*- coding: utf-8 -*-
"""
Created on Mon May 24 13:13:39 2021

@author: lfl
"""
import matplotlib.pyplot as plt


import numpy as np
import pandas as pd
import scipy as sp
import scipy as scy
from matplotlib import cm
import sympy as sy
import csv
import itertools
from scipy.interpolate import interp1d
import scipy.fftpack
import time
import os
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.colors import LightSource
from types import SimpleNamespace
pi=np.pi
import zhinst.utils as ziut
import seaborn as sns; sns.set() # styling
from matplotlib.ticker import FormatStrFormatter
import imageio
sns.set_style('ticks')
plt.rcParams['figure.dpi'] = 300

def spec_plot(freq,I,Q,readout_power=-30,qubit_drive_amp=0.2):

    freq = freq*1e9
    I = I*1e3
    Q = Q*1e3
    mag = np.abs(I*I.conjugate()+Q*Q.conjugate())

    phase = np.unwrap(np.angle(I+1j*Q))
    sigma = np.std(I)
    peak = scy.signal.find_peaks(I,height=np.mean(I)+sigma,distance=100)
    # print(peak[0])
    peaks = peak[0]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(freq,phase,'-o', markersize = 3, c='C0')
    ax.set_xlabel('Frequency (GHz)')
    ax.set_ylabel('Phase (rad)')
    textstr = "$P_r = %.1f$ dBm\n Qubit Wfm Amp = %.1f mV" %(readout_power,qubit_drive_amp*1e3)
    plt.gcf().text(1, 0.25, textstr, fontsize=14)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(freq,Q,'-o', markersize = 3, c='C0')
    ax.set_xlabel('Frequency (GHz)')
    ax.set_ylabel('$I_{ON}-I_{OFF}$ (mV)')
    textstr = "$P_r = %.1f$ dBm\n Qubit Wfm Amp = %.1f mV" %(readout_power,qubit_drive_amp*1e3)
    plt.gcf().text(1, 0.25, textstr, fontsize=14)

    txt = "\n"
    for i in peaks:
        txt += "%.4f GHz \n" %(freq[i]*1e-9)
    print('Peaks are at: %s'%txt)

def fit_data(x_vector,y_vector,sequence='rabi',dt=0.01,qubitDriveFreq=3.8e9,amplitude_hd=1,mu=0,sigma=0,B0=0,nu=0,tauk=0,
                 pi2Width=50e-9,piWidth_Y=0,fitting = 1, plot=1,save_fig = 0,iteration=1,nAverages=1,sampling_rate=2.4e9,
                 integration_length=2e-6,AC_freq=5e9,cav_resp_time=5e-6,sweep=0,stepSize=5e-6,Tmax=5e-6,measPeriod=5e-6,nSteps=101,
                 prePulseLength=1500e-9,postPulseLength=1500e-9,threshold=10e-3,active_reset=False,rr_IF=5e6,source=0):

    '''
    fit experiment data

    sequence:          'Rabi', 'T1' or 'T2'
    complex_amplitude:  complex amplitude from the measurement
    x_vector:           x data
    fitting:     0:     do not fit
                  1:     do fit
    save_fig:    0:     do not save
                  1:     do save
    '''
    x_vector = x_vector*1e6
    abs_camplitude = y_vector*1e3

    if sequence == "rabi":
        amp = (max(abs_camplitude)-min(abs_camplitude))/2
        offset = np.mean(abs_camplitude)
        period = 1e3/(extract_freq(x_vector*1e3, abs_camplitude, dt,plot=0))
        # period = 110
        print('Period Initial Guess: %.1f ns'%(period))
        phase = 0
        p0 = [amp,period,phase,offset]
        fitted_pars, covar = scy.optimize.curve_fit(rabi, x_vector*1e3, abs_camplitude,p0=p0,xtol=1e-6,maxfev=3000)
        pi_pulse = np.round(1/2*fitted_pars[1])
        error = np.sqrt(abs(np.diag(covar)))
        print("Pi pulse duration is %.1f ns"%(pi_pulse))

        # return pi_pulse,error

    elif sequence == "ramsey":
        amp = abs_camplitude[0]-abs_camplitude[-1]
        offset = np.mean(abs_camplitude)
        f = extract_freq(x_vector, abs_camplitude,dt)
        print('Initial Guess for Freq:%.4f MHz'%(f))
        if x_vector[-1] > 10:
            tau = 10
        else:
            tau = 0.2
        phase = 0
        p0 = [amp,f,phase,tau,offset]
        try:
            fitted_pars, covar = scy.optimize.curve_fit(ramsey, x_vector, abs_camplitude,p0=p0,xtol=1e-6,maxfev=6000)
            detuning = fitted_pars[1]
            T_phi = fitted_pars[3]
            error = np.sqrt(abs(np.diag(covar)))
        except:
            print('fitting failed')
            fitted_pars = np.zeros(5)
            detuning = 0
            T_phi = 0
            error = 20*np.ones(5)

        # return detuning,T_phi,error

    elif sequence == "echo":
        amp = abs_camplitude[0]-abs_camplitude[-1]
        offset = np.mean(abs_camplitude)
        if x_vector[-1] < 2:
            tau = 0.6
        else:
            tau = 1.5
        p0 = [amp,tau,offset]
        try:
            fitted_pars, covar = scy.optimize.curve_fit(decay, x_vector, abs_camplitude,p0=p0,xtol=1e-6,maxfev=6000)
            T_2 = fitted_pars[1]
            error = np.sqrt(abs(np.diag(covar)))
            # if sum(error) > 2: # try fitting again, this time excluding the first few points
            #     fitted_pars, covar = scy.optimize.curve_fit(decay, x_vector[10:], abs_camplitude[10:],p0=p0,xtol=1e-6,maxfev=6000)
            #     T_2 = fitted_pars[1]
            #     error = np.sqrt(abs(np.diag(covar)))
        except:
            print('fitting failed')
            fitted_pars = np.zeros(3)
            T_2 = 0
            error = 20*np.ones(3)

        # return T_2,error

    elif sequence == "T1":
        amp = abs_camplitude[0]-abs_camplitude[-1]
        offset = np.mean(abs_camplitude)
        tau = 2
        p0 = [amp,tau,offset]
        try:
            fitted_pars, covar = scy.optimize.curve_fit(decay, x_vector, abs_camplitude,p0=p0,xtol=1e-6,maxfev=3000)
            T_1 = fitted_pars[1]
            error = np.sqrt(abs(np.diag(covar)))
        except:
            print('fitting failed')
            fitted_pars = np.zeros(3)
            T_1 = 0
            error = 20*np.ones(3)

    return fitted_pars,error




def plot_data(awg,x_vector,y_vector,sequence='rabi',dt=0.01,qubitDriveFreq=3.8e9,amplitude_hd=1,mu=0,sigma=0,B0=0,nu=0,tauk=0,
                 pi2Width=50e-9,piWidth_Y=0,fitting = 1, plot=1,save_fig = 0,iteration=1,nAverages=1,sampling_rate=2.4e9,
                 integration_length=2e-6,AC_freq=5e9,cav_resp_time=5e-6,sweep=0,stepSize=5e-6,Tmax=5e-6,measPeriod=5e-6,nSteps=101,
                 prePulseLength=1500e-9,postPulseLength=1500e-9,threshold=10e-3,active_reset=False,fitted_pars=np.zeros(5),pi_pulse=50,
                 fit_single_par_point=0,plot_mode=0,rr_IF=5e6,source=0):

    if fit_single_par_point == 0:
        x_vector = x_vector*1e6
        y_vector = y_vector*1e3
    elif fit_single_par_point == 1:
        x_vector[0] = x_vector[0]*1e6
        y_vector[0] = y_vector[0]*1e3
        x_vector[1] = x_vector[1]*1e6
        y_vector[1] = y_vector[1]*1e3

    if sequence == "rabi":
        fig, ax = plt.subplots()
        ax.plot(x_vector*1e3, y_vector, '-o', markersize = 3, c='C0')
        ax.set_ylabel('Digitizer Voltage (mV)')
        ax.set_xlabel('Pulse Duration (ns)')
        ax.plot(x_vector*1e3,rabi(x_vector*1e3, fitted_pars[0], fitted_pars[1], fitted_pars[2],fitted_pars[3]),'r')
        ax.set_title('Rabi Measurement %03d'%(iteration))
        textstr = '$\omega_d$ = %.4f GHz\n$A_d$ = %.2f V\n$T_{\pi/2}$ = %.1f ns\n$\mu$ = %d mV\n$\omega_{AC}$ = %.4f GHz\n$\hatn$ = %d'%(qubitDriveFreq*1e-9,amplitude_hd,round(fitted_pars[1]/4,1),mu*1e3,AC_freq*1e-9,nAverages)
        plt.gcf().text(0.95, 0.15, textstr, fontsize=14)

    elif sequence == "ramsey":
        if B0 != 0 or mu != 0 or sigma != 0: # if noise is injected, sz or sx, plot data and noise instances
            if fit_single_par_point == 1: # for plotting data averaged over many noise realizations post par sweep
                fig = plt.figure(figsize=((14,8)))
                ax1 = fig.add_subplot(111)
                # Plot data
                fontSize = 22
                tickSize = 20
                markersize= 10
                linewidth = 6
                # plot data
                ax1.plot(x_vector[0], y_vector[0], 'o', markersize = markersize, c='C0',label='$B_0 = 0$')
                if nu != 0:
                    ax1.plot(x_vector[1], y_vector[1], 'o', markersize = markersize, c='k', label='$B_0$ = %.1f mV | $\\nu = 2\pi\\times$%.1f kHz | $\\tau_k/\\tau_0$ = %.3f $\mu$s'%(B0*1e3,nu,tauk/fitted_pars[0][3]))
                elif nu == 0:
                    ax1.plot(x_vector[1], y_vector[1], 'o', markersize = markersize, c='k', label='$B_0 = %.1f mV$ | $\\tau_k/\\tau_0 = %.3f \mu s$'%(B0*1e3,tauk/fitted_pars[0][3]))
                ax1.set_ylabel('Digitizer Voltage (mV)',fontsize=fontSize)
                ax1.set_xlabel('Pulse Separation ($\mu$s)',fontsize=fontSize)
                for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
                    label.set_fontsize(tickSize)
                # plot fits
                ax1.plot(x_vector[0],ramsey(x_vector[0], fitted_pars[0][0], fitted_pars[0][1], fitted_pars[0][2],fitted_pars[0][3],fitted_pars[0][4]),'r',linewidth=linewidth,label='$\\tau_0 = %.1f \mu s$'%(fitted_pars[0][3]))
                if fitted_pars[1][1] < 0.01:
                    fitted_pars[1][1] = 0
                ax1.plot(x_vector[1],ramsey(x_vector[1], fitted_pars[1][0], fitted_pars[1][1], fitted_pars[1][2],fitted_pars[1][3],fitted_pars[1][4]),'g',linewidth=linewidth,label='$\omega = 2\pi\\times%.1f$ MHz, $\\tau/\\tau_0$ = %.1f $\mu$s'%(fitted_pars[1][1],fitted_pars[1][3]/fitted_pars[0][3]))
                ax1.legend(loc='upper right',prop={'size':22})
                # ax1.set_xlim([0,x_vector[0][-1]])
                textstr = '$T_{\pi/2}$=%.1f ns\n$\omega_d$ = %.5f GHz\n$A_d$ = %.3f V\n$\Delta$=%.3f MHz\n$\\tau$/$\\tau_0$=%.1f\n$\mu$ = %d mV\n$\omega_{AC}$ = %.4f GHz\n$\sigma$ = %.1f mV\n$\hatn$ = %d'%(pi2Width*1e9,qubitDriveFreq*1e-9,amplitude_hd,fitted_pars[0][1],fitted_pars[1][3]/fitted_pars[0][3],mu*1e3,AC_freq*1e-9,sigma*1e3,nAverages)
                # plt.gcf().text(0.95, 0.15, textstr, fontsize=fontSize)
            elif fit_single_par_point == 0: #for plotting non-par sweep data and par sweep data during the sweep
                fig = plt.figure(figsize=(20,18))
                ax1 = fig.add_subplot(2,2,(1,2))
                # Plot data
                fontSize = 36
                tickSize = 24
                markersize= 12
                linewidth = 6
                ax1.plot(x_vector, y_vector, 'o', markersize = markersize, c='C0')
                ax1.set_ylabel('Digitizer Voltage (mV)',fontsize=fontSize)
                for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
                    label.set_fontsize(tickSize)
                # Retrieve and plot telegraph noise instrance
                waveforms = ziut.parse_awg_waveform(awg.get('/dev8233/awgs/0/waveform/waves/0')['dev8233']['awgs']['0']['waveform']['waves']['0'][0]['vector'],channels=2)
                t_arr = np.linspace(0,x_vector[-1],len(waveforms[0]))
                ax2 = fig.add_subplot(2,2,3)
                ax2.plot(t_arr,1e3*waveforms[0],linewidth=linewidth,color='b')
                ax3 = fig.add_subplot(2,2,4)
                ax3.plot(t_arr,1e3*waveforms[1],linewidth=linewidth,color='r')
                for label in (ax2.get_xticklabels() + ax2.get_yticklabels()):
                    label.set_fontsize(tickSize)
                for label in (ax3.get_xticklabels() + ax3.get_yticklabels()):
                    label.set_fontsize(tickSize)
                ax2.set_xlabel('Pulse Separation ($\mu$s)',fontsize=fontSize)
                ax2.set_ylabel('Amplitude (mV)',fontsize=fontSize)
                ax2.set_title('$\sigma_x$ waveform',fontsize=fontSize)
                ax3.set_xlabel('Pulse Separation ($\mu$s)',fontsize=fontSize)
                ax3.set_title('$\sigma_z$ waveform',fontsize=fontSize)
                if tauk != 0:
                    textstr = '$T_{\pi/2}$=%.1f ns\n$\omega_d$ = %.6f GHz\n$A_d$ = %.3f V\n$\Delta$=%.3f MHz\n$T_2^*$=%.2f $\mu$s\n$\mu$ = %d mV\n$\omega_{AC}$ = %.4f GHz\n$\sigma$ = %.1f mV\n$B_0$ = %.2f mV\n$\\tau_k$ = %.3f $\mu s$\n$\\nu$ = %.1f kHz\n$\hatn$ = %d'%(pi2Width*1e9,qubitDriveFreq*1e-9,amplitude_hd,fitted_pars[1],fitted_pars[3],mu*1e3,AC_freq*1e-9,sigma*1e3,B0*1e3,tauk,nu,nAverages)
                elif nu == 0:
                    textstr = '$T_{\pi/2}$=%.1f ns\n$\omega_d$ = %.6f GHz\n$A_d$ = %.3f V\n$\Delta$=%.3f MHz\n$T_2^*$=%.2f $\mu$s\n$\mu$ = %d mV\n$\omega_{AC}$ = %.4f GHz\n$\sigma$ = %.1f mV\n$B_0$ = %.2f mV\n$\\tau_k$ = %.3f $\mu s$\n$\hatn$ = %d'%(pi2Width*1e9,qubitDriveFreq*1e-9,amplitude_hd,fitted_pars[1],fitted_pars[3],mu*1e3,AC_freq*1e-9,sigma*1e3,B0*1e3,nu,nAverages)
                if plot_mode == 0:
                    # for plotting single realization
                    # ax1.plot(x_vector,ramsey(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2],fitted_pars[3],fitted_pars[4]),'r',linewidth=linewidth)
                    plt.gcf().text(0.925, 0.25, textstr, fontsize=fontSize,bbox=dict(boxstyle='round,rounding_size=1.25',facecolor='silver',alpha=0.5))
                    ax1.set_title('Ramsey Measurement %03d' %(iteration),size=fontSize)
                elif plot_mode == 1:
                    textstr = '$T_2^*$=%.2f $\mu$s\n$\sigma$ = %.1f mV'%(fitted_pars[3],sigma*1e3)
                    plt.gcf().text(0.925, 0.45, textstr, fontsize=fontSize+10,bbox=dict(boxstyle='round,rounding_size=1.25',facecolor='silver',alpha=0.5))
                elif plot_mode == 2:
                    textstr = '$\omega$=%.3f MHz\n$T_2^*$=%.2f $\mu$s\n$\sigma$ = %.1f mV\n$B_0$ = %.2f mV\n$\\tau_k$ = %.3f $\mu s$'%(fitted_pars[1],fitted_pars[3],sigma*1e3,B0*1e3,nu)
                    plt.gcf().text(0.925, 0.45, textstr, fontsize=fontSize+10,bbox=dict(boxstyle='round,rounding_size=1.25',facecolor='silver',alpha=0.5))
        elif mu == 0 and B0 == 0 and sigma == 0: # plot data without any noise instances
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            fontSize = 16
            tickSize = 11
            markersize= 6
            linewidth = 2
            ax1.plot(x_vector, y_vector, 'o', markersize = markersize, c='C0')
            ax1.set_ylabel('Digitizer Voltage (mV)',fontsize=fontSize)
            ax1.set_xlabel('Pulse Separation ($\mu$s)')
            for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
                label.set_fontsize(tickSize)
            ax1.plot(x_vector,ramsey(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2],fitted_pars[3],fitted_pars[4]),'r',linewidth=linewidth)
            textstr = '$T_{\pi/2}$=%.1f ns\n$\omega_d$ = %.6f GHz\n$A_d$ = %.3f V\n$\Delta$=%.3f MHz\n$T_2^*$=%.2f $\mu$s\n$\mu$ = %d mV\n$\omega_{AC}$ = %.4f GHz\n$\sigma$ = %.1f mV\n$B_0$ = %.1f mV\n$\\tau_k$ = %.3f $\mu s$\n$\hatn$ = %d'%(pi2Width*1e9,qubitDriveFreq*1e-9,amplitude_hd,fitted_pars[1],fitted_pars[3],mu*1e3,AC_freq*1e-9,sigma*1e3,B0*1e3,nu,nAverages)
            plt.gcf().text(0.95, 0.15, textstr, fontsize=fontSize)
            ax1.set_title('Ramsey Measurement %03d' %(iteration),size=fontSize)

    elif sequence == "echo":
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(x_vector, y_vector, '-o', markersize = 3, c='C0')
        ax.set_ylabel('Digitizer Voltage (mV)')
        ax.set_xlabel('Pulse Separation ($\mu$s)')
        ax.plot(x_vector,decay(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2]),'r')
        textstr = '$T_{\pi/2}$=%.1f ns\n$T_{\pi(Y)} = %.1f$ ns\n$\omega_d$ = %.4f GHz\n$A_d$ = %.2f V\n$T_2$=%.2f $\mu$s\n$\mu$ = %.3f V\n$\omega_{AC}$ = %.4f GHz\n$\sigma$ = %.3f V\n$B_0$ = %.3f V\n$\\tau_k$ = %.2f $\mu s$\n$\hatn$ = %d'%(pi2Width*1e9,piWidth_Y*1e9,qubitDriveFreq*1e-9,amplitude_hd,fitted_pars[1],mu,AC_freq*1e-9,sigma,B0,nu,nAverages)
        ax.set_title('Echo Measurement %03d' %(iteration))
        plt.gcf().text(0.95, 0.15, textstr, fontsize=14)

    elif sequence == "T1":
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(x_vector, y_vector, '-o', markersize = 3, c='C0')
        ax.set_ylabel('Digitizer Voltage (mV)')
        ax.set_xlabel('Delay ($\mu$s)')
        ax.plot(x_vector,decay(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2]),'r')
        textstr = '$T_{\pi/2}$=%.1f ns\n$\omega_d$ = %.4f GHz\n$A_d$ = %.2f V\n$T_1$=%.2f $\mu$s\n$\mu$ = %.3f V\n$\omega_{AC}$ = %.4f GHz\n$\sigma$ = %.3f V\n$B_0$ = %.2f V\n$\\tau_k$ = %.2f $\mu s$\n$\hatn$ = %d'%(pi2Width*1e9,qubitDriveFreq*1e-9,amplitude_hd,fitted_pars[1],mu,AC_freq*1e-9,sigma,B0,nu,nAverages)
        ax.set_title('T1 Measurement %03d' %(iteration))
        plt.gcf().text(0.95, 0.15, textstr, fontsize=14)

    return fig

def sweep_plot(x_data,y_data,z_data,data='tau'):
    fig,ax = plt.subplots(subplot_kw={"projection": "3d"},dpi=300)
    # fig = plt.figure(figsize=(3,3),constrained_layout=True)
    # ax = fig.add_subplot(111,projection="3d")
    if len(x_data) < len(y_data):
        y_data_mod = y_data
        x_data_mod = np.append(np.zeros(np.abs(len(y_data)-len(x_data))),x_data,axis=0)
        z_data_mod = np.append(np.ones((np.abs(len(y_data)-len(x_data)),len(y_data))),z_data,axis=0)
    elif len(x_data) > len(y_data):
        x_data_mod = x_data
        y_data_mod = np.append(np.zeros(np.abs(len(y_data)-len(x_data))),y_data,axis=0)
        z_data_mod = np.append(np.ones((len(x_data),np.abs(len(y_data)-len(x_data)))),z_data,axis=1)
    else:
        x_data_mod = x_data
        y_data_mod = y_data
        z_data_mod = z_data

    X,Y = np.meshgrid(x_data_mod,y_data_mod)
    surf = ax.plot_surface(X,Y,z_data_mod,cmap=plt.cm.jet,linewidth=0,shade=False)

    labelsize = 16
    tickLabelSize = 6
    labelpad = 2
    # ax.set_zlim(0,2)

    ax.set_xlabel("$\\nu$ (MHz)",labelpad=0)
    # ax.set_xlabel("$B_0 (mV)$",labelpad=0)
    ax.set_ylabel("$\\tau_k (\mu s)$",labelpad=0)
    if data == 'tau':
        label = "$\\tau$ ($\mu$s)"
    elif data == 'omega':
        label = "$\omega$ (MHz)"
    ax.set_zlabel(label,labelpad=-12)
    ax.view_init(elev=30,azim=-120)

    plt.xticks(rotation=0)
    plt.yticks(rotation=60)
    ax.zaxis.set_ticklabels([])
    ax.set_xlim((min(x_data),max(x_data)))
    # ax.set_xlim((min(x_data*1e3),max(x_data*1e3)))
    ax.set_ylim((0.1,10))
    ax.set_zlim((1,8))
    cbar = fig.colorbar(surf,cmap=cm.coolwarm, ticks=np.around(np.linspace(0,ax.get_zlim()[1],12)),pad=0.1,values=np.linspace(0,ax.get_zlim()[1],1000),shrink=0.5,location='top')
    cbar.set_label(label='($\mu$s)',y=-5,labelpad=-45)
    ax.tick_params(axis='both',pad=-2)
    plt.show()

def plot_slice(par1,par2,par3,ydata,data='tau'):
    '''
    Plots a slice of a 3d plot

    Parameters
    ----------
    par1 : double
        B0
    par2 : double
        nu
    par3: double
        tau_k
    ydata : double (1D array)
    data_type : string, optional
       Defines whether data are omega or tau. The default is 'tau'.


    '''
    fontsize = 22
    tickSize = 20
    fig = plt.figure()
    ax = fig.add_subplot(111)

    if len(par1) > 1 and len(par2) == len(par3):
        ax.plot(par1*1e3,ydata,label='$\\tau_k$ = %.1f $\mu s$'%(par3),c='r',linewidth=4)
        ax.set_xlabel('B_0 (mV)',fontsize=fontsize)
    elif len(par2) > 1 and len(par1) == len(par3):
        ax.plot(par2,ydata,label='$B_0$ = %.2f mV\n$\\tau_k$ = %.1f $\mu s$'%(round(par1*1e3),par3),linewidth=2,c='b')
        ax.set_xlabel('$\\nu (kHz)',fontsize=fontsize)
    elif len(par3) > 1 and len(par1) == len(par2):
        ax.plot(par3,ydata,label='$B_0$ = %.2f mV, $\\nu = %d kHz'%(round(par1*1e3),par2),linewidth=2,c='b')
        ax.set_xlabel('$\\tau_k$ ($\mu$s)',fontsize=fontsize)
    else:
        raise Exception('Incorrect parameter combination!')

    ax.set_ylabel('$\\tau/\\tau_0$',fontsize=fontsize)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.legend(loc='center right',fontsize=fontsize-6)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(tickSize)
    plt.show()

def fit_sweep(par1,par2,par3,sweep='sweep_001'):

    T2_b_arr = np.zeros((len(par1),len(par2)))
    detun_b_arr = np.zeros((len(par1),len(par2)))
    T2_arr = np.zeros((len(par1),len(par2)))
    detun_arr = np.zeros((len(par1),len(par2)))

    for i in range(len(par2)):
        for j in range(len(par1)):
            start = time.time()
            # t_b, data_b , t , data, par_dictionary = extract_data(sweep,par1=par1[j], par2 = par2[i])
            # t_b, data_b , t , data, par_dictionary = extract_data(sweep,par1=par3[0], par2 = par2[i],par3=par1[j]*1e3)
            print('Now fitting B0 = %.4f | nu = %.3f | tau = %.3f' %(par1[j],par3[0],par2[i]))
            # fitted_pars, error = fit_data(x_vector=t, y_vector=np.mean(data,axis=0),dt=t[-1]/data.shape[1], **par_dictionary)
            fig,fitted_pars = plot_single_par_point(par1=par1[j], par2=par2[i], par3=par3[0],sweep=sweep,plot=1)
            plt.savefig(os.path.join('E:\\generalized-markovian-noise\\CandleQubit_6\\sweep_data\\ramsey\\'+sweep+'\\plot_images\\plot_B0_%d_uV_nu_%d_kHz_tau_%d_ns.png' %(round(par1[j]*1e6),round(par3[0]*1e3),round(par2[i]*1e3))) , bbox_inches='tight')
            plt.close(fig)
            T2_arr[j,i] = fitted_pars[1][3]
            detun_arr[j,i] = fitted_pars[1][1]
            # fit_beats(sequence='ramsey', dt=Tmax/data.shape[1], plot=plot,x_vector=t, y_vector=-np.mean(data,axis=0),AC_pars=[0.6,0.08],RT_pars=[par1[j],par2[i]])
            # fitted_pars, error = fit_data(x_vector=t, y_vector=np.mean(data_b,axis=0),dt=t_b[-1]/data_b.shape[1], **par_dictionary)
            T2_b_arr[j,i] = fitted_pars[0][3]
            detun_b_arr[j,i] = fitted_pars[0][1]
            end = time.time()
            print('fitting one point took:%.2f sec'%(end-start))

    return detun_arr, T2_arr, detun_b_arr, T2_b_arr

# def fit_single_instance(par1,par2,sweep='sweep_001',plot=1):
#     # start = time.time()
#     t_b, data_b , t , data = extract_data(sweep,par1=par1, par2 = par2)
#     # data = data[instance,:]
#     # detun , T2, error = pulse_plot1d(sequence='ramsey', dt=20/data.shape[1], plot=plot,x_vector=t, y_vector=np.mean(data,axis=0),AC_pars=[0.6,0.08],RT_pars=[par1,par2])
#     fit_beats(sequence='ramsey', dt=5/data.shape[1], plot=plot,x_vector=t, y_vector=-np.mean(data,axis=0),AC_pars=[0.6,0.08],RT_pars=[par1,par2])
#     end = time.time()

    # return detun, T2
def extract_data(sweep,B0,nu,tauk,meas_device='CandleQubit_6',nMeasurements=128,nBackMeasurements=128):

    filename = 'B0_%d_uV_nu_%d_Hz_tau_%d_ns' %(round(B0*1e6),round(nu*1e3),round(tauk*1e3))
    datafile = "E:\\generalized-markovian-noise\\%s\\sweep_data\\ramsey\\%s\\data\\data_%s.csv"%(meas_device,sweep,filename)
    # get data
    tdata_background = pd.read_csv(datafile,on_bad_lines='skip',skiprows=3,header=None,nrows=1).to_numpy(np.float64)[0]
    tdata = pd.read_csv(datafile,on_bad_lines='skip',skiprows=5,header=None,nrows=1).to_numpy(np.float64)[0]
    ydata_background = pd.read_csv(datafile,on_bad_lines='skip',skiprows=7,header=None,nrows=nBackMeasurements).to_numpy(np.float64)
    ydata = pd.read_csv(datafile,on_bad_lines='skip',skiprows=7+nBackMeasurements+1,header=None,nrows=nMeasurements).to_numpy(np.float64)

    # get exp parameters
    pars = pd.read_csv(datafile,on_bad_lines='skip',header=[0],nrows=1)
    keys = pars.keys()
    values = pars.values
    dictionary = dict(zip(keys,values[0]))

    return tdata_background, ydata_background, tdata, ydata, dictionary

def plot_single_par_point(par1,par2,par3,sweep,meas_device='CandleQubit_6',plot=1):
    '''
    DESCRIPTION: Plots the averaged ramsey trace over the different noise realizations
    '''
    tdata_background,ydata_background, tdata, ydata, exp_pars = extract_data(sweep=sweep, B0=par1, nu=par2,tauk=par3,meas_device=meas_device)
    # average
    ydata_avg_background = np.mean(ydata_background,axis=0)
    ydata_avg = np.mean(ydata,axis=0)
    #fit
    fitted_pars,error = fit_data(x_vector=tdata, y_vector=ydata_avg,sequence='ramsey',dt=tdata[-1]/len(tdata),**exp_pars)
    fitted_pars_background,error_b = fit_data(x_vector=tdata_background, y_vector=ydata_avg_background,sequence='ramsey',dt=tdata_background[-1]/len(tdata_background),**exp_pars)
    #plot
    if plot == 1:
        fig = plot_data(awg=None,x_vector=[tdata_background,tdata],y_vector=[ydata_avg_background,ydata_avg],sequence='ramsey',fitted_pars=[fitted_pars_background,fitted_pars],**exp_pars,fit_single_par_point=1)

    return fig,[fitted_pars_background,fitted_pars]

def rabi(x, amp,period,phase,offset):
    return amp*np.cos(2*pi*x/period+phase)+offset

def ramsey(x,amp,f,phase,tau,offset):
    return amp*np.cos(2*pi*f*x+phase)*np.exp(-x/tau)+offset

def beats(x,amp,f1,f2,phase1,phase2,tau,offset):
    return amp*np.cos(pi*(f1+f2)*x+phase1)*np.cos(pi*(f2-f1)*x+phase2)*np.exp(-x/tau)+offset

def decay(x,amp,tau,offset):
    return amp*np.exp(-x/tau)+offset

def extract_freq(t_vector,y_vector,dt,plot=0):
    N = len(t_vector)
    dt = dt*1e6
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

def Volt2dBm(data):

    return 10*np.log10(1e3*data**2/50)

def Watt2dBm(x):
    '''
    converts from units of Watts to dBm
    '''
    return 10.*np.log10(x*1000.)

def plot_single_shot(data_OFF,data_pi):
    states = []
    [states.append('|g>') for i in range(len(data_OFF.real))]
    [states.append('|e>') for i in range(len(data_pi.real))]
    data = {
            'I (mV)':   np.hstack((data_OFF.real*1e3,data_pi.real*1e3)),
            'Q (mV)':   np.hstack((data_OFF.imag*1e3,data_pi.imag*1e3)),
            'states':   states
                }
    dataF = pd.DataFrame(data=data,dtype=np.float32)
    plot = sns.jointplot(data=dataF, x='I (mV)',y='Q (mV)',hue='states')
    # plot.fig.set_figwidth(7)
    # plot.fig.set_figheight(4)

def plot_noise(length=1000,gif_make=0,wfm1=[],wfm2=[]):
    mu = 325
    A_d = 200
    B0 = 100
    fontsize = 40
    linewidth = 10
    ticksize = 30
    if gif_make == 0:
        AC_stark_waveform = np.concatenate((0*np.ones(200),mu*np.ones(500),mu*np.ones(100),mu*np.ones(length)+np.random.normal(loc=0.325, scale=30, size=length),mu*np.ones(100),mu*np.ones(250),0*np.ones(200)))
        qubit_channel_waveform = np.concatenate((0*np.ones(200),A_d*np.zeros(500),A_d*np.ones(100),B0*gen_tel_noise(length, tau=0.75e6, dt=0.01),A_d*np.ones(100),0*np.ones(250),0*np.ones(200)))
    elif gif_make == 1:
        AC_stark_waveform = np.concatenate((wfm1[:800],wfm1[800:length+800],wfm1[-550:]))
        print(len(AC_stark_waveform))
        qubit_channel_waveform = np.concatenate((wfm2[:800],wfm2[800:length+800],wfm2[-550:]))
    t = np.linspace(0,1e-6,len(qubit_channel_waveform))
    fig = plt.figure(figsize=(20,16))
    ax = fig.add_subplot(111)
    ax.plot(t,AC_stark_waveform,c='r',linewidth=linewidth,label='AC Stark Channel Waveform ($\sigma_z$)')
    ax.set_xlabel('time ($\mu$s)',fontsize=fontsize)
    ax.plot(t,qubit_channel_waveform,linewidth=linewidth,label='Qubit Channel Waveform ($\sigma_x$)')
    ax.set_ylabel('Voltage (mV)',fontsize=fontsize)
    ax.set_ylim([-250,600])
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(ticksize)
    ax.legend(loc='upper center',prop={'size':30})
    plt.show()
    return fig,AC_stark_waveform,qubit_channel_waveform

def create_wfm_gif(Lmax=1000):
    numImag = 50
    filenames = []
    fig,AC_stark_waveform,qubit_channel_waveform = plot_noise(length=Lmax,gif_make=0)
    path = 'G:\Shared drives\LFL\Projects\Generalized Markovian noise\MarchMeeting2022\gif\\'
    for i in range(numImag):
        fig,wfm1,wfm2 = plot_noise(int(i*Lmax/numImag),gif_make=1,wfm1=AC_stark_waveform,wfm2=qubit_channel_waveform)
        filename = f'{i}.png'
        filenames.append(filename)
        plt.savefig('G:\Shared drives\LFL\Projects\Generalized Markovian noise\MarchMeeting2022\gif\\'+filename,  bbox_inches='tight')
        plt.close(fig)

   # build gif
    with imageio.get_writer(os.path.join(path,'ramsey_sequence.gif'), mode='I') as writer:
        for filename in filenames:
            image = imageio.imread(os.path.join(path,filename))
            writer.append_data(image)

    # Remove files
    for filename in set(filenames):
        os.remove(filename)

def gen_tel_noise(numPoints,tau,dt):

    signal = np.ones(numPoints)*(-1)**np.random.randint(0,2)
    for i in range(1,numPoints-1):
        if np.random.rand() < 1/(2*tau*1e-6/dt)*np.exp(-1/(2*tau*1e-6/dt)):
            signal[i+1] = - signal[i]
        else:
            signal[i+1] = signal[i]
    return signal

# def extract_data(sweep,par1,par2,par3,meas_device='CandleQubit_6'):
#     # filename = 'B0_%d_uV_nu_%d_kHz_tau_%d_ns' %(round(par1*1e6),round(par3),round(par2*1e3))
#     filename = 'B0_%d_uV_nu_%d_kHz_tau_%d_ns' %(round(par1*1e6),round(par2*1e3),round(par3*1e3))
#     # filename = 'RTN_B0_%d_uV_tau_%d_ns' %(round(par1*1e6),round(par2*1e3))
#     with open("E:\\generalized-markovian-noise\\%s\\sweep_data\\ramsey\\%s\\data\\data_%s.csv"%(meas_device,sweep,filename)) as datafile:
#         csv_reader = csv.reader(datafile,delimiter=',')
#         file_data = list(csv_reader)

#         for row in file_data:
#             if row[0] == "Background Time Data":
#                 tdata_background = np.array(file_data[file_data.index(row)+1],dtype=float)
#             if row[0] == "Time Data":
#                 datafile.seek(0)
#                 tdata= np.array(file_data[file_data.index(row)+1],dtype=float)
#             if row[0] == "Background Data: Channel 1":
#                 line_start_background = file_data.index(row) + 1
#             if row[0] == "Background Data: Channel 2":
#                 line_end_background = file_data.index(row) - 1
#             if row[0] == "Data: Channel 1":
#                 line_start = file_data.index(row) + 1
#             if row[0] == "Data: Channel 2":
#                 line_end = file_data.index(row) - 1


#         datafile.seek(0)
#         # extract traces
#         ydata_background = np.zeros((line_end_background-line_start_background+1,len(tdata_background)))
#         line = line_start_background
#         while line <= line_end_background:
#             datafile.seek(0)
#             # print(np.array(next(itertools.islice(csv_reader, line,None)),dtype=np.float32))
#             ydata_background[line-line_start_background,:] = np.array(next(itertools.islice(csv_reader, line,None)),dtype=np.float32)
#             line += 1
#             # print(ydata_background)

#         datafile.seek(0)
#         trace = np.array(next(itertools.islice(csv_reader, line_start,None)),dtype=np.float32)
#         # ydata = np.zeros((line_end-line_start+1,len(trace)))
#         ydata = np.zeros((line_end-line_start+1,len(tdata)))
#         line = line_start
#         while line <= line_end:
#             datafile.seek(0)
#             # data = np.array(next(itertools.islice(csv_reader, line,None)),dtype=np.float32)
#             # ydata[line-line_start,:] = data[:148]
#             ydata[line-line_start,:] = np.array(next(itertools.islice(csv_reader, line,None)),dtype=np.float32)[:]
#             line += 1

#         #extract data point parameters (qubitDriveFreq, ACstarkFreq, etc.)
#         datafile.seek(0)
#         keys = file_data[0]
#         values =  file_data[1]
#         dictionary = dict.fromkeys(keys,0)

#         i = 0
#         for key in dictionary:
#             if key == 'AC_pars' or key == 'RT_pars':
#                 pars = np.zeros((3,1))
#                 values[i] = values[i].replace('[','')
#                 values[i] = values[i].replace(']','')
#                 values[i] = values[i].replace(',','')
#                 pars[0] = float(values[i][1:values[i].find(' ')])
#                 pars[1] = float(values[i][values[i].find(' ')+1:values[i].find(' ',values[i].find(' ')+1)])
#                 pars[2] = float(values[i][values[i].find(' ',values[i].find(' ')+1):])
#                 dictionary[key] = pars
#             try:
#                 dictionary[key] = float(values[i])
#             except:
#                 dictionary[key] = values[i]
#             i += 1


#     return tdata_background, ydata_background, tdata, ydata, dictionary
# def plot_data(awg,x_vector,y_vector,sequence='rabi',dt=0.01,qubitDriveFreq=3.8e9,amplitude_hd=1,AC_pars=[0,0],RT_pars=[0,0],
#                  pi2Width=50e-9,piWidth_Y=0,fitting = 1, plot=1,save_fig = 0,iteration=1,nAverages=1,sampling_rate=2.4e9,
#                  integration_length=2e-6,AC_freq=5e9,cav_resp_time=5e-6,sweep=0,stepSize=5e-6,Tmax=5e-6,measPeriod=5e-6,nSteps=101,
#                  prePulseLength=1500e-9,postPulseLength=1500e-9,threshold=10e-3,active_reset=False,fitted_pars=np.zeros(5),pi_pulse=50,
#                  fit_single_par_point=0,plot_mode=0,rr_IF=5e6,source=0):

#     if fit_single_par_point == 0:
#         x_vector = x_vector*1e6
#         y_vector = y_vector*1e3
#     elif fit_single_par_point == 1:
#         x_vector[0] = x_vector[0]*1e6
#         y_vector[0] = y_vector[0]*1e3
#         x_vector[1] = x_vector[1]*1e6
#         y_vector[1] = y_vector[1]*1e3

#     if sequence == "rabi":
#         fig, ax = plt.subplots()
#         ax.plot(x_vector*1e3, y_vector, '-o', markersize = 3, c='C0')
#         ax.set_ylabel('Digitizer Voltage (mV)')
#         ax.set_xlabel('Pulse Duration (ns)')
#         ax.plot(x_vector*1e3,rabi(x_vector*1e3, fitted_pars[0], fitted_pars[1], fitted_pars[2],fitted_pars[3]),'r')
#         ax.set_title('Rabi Measurement %03d'%(iteration))
#         textstr = '$\omega_d$ = %.4f GHz\n$A_d$ = %.2f V\n$T_{\pi/2}$ = %.1f ns\n$\mu$ = %d mV\n$\omega_{AC}$ = %.4f GHz\n$\hatn$ = %d'%(qubitDriveFreq*1e-9,amplitude_hd,round(fitted_pars[1]/4,1),mu*1e3,AC_freq*1e-9,nAverages)
#         plt.gcf().text(0.95, 0.15, textstr, fontsize=14)

#     elif sequence == "ramsey":
#         if RT_pars[0] != 0 or AC_pars[0] != 0 or AC_pars[1] != 0: # if noise is injected, sz or sx, plot data and noise instances
#             if sweep == 1 and fit_single_par_point == 1: # for plotting data averaged over many noise realizations post par sweep
#                 fig = plt.figure(figsize=((14,8)))
#                 ax1 = fig.add_subplot(111)
#                 # Plot data
#                 fontSize = 22
#                 tickSize = 20
#                 markersize= 10
#                 linewidth = 6
#                 # plot data
#                 ax1.plot(x_vector[0], y_vector[0], 'o', markersize = markersize, c='C0',label='$B_0 = 0$')
#                 nu = float(RT_pars[RT_pars.find(' ',RT_pars.find(' ')+1):])
#                 if nu != 0:
#                     ax1.plot(x_vector[1], y_vector[1], 'o', markersize = markersize, c='k', label='$B_0$ = %.1f mV | $\\nu = 2\pi\\times$%.1f MHz | $\\tau_k/\\tau_0$ = %.3f $\mu$s'%(float(RT_pars[1:RT_pars.find(' ')])*1e3,float(RT_pars[RT_pars.find(' ',RT_pars.find(' ')+1):]),float(RT_pars[RT_pars.find(' ')+1:RT_pars.find(' ',RT_pars.find(' ')+1)])/fitted_pars[0][3]))
#                 elif nu == 0:
#                     ax1.plot(x_vector[1], y_vector[1], 'o', markersize = markersize, c='k', label='$B_0 = %.1f mV$ | $\\tau_k/\\tau_0 = %.3f \mu s$'%(float(RT_pars[1:RT_pars.find(' ')])*1e3,float(RT_pars[RT_pars.find(' ')+1:RT_pars.find(' ',RT_pars.find(' ')+1)])/fitted_pars[0][3]))
#                 # ax1.plot(x_vector[1], y_vector[1], 'o', markersize = markersize, c='k', label='$B_0 = %.1f mV$ | $\\tau_k/\\tau_0$ = %.3f $\mu$s'%(float(RT_pars[1:RT_pars.find(' ')])*1e3,float(RT_pars[RT_pars.find(' ')+1:])/fitted_pars[0][3]))
#                 ax1.set_ylabel('Digitizer Voltage (mV)',fontsize=fontSize)
#                 ax1.set_xlabel('Pulse Separation ($\mu$s)',fontsize=fontSize)
#                 for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
#                     label.set_fontsize(tickSize)
#                 # plot fits
#                 ax1.plot(x_vector[0],ramsey(x_vector[0], fitted_pars[0][0], fitted_pars[0][1], fitted_pars[0][2],fitted_pars[0][3],fitted_pars[0][4]),'r',linewidth=linewidth,label='$\\tau_0 = %.1f \mu s$'%(fitted_pars[0][3]))
#                 if fitted_pars[1][1] < 0.01:
#                     fitted_pars[1][1] = 0
#                 ax1.plot(x_vector[1],ramsey(x_vector[1], fitted_pars[1][0], fitted_pars[1][1], fitted_pars[1][2],fitted_pars[1][3],fitted_pars[1][4]),'g',linewidth=linewidth,label='$\omega = 2\pi\\times%.1f$ MHz, $\\tau/\\tau_0$ = %.1f $\mu$s'%(fitted_pars[1][1],fitted_pars[1][3]/fitted_pars[0][3]))
#                 ax1.legend(loc='upper right',prop={'size':22})
#                 # ax1.set_xlim([0,x_vector[0][-1]])
#                 textstr = '$T_{\pi/2}$=%.1f ns\n$\omega_d$ = %.5f GHz\n$A_d$ = %.3f V\n$\Delta$=%.3f MHz\n$\\tau$/$\\tau_0$=%.1f\n$\mu$ = %d mV\n$\omega_{AC}$ = %.4f GHz\n$\sigma$ = %.1f mV\n$\hatn$ = %d'%(pi2Width*1e9,qubitDriveFreq*1e-9,amplitude_hd,fitted_pars[0][1],fitted_pars[1][3]/fitted_pars[0][3],float(AC_pars[1:AC_pars.find(' ')])*1e3,AC_freq*1e-9,float(AC_pars[AC_pars.find(' ')+1:])*1e3,nAverages)
#                 # plt.gcf().text(0.95, 0.15, textstr, fontsize=fontSize)
#             elif (sweep == 0 or sweep == 1) and fit_single_par_point == 0: #for plotting non-par sweep data and par sweep data during the sweep
#                 fig = plt.figure(figsize=(20,18))
#                 ax1 = fig.add_subplot(2,2,(1,2))
#                 # Plot data
#                 fontSize = 36
#                 tickSize = 24
#                 markersize= 12
#                 linewidth = 6
#                 ax1.plot(x_vector, y_vector, 'o', markersize = markersize, c='C0')
#                 ax1.set_ylabel('Digitizer Voltage (mV)',fontsize=fontSize)
#                 for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
#                     label.set_fontsize(tickSize)
#                 # Retrieve and plot telegraph noise instrance
#                 waveforms = ziut.parse_awg_waveform(awg.get('/dev8233/awgs/0/waveform/waves/0')['dev8233']['awgs']['0']['waveform']['waves']['0'][0]['vector'],channels=2)
#                 t_arr = np.linspace(0,x_vector[-1],len(waveforms[0]))
#                 ax2 = fig.add_subplot(2,2,3)
#                 ax2.plot(t_arr,1e3*waveforms[0],linewidth=linewidth,color='b')
#                 ax3 = fig.add_subplot(2,2,4)
#                 ax3.plot(t_arr,1e3*waveforms[1],linewidth=linewidth,color='r')
#                 for label in (ax2.get_xticklabels() + ax2.get_yticklabels()):
#                     label.set_fontsize(tickSize)
#                 for label in (ax3.get_xticklabels() + ax3.get_yticklabels()):
#                     label.set_fontsize(tickSize)
#                 ax2.set_xlabel('Pulse Separation ($\mu$s)',fontsize=fontSize)
#                 ax2.set_ylabel('Amplitude (mV)',fontsize=fontSize)
#                 ax2.set_title('$\sigma_x$ waveform',fontsize=fontSize)
#                 ax3.set_xlabel('Pulse Separation ($\mu$s)',fontsize=fontSize)
#                 ax3.set_title('$\sigma_z$ waveform',fontsize=fontSize)
#                 if tauk != 0:
#                     textstr = '$T_{\pi/2}$=%.1f ns\n$\omega_d$ = %.6f GHz\n$A_d$ = %.3f V\n$\Delta$=%.3f MHz\n$T_2^*$=%.2f $\mu$s\n$\mu$ = %d mV\n$\omega_{AC}$ = %.4f GHz\n$\sigma$ = %.1f mV\n$B_0$ = %.2f mV\n$\\tau_k$ = %.3f $\mu s$\n$\\nu$ = %.1f MHz\n$\hatn$ = %d'%(pi2Width*1e9,qubitDriveFreq*1e-9,amplitude_hd,fitted_pars[1],fitted_pars[3],AC_pars[0]*1e3,AC_freq*1e-9,AC_pars[1]*1e3,RT_pars[0]*1e3,RT_pars[1],RT_pars[2],nAverages)
#                 elif RT_pars[2] == 0:
#                     textstr = '$T_{\pi/2}$=%.1f ns\n$\omega_d$ = %.6f GHz\n$A_d$ = %.3f V\n$\Delta$=%.3f MHz\n$T_2^*$=%.2f $\mu$s\n$\mu$ = %d mV\n$\omega_{AC}$ = %.4f GHz\n$\sigma$ = %.1f mV\n$B_0$ = %.2f mV\n$\\tau_k$ = %.3f $\mu s$\n$\hatn$ = %d'%(pi2Width*1e9,qubitDriveFreq*1e-9,amplitude_hd,fitted_pars[1],fitted_pars[3],AC_pars[0]*1e3,AC_freq*1e-9,AC_pars[1]*1e3,RT_pars[0]*1e3,RT_pars[1],nAverages)
#                 if plot_mode == 0:
#                     # for plotting single realization
#                     ax1.plot(x_vector,ramsey(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2],fitted_pars[3],fitted_pars[4]),'r',linewidth=linewidth)
#                     plt.gcf().text(0.925, 0.25, textstr, fontsize=fontSize,bbox=dict(boxstyle='round,rounding_size=1.25',facecolor='silver',alpha=0.5))
#                     ax1.set_title('Ramsey Measurement %03d' %(iteration),size=fontSize)
#                 elif plot_mode == 1:
#                     textstr = '$T_2^*$=%.2f $\mu$s\n$\sigma$ = %.1f mV'%(fitted_pars[3],AC_pars[1]*1e3)
#                     plt.gcf().text(0.925, 0.45, textstr, fontsize=fontSize+10,bbox=dict(boxstyle='round,rounding_size=1.25',facecolor='silver',alpha=0.5))
#                 elif plot_mode == 2:
#                     textstr = '$\omega$=%.3f MHz\n$T_2^*$=%.2f $\mu$s\n$\sigma$ = %.1f mV\n$B_0$ = %.2f mV\n$\\tau_k$ = %.3f $\mu s$'%(fitted_pars[1],fitted_pars[3],AC_pars[1]*1e3,RT_pars[0]*1e3,RT_pars[1])
#                     plt.gcf().text(0.925, 0.45, textstr, fontsize=fontSize+10,bbox=dict(boxstyle='round,rounding_size=1.25',facecolor='silver',alpha=0.5))
#         elif AC_pars[0] == 0 and RT_pars[0] == 0 and AC_pars[1] == 0: # plot data without any noise instances
#             fig = plt.figure()
#             ax1 = fig.add_subplot(111)
#             fontSize = 16
#             tickSize = 11
#             markersize= 6
#             linewidth = 2
#             ax1.plot(x_vector, y_vector, 'o', markersize = markersize, c='C0')
#             ax1.set_ylabel('Digitizer Voltage (mV)',fontsize=fontSize)
#             ax1.set_xlabel('Pulse Separation ($\mu$s)')
#             for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
#                 label.set_fontsize(tickSize)
#             ax1.plot(x_vector,ramsey(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2],fitted_pars[3],fitted_pars[4]),'r',linewidth=linewidth)
#             textstr = '$T_{\pi/2}$=%.1f ns\n$\omega_d$ = %.6f GHz\n$A_d$ = %.3f V\n$\Delta$=%.3f MHz\n$T_2^*$=%.2f $\mu$s\n$\mu$ = %d mV\n$\omega_{AC}$ = %.4f GHz\n$\sigma$ = %.1f mV\n$B_0$ = %.1f mV\n$\\tau_k$ = %.3f $\mu s$\n$\hatn$ = %d'%(pi2Width*1e9,qubitDriveFreq*1e-9,amplitude_hd,fitted_pars[1],fitted_pars[3],AC_pars[0]*1e3,AC_freq*1e-9,AC_pars[1]*1e3,RT_pars[0]*1e3,RT_pars[1],nAverages)
#             plt.gcf().text(0.95, 0.15, textstr, fontsize=fontSize)
#             ax1.set_title('Ramsey Measurement %03d' %(iteration),size=fontSize)

#     elif sequence == "echo":
#         fig = plt.figure()
#         ax = fig.add_subplot(111)
#         ax.plot(x_vector, y_vector, '-o', markersize = 3, c='C0')
#         ax.set_ylabel('Digitizer Voltage (mV)')
#         ax.set_xlabel('Pulse Separation ($\mu$s)')
#         ax.plot(x_vector,decay(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2]),'r')
#         textstr = '$T_{\pi/2}$=%.1f ns\n$T_{\pi(Y)} = %.1f$ ns\n$\omega_d$ = %.4f GHz\n$A_d$ = %.2f V\n$T_2$=%.2f $\mu$s\n$\mu$ = %.3f V\n$\omega_{AC}$ = %.4f GHz\n$\sigma$ = %.3f V\n$B_0$ = %.3f V\n$\\tau_k$ = %.2f $\mu s$\n$\hatn$ = %d'%(pi2Width*1e9,piWidth_Y*1e9,qubitDriveFreq*1e-9,amplitude_hd,fitted_pars[1],AC_pars[0],AC_freq*1e-9,AC_pars[1],RT_pars[0],RT_pars[1],nAverages)
#         ax.set_title('Echo Measurement %03d' %(iteration))
#         plt.gcf().text(0.95, 0.15, textstr, fontsize=14)

#     elif sequence == "T1":
#         fig = plt.figure()
#         ax = fig.add_subplot(111)
#         ax.plot(x_vector, y_vector, '-o', markersize = 3, c='C0')
#         ax.set_ylabel('Digitizer Voltage (mV)')
#         ax.set_xlabel('Delay ($\mu$s)')
#         ax.plot(x_vector,decay(x_vector, fitted_pars[0], fitted_pars[1], fitted_pars[2]),'r')
#         textstr = '$T_{\pi/2}$=%.1f ns\n$\omega_d$ = %.4f GHz\n$A_d$ = %.2f V\n$T_1$=%.2f $\mu$s\n$\mu$ = %.3f V\n$\omega_{AC}$ = %.4f GHz\n$\sigma$ = %.3f V\n$B_0$ = %.2f V\n$\\tau_k$ = %.2f $\mu s$\n$\hatn$ = %d'%(pi2Width*1e9,qubitDriveFreq*1e-9,amplitude_hd,fitted_pars[1],AC_pars[0],AC_freq*1e-9,AC_pars[1],RT_pars[0],RT_pars[1],nAverages)
#         ax.set_title('T1 Measurement %03d' %(iteration))
#         plt.gcf().text(0.95, 0.15, textstr, fontsize=14)

#     return fig

# def fit_data(x_vector,y_vector,sequence='rabi',dt=0.01,qubitDriveFreq=3.8e9,amplitude_hd=1,AC_pars=[0,0],RT_pars=[0,0],
#                  pi2Width=50e-9,piWidth_Y=0,fitting = 1, plot=1,save_fig = 0,iteration=1,nAverages=1,sampling_rate=2.4e9,
#                  integration_length=2e-6,AC_freq=5e9,cav_resp_time=5e-6,sweep=0,stepSize=5e-6,Tmax=5e-6,measPeriod=5e-6,nSteps=101,
#                  prePulseLength=1500e-9,postPulseLength=1500e-9,threshold=10e-3,active_reset=False,rr_IF=5e6,source=0):

#     '''
#     fit experiment data

#     sequence:          'Rabi', 'T1' or 'T2'
#     complex_amplitude:  complex amplitude from the measurement
#     x_vector:           x data
#     fitting:     0:     do not fit
#                   1:     do fit
#     save_fig:    0:     do not save
#                   1:     do save
#     '''
#     x_vector = x_vector*1e6
#     abs_camplitude = y_vector*1e3

#     if sequence == "rabi":
#         amp = (max(abs_camplitude)-min(abs_camplitude))/2
#         offset = np.mean(abs_camplitude)
#         period = 1e3/(extract_freq(x_vector*1e3, abs_camplitude, dt,plot=0))
#         # period = 110
#         print('Period Initial Guess: %.1f ns'%(period))
#         phase = 0
#         p0 = [amp,period,phase,offset]
#         fitted_pars, covar = scy.optimize.curve_fit(rabi, x_vector*1e3, abs_camplitude,p0=p0,xtol=1e-6,maxfev=3000)
#         pi_pulse = np.round(1/2*fitted_pars[1])
#         error = np.sqrt(abs(np.diag(covar)))
#         print("Pi pulse duration is %.1f ns"%(pi_pulse))

#         # return pi_pulse,error

#     elif sequence == "ramsey":
#         amp = abs_camplitude[0]-abs_camplitude[-1]
#         offset = np.mean(abs_camplitude)
#         f = extract_freq(x_vector, abs_camplitude,dt)
#         print('Initial Guess for Freq:%.4f MHz'%(f))
#         if x_vector[-1] > 10:
#             tau = 10
#         else:
#             tau = 0.2
#         phase = 0
#         p0 = [amp,f,phase,tau,offset]
#         try:
#             fitted_pars, covar = scy.optimize.curve_fit(ramsey, x_vector, abs_camplitude,p0=p0,xtol=1e-6,maxfev=6000)
#             detuning = fitted_pars[1]
#             T_phi = fitted_pars[3]
#             error = np.sqrt(abs(np.diag(covar)))
#         except:
#             print('fitting failed')
#             fitted_pars = np.zeros(5)
#             detuning = 0
#             T_phi = 0
#             error = 20*np.ones(5)

#         # return detuning,T_phi,error

#     elif sequence == "echo":
#         amp = abs_camplitude[0]-abs_camplitude[-1]
#         offset = np.mean(abs_camplitude)
#         if x_vector[-1] < 2:
#             tau = 0.6
#         else:
#             tau = 1.5
#         p0 = [amp,tau,offset]
#         try:
#             fitted_pars, covar = scy.optimize.curve_fit(decay, x_vector, abs_camplitude,p0=p0,xtol=1e-6,maxfev=6000)
#             T_2 = fitted_pars[1]
#             error = np.sqrt(abs(np.diag(covar)))
#             # if sum(error) > 2: # try fitting again, this time excluding the first few points
#             #     fitted_pars, covar = scy.optimize.curve_fit(decay, x_vector[10:], abs_camplitude[10:],p0=p0,xtol=1e-6,maxfev=6000)
#             #     T_2 = fitted_pars[1]
#             #     error = np.sqrt(abs(np.diag(covar)))
#         except:
#             print('fitting failed')
#             fitted_pars = np.zeros(3)
#             T_2 = 0
#             error = 20*np.ones(3)

#         # return T_2,error

#     elif sequence == "T1":
#         amp = abs_camplitude[0]-abs_camplitude[-1]
#         offset = np.mean(abs_camplitude)
#         tau = 2
#         p0 = [amp,tau,offset]
#         try:
#             fitted_pars, covar = scy.optimize.curve_fit(decay, x_vector, abs_camplitude,p0=p0,xtol=1e-6,maxfev=3000)
#             T_1 = fitted_pars[1]
#             error = np.sqrt(abs(np.diag(covar)))
#         except:
#             print('fitting failed')
#             fitted_pars = np.zeros(3)
#             T_1 = 0
#             error = 20*np.ones(3)

#     return fitted_pars,error
# def fit_beats(sequence,x_vector,y_vector,dt=0.01,qubitDriveFreq=3.8e9,amplitude_hd=1,AC_pars=[0,0],RT_pars=[0,0],pi2Width=50e-9,fitting = 1, plot=1,save_fig = 0,iteration=1):
#     x_vector = x_vector*1e6
#     abs_camplitude = np.abs(y_vector*1e3)
#     amp = abs_camplitude[0]-abs_camplitude[-1]
#     offset = np.mean(abs_camplitude)
#     f1 = 15
#     f2 = 1
#     tau = 15
#     phi1 = 0
#     phi2 = 0
#     p0 = [amp,f1,f2,phi1,phi2,tau,offset]
#     lb = [-1000,-10,-10,-np.pi,-np.pi,0,-np.inf]
#     ub = [1000,20,20,np.pi,np.pi,30,np.inf]
#     fitted_pars, covar = scy.optimize.curve_fit(beats, x_vector, abs_camplitude,p0=p0,bounds=(lb,ub),xtol=1e-12,maxfev=6000)
#     f1 = best_vals[1]
#     f2 = best_vals[2]
#     T_phi = best_vals[5]
#     error = np.sqrt(abs(np.diag(covar)))
#     fig = plt.figure()
#     ax = fig.add_subplot(111)
#     ax.plot(x_vector, abs_camplitude, '-o', markersize = 3, c='C0')
#     ax.set_ylabel('Digitizer Voltage (mV)')
#     ax.set_xlabel('Pulse Separation ($\mu$s)')
#     ax.plot(x_vector,beats(x_vector,best_vals[0],best_vals[1],best_vals[2],best_vals[3],best_vals[4],best_vals[5],best_vals[6]),'r')
#     textstr = '$A$ = %.2f mV\n$\omega_1$=%.2f MHz\n$\omega_2$=%.2f MHz\n$T_2^*$=%.2f $\mu$s\n$B_0$ = %.2f V\n$\\tau_k$ = %.2f $\mu s$'%(best_vals[0],best_vals[1],best_vals[2],best_vals[5],RT_pars[0],RT_pars[1])
#     ax.set_title('Ramsey Measurement %03d' %(iteration))
#     plt.gcf().text(1, 0.25, textstr, fontsize=14)