# -*- coding: utf-8 -*-
"""
Created on Thu May 11 14:49:10 2023

@author: Evangelos
"""

import time
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
# from plotting import extract_freq

# Retry decorator with exponential backoff
def retry(tries, delay=3, backoff=2):
    
    '''Retries a function or method until it returns True.
     delay sets the initial delay in seconds, and backoff sets the factor by which
    the delay should lengthen after each failure. backoff must be greater than 1,
    or else it isn't really a backoff. tries must be at least 0, and delay
    greater than 0.'''

    if backoff <= 1:
        raise ValueError("backoff must be greater than 1")
    tries = math.floor(tries)
    if tries < 0:
       raise ValueError("tries must be 0 or greater")
    if delay <= 0:
       raise ValueError("delay must be greater than 0")
    def deco_retry(f):
         
        def f_retry(*args, **kwargs):
            mtries, mdelay = tries, delay # make mutable
            rv = f(*args, **kwargs) # first attempt
            while mtries > 0:
                if rv is True: # Done on success
                    return True
                mtries -= 1      # consume an attempt
                time.sleep(mdelay) # wait...
                mdelay *= backoff  # make future wait longer
                rv = f(*args, **kwargs) # Try again
            return False # Ran out of tries :-(
        return f_retry # true decorator -> decorated function
    return deco_retry  # @retry(arg[, ...]) -> true decorator

def odd(n):
    return range(1,n,2)

def even(n):
    return range(0,n,2)

def roundToBase(n_points,base=16):
    '''Make the AWG happy by uploading a wfm whose points are multiple of 16'''
    y = int(base*math.floor(n_points/base))
    if y==0:
        y = int(base*math.floor(n_points/base+1))
        
    return y

def Volt2dBm(data):

    return 10*np.log10(1e3*data**2/50)

def Watt2dBm(x):
    '''
    converts from units of Watts to dBm
    '''
    return 10.*np.log10(x*1000.)

#def gen_arb_wfm(wfm_type,wfm_pars,normalize=False,n_points=1024):
def gen_arb_wfm(wfm_type,wfm_pars={},exp_pars={},qb_pars={},normalize=False,n_points=1024):
    '''a function that generates different types of arbitrary waveforms, 
    based on the wfm parameters
    INPUT:
    -----
    OUTPUT:
    -----
    '''
    exp = exp_pars['exp']
    if exp == 'coherence-stabilization':
    #
        time_arr = np.arange(wfm_pars['t0'],wfm_pars['tmax'],wfm_pars['dt'])
    else:
        # = np.arange(wfm_pars['t0'],wfm_pars['tmax']-wfm_pars['dt']/2,wfm_pars['dt'])
    #print(len(time_arr), n_points)
        time_arr = np.linspace(wfm_pars['t0'],wfm_pars['tmax'],n_points)
    if wfm_type == 'rising':
        gamma = 1/(2*wfm_pars['T2'])
        #gamma = np.pi/(2*wfm_pars['T2'])
        wfm = compute_wfm(time_arr,
                          exp_pars=exp_pars,
                          qb_pars=qb_pars,
                          gamma=gamma,
                          plot = True)
        print('\n rising is',len(time_arr))
        #wfm = compute_wfm(time_arr,gamma)
    elif wfm_type == 'markov':
        wfm_long = np.random.normal(wfm_pars['mu'], wfm_pars['sigma'], 2*len(time_arr))
        wfm = wfm_long[int(len(time_arr)/2):3*int(len(time_arr)/2)]
        Nf = 1/2*wfm_pars["fsAWG"]
        # print(f'Nyquist frequency is {Nf*1e-6:.1f} MHz')
        fc = 50e6
        Wn = fc/Nf
        b,a = signal.butter(8,Wn,btype='lowpass')
        wfm = signal.lfilter(b, a, wfm)
        # extract_freq(time_arr, wfm, wfm_pars['dt']*1e6,plot=1)
    elif wfm_type == 'telegraph':
        wfm = np.cos(2*np.pi*wfm_pars['nu']*1e3*time_arr + 2*np.pi*np.random.rand()) * gen_tel_noise(len(time_arr), wfm_pars['tau'], dt = wfm_pars['tmax']/len(time_arr))

    return wfm

def gen_tel_noise(numPoints,tau,dt):
    '''Generates a single instance of telegraph noise'''
    signal = np.ones(numPoints)*(-1)**np.random.randint(0,2)
    for i in range(1,numPoints-1):
        if np.random.rand() < 1/(2*tau*1e-6/dt)*np.exp(-1/(2*tau*1e-6/dt)):
            signal[i+1] = - signal[i]
        else:
            signal[i+1] = signal[i]
    return signal

def calc_steps(pars,verbose=True):
    """
    Calculates the number of steps in the sequence and number of points in the waveform used in the experiment. The final stepsize of the sequence
    might be different than the one specified due to the fact that the AWG has a granularity of the waveform is 16 samples.

    Args:
        sequence (string, optional): Type of experiment (rabi,ramsey,T1, echo). Defaults to 'ramsey'.
        fsAWG (float, optional): HDAWG sampling rate. Defaults to 1.2e9.
        stepSize (float, optional): dt for sequence. Defaults to 10e-9.
        Tmax (float, optional): maximum experiment array length. For example if Tmax = 5e-6 for ramsey experiment then maximum pulse separation is 5us. Defaults to 5e-6.
        verbose (boolean, optional):  Defaults to 1.

    Returns:
        n_points (int): number of waveform points.
        n_steps (int): number of sequence steps.
        t0 (int): sequence step size in units of samples 1/fsAWG.
        dt (int): starting sequence point in units of samples.

    """
    
    # convert inital time step to samples (integer multiple of 16)
    t0 = roundToBase(pars['fsAWG']*pars['x0'])
    # convert time intervals into samples (integer multiple of 16)
    dt = roundToBase(pars['fsAWG']*pars['dx'])
    # compute number of samples required for the full time (integer multiple of 16)
    n_points = roundToBase((pars['xmax']-pars['x0'])*pars['fsAWG']) # this ensures there is an integer number of time points
    n_steps = int(n_points/dt) 
    
    tmax = dt*n_steps + t0
    #tmax = roundToBase((pars['xmax'])*pars['fsAWG']) #this is the same value as tmax
    
    if verbose:
        print("dt is %.1f ns (%d pts) ==> f_s = %.1f MHz \nn_points = %d | n_steps is %d | Pulse length start = %.1f ns (%d pts)" %(dt/pars['fsAWG']*1e9,dt,1e-6*pars['fsAWG']/dt,n_points,n_steps,t0*1e9/pars['fsAWG'],t0))
    else:
        pass
    
    if n_steps > 1024:
        raise Exception('Error: The maximum number of steps is 1024')

    return n_points,n_steps,t0,tmax,dt

def generate_xarray(pars):
    '''Generate xarray parameters for power-rabi and virtual z-gate experiments'''
    x0 = pars['x0']
    xmax = pars['xmax']
    dx = pars['dx']
    x_array = np.arange(x0,xmax+dx/2,dx)
    n_steps = len(x_array)

    return x0,xmax,dx,x_array,n_steps

def compute_bloch(data,calib_pars):
    '''computes bloch vector components given the calibrated voltages'''
    #print(data)
    v_b = np.zeros((3,data.shape[1]))
    
    for j in range(data.shape[1]):
        #v_b[0,j] = 1-2*(calib_pars[0]-data[0,j])/calib_pars[3]
        #v_b[1,j] = 1-2*(calib_pars[1]-data[1,j])/calib_pars[4]
        #v_b[2,j] = 1-2*(calib_pars[2]-data[2,j])/calib_pars[5]
        for i in range(3):
            v_b[i,j] = 1-2*(calib_pars[2]-data[i,j])/calib_pars[3]
        
    return v_b

def compute_rho(vb):
    
    iden = [[1,0],[0,1]]
    sx = [[0,1],[1,0]]
    sy = [[0,-1j],[1j,0]]
    sz = [[1,0],[0,-1]]
    rho = 1/2*(iden+np.multiply(vb[0],sx)+np.multiply(vb[1],sy)+np.multiply(vb[2],sz))
    
    return rho

def compute_coherence(vb,calib_states,plane='ZX'):
    
    v1 = []
    v2 = []
    coherence = []
    
    for i in range(vb.shape[1]):
        # print(np.array(calib_states)*1e3)
        # vb = compute_bloch(data[:,i], calib_states)
        if plane == 'XY':
            v1.append(vb[0,i])
            v2.append(vb[1,i])
        elif plane == 'ZX':
            v1.append(vb[2,i])
            v2.append(vb[0,i])
        coherence.append(v1[i]**2+v2[i]**2)
    
    return np.array(coherence)

def compute_purity(v_b,calib_states):
    
    purr = []
    
    for i in range(v_b.shape[1]):
        # vb = compute_bloch(data[:,i], calib_states)
        purr.append(1/2*(1+np.linalg.norm(v_b[:,i])**2))
        
    return np.array(purr)

def normalize_data(data):
    
    offset = np.mean(data[2,:])
    amp = (max(data[2,:])-min(data[2,:]))/2
    
    normalized_data = (data-offset)/amp
    
    return normalized_data

# def compute_wfm(time_arr,gamma,plot=True):
#     '''
#     computes amplitude waveform based on input value of decoherence rate gamma

#     Parameters
#     ----------
#     time_arr : TYPE
#         DESCRIPTION.
#     gamma : decoherence rate in Hz
#     plot : TYPE, optional
#         DESCRIPTION. The default is True.

#     Returns
#     -------
#     TYPE
#         DESCRIPTION.

#     '''
#     wfm = np.zeros(len(time_arr))
#     tb = 1/(4*gamma)
#     for i in range(len(time_arr)):
#         value = -gamma/(np.sqrt(1 - time_arr[i]/tb))
#         if math.isnan(value):
#             wfm[i:] = 0
#             #wfm[i:] = wfm[i-1]
#             break
#         else:
#             # print(value*1e-6)
#             wfm[i] = convert_w_to_v(value)
    
#     # print(wfm)
#     if plot:
#         plt.plot(time_arr,wfm,time_arr,np.full(len(time_arr),tb))
#         plt.ylim([-1,1])
        
#     return np.array(wfm)

def compute_wfm(time_arr,gamma,exp_pars={},qb_pars={},plot=True):
    '''
    computes amplitude waveform based on input value of decoherence rate gamma

    INPUT
    ----------
    time_arr (ARRAY): time array for waveform.
    gamma (FLOAT): decoherence rate in Hz.
    exp_pars (DICTIONARY): list of experimental parameters defined in Master Script before running.
    qb_pars (DICTIONARY): list of qubit parameters defined in JSON file.
    plot (BOOLEAN): optional, default set to TRUE


    OUTPUT
    -------
    wfm (ARRAY): array of waveform amplitude at each timestep in time array

    '''
    
    # determine polar and azimuthal axes
    theta = exp_pars['initial-state'][1] * np.pi # POLAR
   
    phi = exp_pars['initial-state'][0] * np.pi # AZIMUTHAL
    
    #Compute vz and vy (Bloch vector components)
    vy = np.sin(theta) * np.sin(phi)
    vz = np.cos(theta)
    
    
    wfm = np.zeros(len(time_arr))
    
    # Compute Breakdown Time
    tb = vy**2 / (4*gamma* vz**2)
    
    for i in range(len(time_arr)):
        
        #compute amplitude
        #value = -gamma/(np.sqrt(1 - time_arr[i]/tb))
        #value = 0.000001*-1/np.pi * np.sqrt(gamma /(tb-time_arr[i]))
        value = -1/2 * np.sqrt(gamma /(tb-time_arr[i]))
        if math.isnan(value):
            #turn off wave when t = tb
            wfm[i:] = 0
            #wfm[i:] = wfm[i-1]
            break
        else:
            # print(value*1e-6)
           # wfm[i] = convert_w_to_v(value)
           wfm[i] = convert_w_to_v(w=value,qb_pars = qb_pars,b = 0)
    # print(wfm)
    if plot:
        plt.plot(time_arr,wfm,time_arr,np.full(len(time_arr),tb))
        plt.ylim([-1,1])
        
    return np.array(wfm)


def convert_w_to_v(w,qb_pars = {},b = 0):
#   convert_w_to_v(w,a=7.3236,b=0):
    '''converts input from angular units to amplitude in volts
    
    INPUT
    ------
    w (FLOAT):
    #a (FLOAT): 
    b (FLOAT):
    OUTPUT
    ------
    '''
    #a is the Rabi frequency (in MHz) of a 1 V signal
    #to find it, take 1 / (2*pi_len/fs * pi_amp)
    
    a = compute_a(qb_pars = qb_pars)
    return (w*1e-6-b)/(2*np.pi*a)

def line(x,a,b):
    return a*x+b

def sort_data(unsorted_data):
    
    L = len(unsorted_data)
    # print(L)
    sorted_data = np.zeros((3,int(L/3)))
    for i in range(L):
        j = i % 3
        k = int(i/3)
        sorted_data[j,k] = unsorted_data[i]
    
    return sorted_data


def compute_a(qb_pars = {}):  
    ''' Computes a that is used in convert_w_to_v'''
    
    return (2 * qb_pars['pi_len']/ 2.4 *1e-3  * qb_pars['pi_amp'])**(-1)
              