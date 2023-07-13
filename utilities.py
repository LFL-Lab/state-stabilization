# -*- coding: utf-8 -*-
"""
Created on Thu May 11 14:49:10 2023

@author: Evangelos
"""

import time
import math
import numpy as np
import matplotlib.pyplot as plt

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
    y = int(base*round(n_points/base))
    if y==0:
        y = int(base*round(n_points/base+1))
        
    return y

def Volt2dBm(data):

    return 10*np.log10(1e3*data**2/50)

def Watt2dBm(x):
    '''
    converts from units of Watts to dBm
    '''
    return 10.*np.log10(x*1000.)

def gen_arb_wfm(wfm_type,wfm_pars,channel='I'):
    
    time_arr = np.arange(wfm_pars['t0'],wfm_pars['tmax']-wfm_pars['dt']/2,wfm_pars['dt'])
    
    if wfm_type == 'rising':
        fun = lambda x : (1/np.sqrt(wfm_pars['tb'] - x)) 
        wfm_I = wfm_pars['amp']*fun(time_arr)
        wfm_Q = wfm_pars['amp']*fun(time_arr)
        plt.plot(time_arr,wfm)
    elif wfm_type == 'markov':
        wfm = np.random.normal(wfm_pars['mu'], wfm_pars['sigma'], wfm_pars['n_points'])
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
    t0 = roundToBase(pars['fsAWG']*pars['x0'])
    dt = roundToBase(pars['fsAWG']*pars['dx'])
    n_points = roundToBase(pars['xmax']*pars['fsAWG']) # this ensures there is an integer number of time points
    n_steps = int((n_points-t0)/dt)+1 # 1 is added to include the first point
    tmax = dt*n_steps
   
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
    v_b = np.zeros((1,3))
    
    for i in range(len(v_b)):
        v_b[i] = 1-2*(data[i]-calib_pars[i])/calib_pars[3]
        
    return v_b

def compute_rho(vb):
    
    iden = [[1,0],[0,1]]
    sx = [[0,1],[1,0]]
    sy = [[0,-1j],[1j,0]]
    sz = [[1,0],[0,-1]]
    rho = 1/2*(iden+np.multiply(vb[0][0],sx)+np.multiply(vb[0][1],sy)+np.multiply(vb[0][2],sz))
    
    return rho

def normalize_data(data):
    
    offset = np.mean(data[2,:])
    amp = (max(data[2,:])-min(data[2,:]))/2
    
    normalized_data = (data-offset)/amp
    
    return normalized_data