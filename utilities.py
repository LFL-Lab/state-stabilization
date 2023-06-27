# -*- coding: utf-8 -*-
"""
Created on Thu May 11 14:49:10 2023

@author: Evangelos
"""

import time
import math

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