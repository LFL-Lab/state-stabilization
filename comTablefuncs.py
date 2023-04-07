#!/usr/bin/env python
# coding: utf-8

# # Import

# In[58]:


import textwrap
# import zhinst.qcodes
import json
import urllib
import jsonschema
import zhinst.ziPython as zp
import zhinst.utils as zu
import zhinst.toolkit as zt
import time
import numpy as np
import UHFQA as qa


# # Predefined functions

# In[2]:

# Please refer to https://docs.zhinst.com/hdawg/commandtable/v2/schema for other settings


def ct_pulse_length(n_wave, pulse_length_start = 32, pulse_length_increment = 16,mu=0,AC_pre_length=0,
                    AC_post_length=0,sequence='rabi',pipulse=50,active_reset=False,noise_rate=0):

    ct = {'header':{'version':'0.2'}, 'table':[]}
    # make table entries for waveform playback during free evolution
    if sequence == 'ramsey':
        pulse_length_start = 0
        for i in range(n_wave-1):
            sample_length = pulse_length_start + (i+2) * pulse_length_increment
            entry = {'index': i,
                      'waveform':{
                          'index': 0,
                          'length':  sample_length,
                            # 'samplingRateDivider': noise_rate,
                      },
                    }
            ct['table'].append(entry)
            i += 1
    else:
        for i in range(n_wave):
            sample_length = pulse_length_start + i * pulse_length_increment
            entry = {'index': i,
                      'waveform':{
                          'index': 0,
                          'length':  sample_length,
                          # 'samplingRateDivider': noise_rate,
                      },
                    }
            ct['table'].append(entry)
            i += 1

    # make pre and post pulses in the case of AC stark noise added to the system
    if mu != 0:
        entry = {'index': n_wave,
                  'waveform':{
                      'index': 1,
                      'length': AC_pre_length
                  },
                  }
        ct['table'].append(entry) # ac pre-pulse
        n_wave += 1

        if sequence != 'rabi':
            entry = {'index': n_wave,
                      'waveform':{
                          'index': 2,
                          'length': AC_post_length
                          },
                      }
            ct['table'].append(entry) # ac post-pulse
            n_wave += 1

            if sequence == 'echo':
                entry = {'index': n_wave,
                          'waveform':{
                              'index': 3,
                              'length': pipulse
                              },
                          }
                ct['table'].append(entry) # ac mid-pulse
                n_wave += 1

        if active_reset == True:
            reset_pulse_len = AC_pre_length + AC_post_length
            entry = {'index': n_wave,
                     'waveform':{
                         'index': 3,
                         'length': reset_pulse_len
                         },
                     }
            ct['table'].append(entry) # reset pulse
            n_wave += 1

    if mu == 0:
        # if sequence == 'echo':
        #     entry = {'index': n_wave,
        #               'waveform':{
        #                   'index': 3,
        #                   'length': mid_pulse_length
        #                   },
        #               }
        #     ct['table'].append(entry) # ac mid-pulse
        entry = {'index': n_wave,
                  'waveform':{
                      'index': 1,
                      'length': int(pipulse/2),
                  },
                  }
        ct['table'].append(entry) # ac pre-pulse
        n_wave += 1

        if active_reset == True:
            entry = {'index': n_wave,
                     'waveform':{
                         'index': 2,
                         'length': pipulse,
                         },
                     }
            ct['table'].append(entry) # reset pulse
            n_wave += 1
    return ct



# Please refer to https://docs.zhinst.com/hdawg/commandtable/v2/schema for other settings
def ct_amplitude_increment(amplitude_start=[0.5,0.5],amplitude_increment = 0.001*np.ones(2)):
    #first entry to just set amplitudes to inital value, no waveform necessary (but wanted here)

    ct = {'header':{'version':'0.2'}, 'table':[]}

    entry = {'index': 0,
               # 'waveform':{
               #     'index': 0
               # },
             'amplitude0':{
                 'value':amplitude_start[0],
                 'increment': False
             },
             'amplitude1':{
                 'value':amplitude_start[1],
                 'increment': False
             }
            }
    ct['table'].append(entry)

    # second entry defines increment and waveform
    entry = {'index': 1,
              'waveform':{
                       'index': 0,
                       'length': 1024
              },
             'amplitude0':{
                 'value':amplitude_increment[0],
                 'increment': True
             },
             'amplitude1':{
                 'value':amplitude_increment[1],
                 'increment': True
             }
            }

    ct['table'].append(entry)
    return ct


# # Run AWG with different pulse length

def example(daq,awg,device_awg):# change pulse length from pulse_length_start to pulse_length_start + n_wave * pulse_length_increment
    pulse_length = 2**14
    pulse_length_start = 320 #
    pulse_length_increment = 1600 # waveform granularity
    n_wave = 10

    if n_wave*pulse_length_increment+pulse_length_start>pulse_length:
        n_wave = np.floor((pulse_length-pulse_length_start)/pulse_length_increment)
        print('n_wave is rounded')

    awg_program = f'''
        const samples = {pulse_length};
        wave w = gauss(samples, 1, samples/2, samples/8); // 'w' can also be '.csv' file
        wave w1 = gauss(samples, 0.5, samples/2, samples/8);
        assignWaveIndex(w, w1, 0);        //For dual-channel waveforms
        var n;
        for(n=0; n < {n_wave}; n++)
        {{
        executeTableEntry(n);
        playZero({pulse_length}*2-{pulse_length_increment}*n); // remain same cycle length
        }}
        '''
    # compile seqc code
    qa.create_and_compile_awg(daq, device_awg, awg_program, 0, 1)
    # generate json file
    ct=ct_pulse_length(n_wave=n_wave, pulse_length_start = pulse_length_start, pulse_length_increment = pulse_length_increment)
    # upload json file
    daq.setVector(f"/{device_awg}/awgs/0/commandtable/data", json.dumps(ct))

    # turn on wave output 1 and 2
    for i in range(2):
        daq.setInt(f'/{device_awg}/sigouts/{i}/on', 1)
    # run AWG
    daq.setInt(f'{device_awg}/awgs/0/enable', 1)
#    print('enable hd_awg')
    time.sleep(20)
    daq.setInt(f'{device_awg}/awgs/0/enable', 0)
#    print('disable hd_awg')
    # get what json file have been loaded
    json.loads(daq.get(f"/{device_awg}/awgs/0/commandtable/data",flat=True)[f"/{device_awg}/awgs/0/commandtable/data"][0]['vector']) #


    # In[103]:


    # change pulse length from pulse_length_start to pulse_length_start + n_wave * pulse_length_increment
#    pulse_length = 1024
#    n_wave = 100
#    amplitude_start = [0.5, 0.5]
#    amplitude_increment = [0.001, -0.001]
#
#    for i in range(2):
#        if amplitude_start[i]+amplitude_increment[i]*n_wave>1 or amplitude_start[i]+amplitude_increment[i]*n_wave<0:
#            n_wave = 100
#            amplitude_start = [0.5, 0.5]
#            amplitude_increment = [0.001, -0.001]
#            print('Amplitude setting is out of range!')
#
#    awg_program = f'''
#        const samples = {pulse_length};
#        wave w = gauss(samples, 0.5, samples/2, samples/8); // 'w' can also be '.csv' file
#        wave w1 = gauss(samples, 0.5, samples/2, samples/8);
#        assignWaveIndex(w, w1, 0);        //For dual-channel waveforms
#        playWave(w,w1);
#        executeTableEntry(0);
#        repeat({n_wave}){{
#        executeTableEntry(1);
#        playZero({pulse_length}); // remain same cycle length
#        }}
#        '''
#    # compile seqc code
#    create_and_compile_awg(daq, device_awg, awg_program, 0, 1)
#    # generate json file
#    ct=ct_amplitude_increment(amplitude_start =amplitude_start, amplitude_increment = amplitude_increment)
#    # upload json file
#    daq.setVector(f"/{device_awg}/awgs/0/commandtable/data", json.dumps(ct))
#
#    # turn on wave output 1 and 2
#    for i in range(2):
#        daq.setInt(f'/{device_awg}/sigouts/{i}/on', 1)
#    # run AWG
##    daq.setInt(f'{device_awg}/awgs/0/enable', 1)
##    print('enable hd_awg')
#    time.sleep(20)
##    daq.setInt(f'{device_awg}/awgs/0/enable', 0)
#    print('disable hd_awg')
#    # get what json file have been loaded
#    json.loads(daq.get(f"/{device_awg}/awgs/0/commandtable/data",flat=True)[f"/{device_awg}/awgs/0/commandtable/data"][0]['vector'])
#    # currently the return value of amplitude setting is inverted, and this have been solved and will be available in the next release
#
    return ct
