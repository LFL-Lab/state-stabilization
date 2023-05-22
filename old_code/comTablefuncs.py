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

def ct_pulse_length(n_wave=10, pulse_length_start = 32, pulse_length_increment = 16,
                    sequence='rabi',pipulse=50,active_reset=False,noise_rate=0):

    ct = {'header': {'version':'1.2'}, 'table':[]}
    # make table entries for waveform playback during free evolution
    if sequence == 'ramsey':
        pulse_length_start = 0
        for i in range(n_wave-1):
            wfm_length = pulse_length_start + (i+2) * pulse_length_increment
            ct = make_entry(ct,wfm_index=0,table_index=i,length=wfm_length)
    else:
        for i in range(n_wave):
            wfm_length = pulse_length_start + i * pulse_length_increment
            ct = make_entry(ct=ct,wfm_index=0,table_index=i,length=wfm_length)

    # make pre and post pulses in the case of AC stark noise added to the system
    if sequence != 'rabi':
        ct = make_entry(ct,wfm_index=1,table_index=n_wave+1,length=int(pipulse/2))

        if active_reset == True:
            ct = make_entry(ct,wfm_index=1,table_index=n_wave+2,length=int(pipulse))

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


def make_entry(ct,wfm_index=0,table_index=0,length=16,amplitude_ch1=1,amplitude_ch2=1):

    entry = {'index': table_index,
              'waveform':{
                  'index': wfm_index,
                  'length':  length,
                  # 'samplingRateDivider': noise_rate,
              },
            }
    ct['table'].append(entry)

    return ct

def make_schema():

    schema = {

    "title": "AWG Command Table Schema",
    "description": "Schema for ZI HDAWG AWG Command Table",
    "definitions": {
        "header": {
            "properties": {
                "version": {
                    "type": "string",
                    "enum": [
                        "1.2"
                    ],
                    "description": "File format version. This version must match with the relevant schema version."
                },
                "partial": {
                    "description": "Set to True for incremental table updates",
                    "type": "boolean",
                    "default": "False"
                },
                "userString": {
                    "description": "User-definable label",
                    "type": "string",
                    "maxLength": 30
                }
            },
            "required": [
                "version"
            ]
        },
        "table": {
            "items": {
                "$ref": "#/definitions/entry"
            },
            "minItems": 0,
            "maxItems": 1024
        },
        "entry": {
            "properties": {
                "index": {
                    "$ref": "#/definitions/tableindex"
                },
                "waveform": {
                    "$ref": "#/definitions/waveform"
                },
                "phase0": {
                    "$ref": "#/definitions/phase"
                },
                "phase1": {
                    "$ref": "#/definitions/phase"
                },
                "amplitude0": {
                    "$ref": "#/definitions/amplitude"
                },
                "amplitude1": {
                    "$ref": "#/definitions/amplitude"
                }
            },
            "additionalProperties": False,
            "required": [
                "index"
            ]
        },
        "tableindex": {
            "type": "integer",
            "minimum": 0,
            "maximum": 1023,
            "exclusiveMinimum": False,
            "exclusiveMaximum": False
        },
        "waveform": {
            "properties": {
                "index": {
                    "$ref": "#/definitions/waveformindex"
                },
                "length": {
                    "$ref": "#/definitions/waveformlength"
                },
                "samplingRateDivider": {
                    "$ref": "#/definitions/samplingratedivider"
                },
                "awgChannel0": {
                    "$ref": "#/definitions/awgchannel"
                },
                "awgChannel1": {
                    "$ref": "#/definitions/awgchannel"
                },
                "precompClear": {
                    "$ref": "#/definitions/precompclear"
                },
                "playZero": {
                    "$ref": "#/definitions/playzero"
                }
            },
            "additionalProperties": False,
            "oneOf": [
                {
                    "required": [
                        "index"
                    ]
                },
                {
                    "required": [
                        "playZero",
                        "length"
                    ]
                }
            ]
        },
        "waveformindex": {
            "description": "Index of the waveform to play as defined with the assignWaveIndex sequencer instruction",
            "type": "integer",
            "minimum": 0,
            "maximum": 65535,
            "exclusiveMinimum": False,
            "exclusiveMaximum": False
        },
        "waveformlength": {
            "description": "The length of the waveform in samples",
            "type": "integer",
            "multipleOf": 16,
            "minimum": 32,
            "exclusiveMinimum": False
        },
        "samplingratedivider": {
            "descpription": "Integer exponent n of the sampling rate divider: 2.4 GSa/s / 2^n, n in range 0 ... 13",
            "type": "integer",
            "minimum": 0,
            "maximum": 13
        },
        "awgchannel": {
            "description": "Assign the given AWG channel to signal output 0 & 1",
            "type": "array",
            "minItems": 1,
            "maxItems": 2,
            "uniqueItems": True,
            "items": [
                {
                    "type": "string",
                    "enum": [
                        "sigout0",
                        "sigout1"
                    ]
                }
            ]
        },
        "precompclear": {
            "description": "Set to True to clear the precompensation filters",
            "type": "boolean",
            "default": False
        },
        "playzero": {
            "description": "Play a zero-valued waveform for specified length of waveform, equivalent to the playZero sequencer instruction",
            "type": "boolean",
            "default": "False"
        },
        "phase": {
            "properties": {
                "value": {
                    "description": "Phase value of the given sine generator in degree",
                    "type": "number"
                },
                "increment": {
                    "description": "Set to True for incremental phase value, or to False for absolute",
                    "type": "boolean",
                    "default": "False"
                }
            },
            "additionalProperties": False,
            "required": [
                "value"
            ]
        },
        "amplitude": {
            "properties": {
                "value": {
                    "description": "Amplitude scaling factor of the given AWG channel",
                    "type": "number",
                    "minimum": -1.0,
                    "maximum": 1.0,
                    "exclusiveMinimum": False,
                    "exclusiveMaximum": False
                },
                "increment": {
                    "description": "Set to True for incremental amplitude value, or to False for absolute",
                    "type": "boolean",
                    "default": "False"
                }
            },
            "additionalProperties": False,
            "required": [
                "value"
            ]
        }
    },
    "properties": {
        "$schema": {
            "type": "string"
        },
        "header": {
            "$ref": "#/definitions/header"
        },
        "table": {
            "$ref": "#/definitions/table"
        }
    },
    "additionalProperties": False,
    "required": [
        "header"
    ]
        }

    return schema