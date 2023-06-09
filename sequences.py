# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 15:53:00 2023

@author: Evangelos Vlachos <evlachos@usc.edu>

"""

def trigger_readout_sequence():
    
    awg_program = '''playWave(1,w_marker);
            playZero(qubit_reset_time,AWG_RATE_1P2MHZ); '''
                     
    return awg_program

def tomographic_pulse_sequence(axis='Z'):
    
    if axis == 'X':
        awg_program = '''playWave(1,2,w_pi2_I,1,2,w_pi2_Q,AWG_RATE_2400MHZ);'''
    elif axis == 'Y':
        awg_program = '''playWave(1,2,w_pi2_Q,1,2,w_pi2_I,AWG_RATE_2400MHZ);'''
    elif axis == 'Z':
        awg_program = ''''''
        
    return awg_program

def finalize_sequence(awg_program,axis='Z'):
    
    awg_program = awg_program.replace('_trigger_readout_',trigger_readout_sequence())
    awg_program = awg_program.replace('_tomographic_pulse_',tomographic_pulse_sequence(axis))
    
    return awg_program

def spec_sequence(on_off=True):

    '''Generate qubit spectroscopy sequence'''
    
    if on_off:
        awg_program = '''       

        // Beginning of the core sequencer program executed on the HDAWG at run time
        wave w_marker = 2*marker(512,1);

        repeat(n_avg) {
            // OFF Measurement
            _trigger_readout_
            // ON measurement
            playWave(1,w_const,2,w_zero);
            _trigger_readout_
                    }'''
    else:
        awg_program = '''       

        // Beginning of the core sequencer program executed on the HDAWG at run time
        wave w_marker = 2*marker(512,1);

        repeat(n_avg) {
            playWave(1,w_const,2,w_zero);
            _trigger_readout_
                    }'''
    
    
    return awg_program

def rabi_sequence(tomography=False):
    '''Generate qubit spectroscopy sequence'''
   
    awg_program = '''
    resetOscPhase();
   
    wave w_marker = 2*marker(512,1);
    var i;
    // Beginning of the core sequencer program executed on the HDAWG at run time
    repeat(n_avg){
        for (i=0; i<n_steps; i++) {
                playWave(1,2,w_gauss_rise_I,1,2,w_zero1);
                executeTableEntry(i);
                playWave(1,2,w_gauss_fall_I,1,2,w_zero2);
                _tomographic_pulse_
                _trigger_readout_
      }
    }'''

    awg_program = awg_program.replace('_trigger_readout_',trigger_readout_sequence())
    
    return awg_program
 
def power_rabi_sequence(tomography=False):
    '''Generate qubit spectroscopy sequence'''
   
    awg_program = '''
    resetOscPhase();
   
    wave w_marker = 2*marker(512,1);
    var i;
    // Beginning of the core sequencer program executed on the HDAWG at run time
    repeat(n_avg){
        for (i=0; i<n_steps; i++) {
                executeTableEntry(i);
                _tomographic_pulse_
                _trigger_readout_
      }
    }'''

    awg_program = awg_program.replace('_trigger_readout_',trigger_readout_sequence())
    
    return awg_program
 
def T1_sequence(tomography=False):
    '''Generate qubit spectroscopy sequence'''
   
    awg_program = '''
    resetOscPhase();
    wave w_marker = 2*marker(512,1);
    var i;
    // Beginning of the core sequencer program executed on the HDAWG at run time
    repeat(n_avg){
        for (i=0; i<n_steps; i++) {
                playWave(1,2,w_pi_I,1,2,w_pi_Q,AWG_RATE_2400MHZ);
                executeTableEntry(i);
                _tomographic_pulse_
                _trigger_readout_
      }
    }'''
    
    awg_program = awg_program.replace('_trigger_readout_',trigger_readout_sequence())
    
    return awg_program
 
def ramsey_sequence(tomography=False):
    '''Generate qubit spectroscopy sequence'''
   
    awg_program = '''
        resetOscPhase();
       
        wave w_marker = 2*marker(512,1);
        var i;
        // Beginning of the core sequencer program executed on the HDAWG at run time
        repeat(n_avg){
            for (i=0; i<n_steps; i++) {
                    playWave(1,2,w_pi2_I,1,2,w_pi2_Q,AWG_RATE_2400MHZ);
                    executeTableEntry(i);
                    playWave(1,2,w_pi2_I,1,2,w_pi2_Q,AWG_RATE_2400MHZ);
                    _tomographic_pulse_
                    _trigger_readout_
          }
        }'''

    
    return awg_program

def echo_sequence(tomography=False):

    awg_program = '''
    resetOscPhase();
    wave w_marker = 2*marker(512,1);
    var i;
    // Beginning of the core sequencer program executed on the HDAWG at run time
    repeat(n_avg){
        for (i=0; i<n_steps; i++) {
                playWave(1,2,w_pi2_I,1,2,w_pi2_Q,AWG_RATE_2400MHZ);
                executeTableEntry(i);
                playWave(1,2,w_pi_I,1,2,w_pi_Q,AWG_RATE_2400MHZ);
                executeTableEntry(i);
                playWave(1,2,w_pi2_I,1,2,w_pi2_Q,AWG_RATE_2400MHZ);
                _tomographic_pulse_
                _trigger_readout_
      }
    }'''

    return awg_program
 
def z_gate_sequence(tomography=False):
    '''Generate qubit spectroscopy sequence'''
   
    awg_program = '''
   
    resetOscPhase();
    wave w_marker = 2*marker(512,1);
    var i;
    // Beginning of the core sequencer program executed on the HDAWG at run time
    repeat(n_avg){
        for (i=0; i<n_steps; i++) {
                resetOscPhase();
                executeTableEntry(0);
                playWave(1,2,w_zero);
                executeTableEntry(i);
                _tomographic_pulse_
                _trigger_readout_
      }
    }'''

    return awg_program
 
def single_shot_sequence():
    '''Generate qubit spectroscopy sequence'''
   
    awg_program = '''
    resetOscPhase();
   
    wave w_marker = 2*marker(512,1);
    // Beginning of the core sequencer program executed on the HDAWG at run time
    repeat(num_samples) {
        // OFF Measurement
        playWave(1,w_marker);
        playZero(qubit_reset_time,AWG_RATE_1P2MHZ);
        playWave(1,2,w_pi_I,1,2,w_pi_Q,AWG_RATE_2400MHZ);
        _trigger_readout_
        }'''

    return awg_program
 
def mixer_calib_sequence():
    '''Generate qubit spectroscopy sequence'''
   
    awg_program = '''
        resetOscPhase();
        wave w_const = amp*ones(N);
        wave w_zeros = zeros(N);

        playWave(1,2,w_const,1,2,w_zeros);
        playHold(1e9);
    '''

    return awg_program
 
def reset_sequence():
    '''Generate qubit spectroscopy sequence'''
   
    awg_program = '''
        waitDigTrigger(1);
        wait(1);
        wait(2000);
        if (getDigTrigger(2) == 0) {
            wait(10);
        } else {
            playWave(1,2,w_pi_I,1,2,w_pi_Q,AWG_RATE_2400MHZ);
            }
        wait(10000);
      '''
     
    return awg_program