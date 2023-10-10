# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 15:53:00 2023

@authors: Evangelos Vlachos <evlachos@usc.edu>
        Malida Hecht <mohecht@usc.edu>
"""

import numpy as np
from scipy.signal.windows import gaussian
from zhinst.toolkit.waveform import Wave, OutputType
from zhinst.toolkit import CommandTable,Sequence,Waveforms
from utilities import roundToBase,gen_arb_wfm
import matplotlib.pyplot as plt


#%% sequence_code_funcs
def gen_seq_code(exp_pars):

    exp = exp_pars['exp']
    try:
        axis = exp_pars['tomographic-axis']
    except:
        pass
    try:
        on_off = exp_pars['on_off']
    except:
        pass
    try:
        active_reset = exp_pars['active_reset']
    except:
        active_reset = False
    if exp == 'coherence-stabilization':
        axis = 'Z'
    
    if exp == 'spectroscopy':
        code = spec_sequence(on_off)
    elif exp == 't-rabi':
        code = time_rabi_sequence()
    elif exp == 'p-rabi':
        code = power_rabi_sequence()
    elif exp == 'calibrate-rabi':
        code = calibrate_rabi_sequence()
    elif exp == 'T1':
        code = T1_sequence()
    elif exp == 'ramsey':
        code = ramsey_sequence()
    elif exp == 'echo':
        code = echo_sequence()
    elif exp == 'tomography':
        code = tomography_sequence()
    elif exp == 'tomography-calibration':
        code = tomography_calibration_sequence()
    elif exp == 'coherence-stabilization':
        code = state_stabilization_sequence()
    elif exp == 'z-gate':
        code = z_gate_sequence()
    elif exp == 'single-shot':
        code = single_shot_sequence()
    elif exp == 'mixer-calibration':
        code = mixer_calib_sequence()
    
    code = finalize_sequence(code,active_reset,axis)
    
    return code  
  
# def prepare_state(awg_program, state='X'):
    
#     # if state == 'X':
#     #     code = '''playWave(1,2,w_prep_I,1,2,w_prep_Q,AWG_RATE_2400MHZ);'''
#     # elif state == 'Y':
#     #     code = '''playWave(1,2,w_prep_Q,1,2,w_prep_Q,AWG_RATE_2400MHZ);'''
#     # elif state == '0':
#     #     code = ''''''
#     # elif state == '1':
#     #     code = '''playWave(1,2,w_pi_I,1,2,w_pi_Q,AWG_RATE_2400MHZ);'''
        
#     return code
 
def trigger_readout_sequence():
    
    awg_program = '''playWave(1,w_marker);'''
                     
    return awg_program

def tomographic_pulse_sequence(axis='Z'):
    
    if axis == 'X' or axis == 'Y':
        awg_program = '''executeTableEntry(n_steps);'''
    elif axis == 'Z':
        awg_program = ''''''
        
    return awg_program

def qubit_reset_code(active_reset):
    
    if active_reset:
        awg_program = qubit_reset_sequence()
    else:
        awg_program =  '''playZero(qubit_reset_time,AWG_RATE_1P2MHZ);'''
    
    return awg_program

def finalize_sequence(awg_program,active_reset=False,axis='Z'):
    
    awg_program = awg_program.replace('_trigger_readout_',trigger_readout_sequence())
    awg_program = awg_program.replace('_tomography_pulse_',tomographic_pulse_sequence(axis))
    awg_program = awg_program.replace('_qubit_reset_',qubit_reset_code(active_reset))
    
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
            wait(5000);
            // ON measurement
            playWave(1,w_const,2,w_zero);
            _trigger_readout_
            _qubit_reset_
                    }'''
    else:
        awg_program = '''       

        // Beginning of the core sequencer program executed on the HDAWG at run time
        wave w_marker = 2*marker(512,1);

        repeat(n_avg) {
            playWave(1,w_const,2,w_zero);
            _trigger_readout_
            _qubit_reset_
                    }'''
    
    
    return awg_program

def time_rabi_sequence():
    '''Generate qubit spectroscopy sequence'''
   
    awg_program = '''
    resetOscPhase();
   
    wave w_marker = 2*marker(512,1);
    var i;
    // Beginning of the core sequencer program executed on the HDAWG at run time
    repeat(n_avg){
        for (i=0; i<n_steps+1; i++) {
                if (i==0) {
                        }
                else {
                executeTableEntry(n_steps+1);
                executeTableEntry(i-1);
                executeTableEntry(n_steps+2);
                }
                _tomography_pulse_
                _trigger_readout_
                _qubit_reset_
      }
    }'''

    awg_program = awg_program.replace('_trigger_readout_',trigger_readout_sequence())
    
    return awg_program

def tomography_calibration_sequence():
    '''Generate qubit spectroscopy sequence'''
   
    awg_program = '''
    resetOscPhase();
   
    wave w_marker = 2*marker(512,1);
    var i;
    // Beginning of the core sequencer program executed on the HDAWG at run time
    repeat(n_avg){
        for (i=0; i<n_steps+1; i++) {
                
        // x tomography
                if (i==0) {
                        }
                else {
                executeTableEntry(n_steps+1);
                executeTableEntry(i-1);
                executeTableEntry(n_steps+2);
                }
                executeTableEntry(n_steps+3);
                _trigger_readout_
                _qubit_reset_

        // y tomography
        
                if (i==0) {
                        }
                else {
                executeTableEntry(n_steps+1);
                executeTableEntry(i-1);
                executeTableEntry(n_steps+2);
                }
                executeTableEntry(n_steps+4);
                _trigger_readout_
                _qubit_reset_
      
        // z tomography
        
                if (i==0) {
                        }
                else {
                executeTableEntry(n_steps+1);
                executeTableEntry(i-1);
                executeTableEntry(n_steps+2);
                }
                executeTableEntry(n_steps+5);
                _trigger_readout_
                _qubit_reset_
      }  
      
    }'''

    # awg_program = awg_program.replace('_trigger_readout_',trigger_readout_sequence())
    
    return awg_program
 
def power_rabi_sequence():
    '''Generate qubit spectroscopy sequence'''
   
    awg_program = '''
    resetOscPhase();
   
    wave w_marker = 2*marker(512,1);
    var i;
    // Beginning of the core sequencer program executed on the HDAWG at run time
    repeat(n_avg){
        for (i=0; i<n_steps; i++) {
                executeTableEntry(i);
                _tomography_pulse_
                _trigger_readout_
                _qubit_reset_
      }
    }'''

    awg_program = awg_program.replace('_trigger_readout_',trigger_readout_sequence())
    
    return awg_program
 
def calibrate_rabi_sequence():
    '''Generate qubit spectroscopy sequence'''
    
    awg_program = '''
    resetOscPhase();
    
    wave w_marker = 2*marker(512,1);
    var i;
    // Beginning of the core sequencer program executed on the HDAWG at run time
    repeat(n_avg){
        for (i=0; i<n_steps; i++) {
                 repeat (i){
                    executeTableEntry(0);
                }
                
                _tomography_pulse_
                _trigger_readout_
                _qubit_reset_
       }
     }'''

    awg_program = awg_program.replace('_trigger_readout_',trigger_readout_sequence())
     
    return awg_program   
 
def T1_sequence():
    '''Generate qubit spectroscopy sequence'''
   
    awg_program = '''
    resetOscPhase();
    wave w_marker = 2*marker(512,1);
    var i;
    // Beginning of the core sequencer program executed on the HDAWG at run time
    repeat(n_avg){
        for (i=0; i<n_steps; i++) {
                executeTableEntry(n_steps);
                executeTableEntry(i);
                _tomography_pulse_
                _trigger_readout_
                _qubit_reset_
      }
    }'''
    
    awg_program = awg_program.replace('_trigger_readout_',trigger_readout_sequence())
    
    return awg_program
 
def ramsey_sequence():
    '''Generate qubit spectroscopy sequence'''
   
    awg_program = '''
        resetOscPhase();
       
        wave w_marker = 2*marker(512,1);
        var i;
        // Beginning of the core sequencer program executed on the HDAWG at run time
        repeat(n_avg){
            for (i=0; i<n_steps; i++) {
                    executeTableEntry(n_steps+1);
                    executeTableEntry(i);
                    executeTableEntry(n_steps+1);
                    _tomography_pulse_
                    _trigger_readout_
                    _qubit_reset_
          }
        }'''

    
    return awg_program

def echo_sequence():

    awg_program = '''
    resetOscPhase();
    wave w_marker = 2*marker(512,1);
    var i;
    // Beginning of the core sequencer program executed on the HDAWG at run time
    repeat(n_avg){
        for (i=0; i<n_steps; i++) {
                executeTableEntry(n_steps+1);// pi /2 pulse
                executeTableEntry(i); //zero pulse for t/2 time
                executeTableEntry(n_steps+2); // pi pulse
                executeTableEntry(i); // zero for t/2
                executeTableEntry(n_steps+1); // pi/2
                _trigger_readout_
                _qubit_reset_
      }
    }'''

    return awg_program
 
def z_gate_sequence():
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
                _trigger_readout_
                _qubit_reset_
      }
    }'''

    return awg_program

def tomography_sequence():
    
    awg_program = '''
   
    resetOscPhase();
    wave w_marker = 2*marker(512,1);
    var i;
    // Beginning of the core sequencer program executed on the HDAWG at run time
    repeat(n_avg){
        for (i=0; i<n_steps; i++) {
                _prepare_state_
               executeTableEntry(i);
                _trigger_readout_
                _qubit_reset_
      }
    }'''
    
    

    return awg_program



def state_stabilization_sequence():
    
    awg_program = '''
   
    resetOscPhase();
    wave w_marker = 2*marker(512,1);
    var i;
    // Beginning of the core sequencer program executed on the HDAWG at run time
    repeat(n_avg){
        wait(100000);
        _trigger_readout_ // gnd state measurement for calibration
        _qubit_reset_
        for (i=0; i<n_steps; i++) {
                executeTableEntry(n_steps); // prep pulse
                wait(50);
                executeTableEntry(i); // evolution
                wait(50);
                executeTableEntry(n_steps+1); // tomography x
                _trigger_readout_
                _qubit_reset_
                executeTableEntry(n_steps); //preparation pulse
                wait(50);
                executeTableEntry(i);
                wait(50);
                executeTableEntry(n_steps+2); // tomography y
                _trigger_readout_
                _qubit_reset_
                executeTableEntry(n_steps); // preparation pulse
                wait(50);
                executeTableEntry(i); // evolution
                wait(50);
                executeTableEntry(n_steps+3); // tomography z
                _trigger_readout_
                _qubit_reset_
                }
                executeTableEntry(n_steps+4); // pi pulse
                _trigger_readout_ // exc state meas for calib
                _qubit_reset_  
      }
    '''
    return awg_program

# def state_stabilization_sequence():
    
#     awg_program = '''
   
#     resetOscPhase();
#     wave w_marker = 2*marker(512,1);
#     var i;
#     // Beginning of the core sequencer program executed on the HDAWG at run time
#     repeat(n_avg){
#         wait(100000);
#         _trigger_readout_ // gnd state measurement for calibration
#         _qubit_reset_
#         for (i=0; i<n_steps; i++) {
#                 wait(100);
#                 executeTableEntry(n_steps);
#                 wait(50);
#                 executeTableEntry(i);
#                 wait(50);
#                 executeTableEntry(n_steps+1);
#                 _trigger_readout_
#                 _qubit_reset_
#                 }
#                 executeTableEntry(n_steps+2); //do pi pulse
#                 _trigger_readout_ // exc state meas for calib
#                 _qubit_reset_  
#       }
#     '''
    
#     return awg_program

def single_shot_sequence():
    '''Generate qubit spectroscopy sequence'''
   
    awg_program = '''
    resetOscPhase();
   
    wave w_marker = 2*marker(512,1);
    // Beginning of the core sequencer program executed on the HDAWG at run time
    repeat(512) {
        // OFF Measurement
        _trigger_readout_
        _qubit_reset_
        executeTableEntry(0);
        _trigger_readout_
        _qubit_reset_
        }'''

    return awg_program
 
def mixer_calib_sequence():
    '''Generate qubit spectroscopy sequence'''
   
    awg_program = '''
        resetOscPhase();
        
        while (true) {
        executeTableEntry(0);
        }
        //playWave(1,2,w_gauss,1,2,w_zeros);
        //playHold(1e9);
    '''

    return awg_program
 
def qubit_reset_sequence():
    '''Generate qubit spectroscopy sequence'''
   
    awg_program = '''
        waitDigTrigger(1);
        wait(1000);
        if (getDigTrigger(2) == 0) {
            wait(300);
        } else {
            executeTableEntry(1023);
            wait(20000);
            }
        wait(20000);
      '''
     
    return awg_program
#executeTableEntry(1023);
#%% waveform_funcs
def make_wave(pulse_type='gauss', wave_name='wave', wfm_pars = {}, exp_pars ={},qb_pars={},amplitude = 1.0, pulse_length = 16 , output_order='1' ):
    '''
    Generates string to be used for generating waveforms to be uploaded to the AWG 
    
    Input:             Type:            Description:
    pulse_type:        string           type of pulse: Gaussian, Constant, Pi, Pi/2, Rise and 
                                        falling Gaussian.
                                        'zero' or 'constant' will generate constant pulses,
                                        'rise' and 'fall' create a rising or falling half gaussian
                                        all other options will generate a gaussian type pulse
    wave_name:         string          name that will be uploaded to awg to define the waveform
    amplitude:         float           amplitude of waveform, will be defined 
    pulse_length:      int             pulse length in number of samples; for
                                        Zurich Instruments, it's in base_theta 16. 
    output_order:      float


    Output:           Type:             Description:
    wave:             string            wave string description to be uploaded to Zurich
    '''
    # (╯°□°）╯︵ ┻━┻ 
    # generate pulse type
    if pulse_type == 'zero':
        pulse = np.zeros(pulse_length)
    elif pulse_type == 'constant':
        pulse = amplitude * np.ones(pulse_length)
    elif pulse_type == 'rise':
        gauss_pulse = amplitude * gaussian(pulse_length , pulse_length/5)
        pulse = gauss_pulse[:int(pulse_length/2)]
    elif pulse_type == 'fall':
        gauss_pulse = amplitude * gaussian(pulse_length, pulse_length/5)
        pulse = gauss_pulse[int(pulse_length/2):]
    elif pulse_type == 'arb_I':
        #pulse = amplitude*gen_arb_wfm('rising',wfm_pars,normalize=True)
        pulse = amplitude *  gen_arb_wfm('rising',
                                         wfm_pars=wfm_pars,
                                         exp_pars= exp_pars,
                                         qb_pars=qb_pars,
                                         normalize=True)
    elif pulse_type == 'arb_Q':
        #pulse = amplitude*gen_arb_wfm('markov',wfm_pars,n_points=pulse_length)
        pulse = amplitude* gen_arb_wfm('markov',
                                       wfm_pars=wfm_pars,
                                       exp_pars=exp_pars,
                                       qb_pars=qb_pars,
                                       n_points=pulse_length)
    else:
        # This should work for the rest of them ಠ_ಠ
        pulse = amplitude * gaussian(pulse_length , pulse_length/5)

    # output parsing 
    if output_order == '1' :
        out = OutputType.OUT1
    elif output_order == '2':
        out = OutputType.OUT2
    elif output_order == '12':
        out = OutputType.OUT1|OutputType.OUT2
    elif output_order == '21':
        out = OutputType.OUT2|OutputType.OUT1
    elif output_order == '11':
        out = OutputType.OUT1|OutputType.OUT1
    elif output_order == '22':
        out = OutputType.OUT2|OutputType.OUT2

    return pulse, wave_name, out

def setup_waveforms(sequence,wfm_pars={},exp_pars={},qb_pars={},n_points=1024):
    ''' Creates the waveforms Necessary for Each Experiment'''
    exp = exp_pars['exp']
    sequence.waveforms = Waveforms()

    if exp == 'spectroscopy':
        qubit_drive_dur = roundToBase(exp_pars['satur_dur']*exp_pars['fsAWG'])
        amp = 2 * exp_pars['amp_q']
        sequence.waveforms[0] = (
            Wave(*make_wave(pulse_type = 'constant',
                        wave_name = 'w_const',
                        amplitude = amp,
                        pulse_length = qubit_drive_dur,
                        output_order = '1')), 
            Wave(*make_wave(pulse_type = 'zero',
                        wave_name = 'w_zero',
                        amplitude = amp,
                        pulse_length = qubit_drive_dur,
                        output_order = '2')) )

    elif exp == 't-rabi' or exp == 'tomography-calibration':
        N = 64
        amp = exp_pars['amp_q']       
    
        sequence.waveforms[0] = (Wave(*make_wave(pulse_type = 'rise',
                        wave_name = 'w_gauss_rise_I',
                        amplitude = amp,
                        pulse_length = N,
                        output_order = '12')), 
            Wave(*make_wave(pulse_type = 'zero',
                        wave_name = 'w_zero1',
                        amplitude = 0,
                        pulse_length = int(N/2), ## Not sure about this one... 
                        output_order = '12')))
    
        sequence.waveforms[1] = (
            Wave(*make_wave(pulse_type = 'fall',
                        wave_name = 'w_gauss_fall_I',
                        amplitude = amp,
                        pulse_length = N,
                        output_order = '12')), 
            Wave(*make_wave(pulse_type = 'zero',
                        wave_name = 'w_zero2',
                        amplitude = 0,
                        pulse_length = int(N/2), ## Not sure about this one... 
                        output_order = '12'))
    )
        sequence.waveforms[2] = (
        Wave(*make_wave(pulse_type = 'constant',
                        wave_name = 'w_const_I',
                        amplitude = amp,
                        pulse_length = n_points,
                        output_order = '12')), 
        Wave(*make_wave(pulse_type = 'zero',
                        wave_name = 'w_const_Q',
                        amplitude = 0,
                        pulse_length = n_points,
                        output_order = '12'))
    )
        
        sequence.waveforms[3] = (
        Wave(*make_wave(pulse_type = 'gaussian',
                        wave_name = 'w_gauss_I',
                        amplitude = 1,
                        pulse_length = qb_pars['pi_len'],
                        output_order = '12')), 
        Wave(*make_wave(pulse_type = 'zero',
                        wave_name = 'w_gauss_Q',
                        amplitude = 0,
                        pulse_length = qb_pars['pi_len'],
                        output_order = '12'))
        )
        
        

    elif exp == 'p-rabi':
        N = qb_pars['gauss_len']

        sequence.waveforms[0] = (
        Wave(*make_wave(pulse_type = 'gaussian',
                        wave_name = 'w_gauss',
                        amplitude = 1,
                        pulse_length = N,
                        output_order = '12')), 
        Wave(*make_wave(pulse_type = 'zero',
                        wave_name = 'w_zero',
                        amplitude = 0,
                        pulse_length = N,
                        output_order = '12'))
    )

    elif exp == 'T1':
        N = qb_pars['pi_len']
        amp = qb_pars['pi_amp']

        sequence.waveforms[0] = (
            Wave(*make_wave(pulse_type = 'pi',
                        wave_name = 'w_pi_I',
                        amplitude = amp,
                        pulse_length = N,
                        output_order = '12')), 
            Wave(*make_wave(pulse_type = 'zero',
                        wave_name = 'w_pi_Q',
                        amplitude = 0,
                        pulse_length = N,
                        output_order = '12'))
    )

        sequence.waveforms[1] = (
            Wave(*make_wave(pulse_type = 'zero',
                        wave_name = 'w_zero_I',
                        amplitude = 0,
                        pulse_length =n_points,
                        wfm_pars=wfm_pars,
                        output_order = '12')), 
            Wave(*make_wave(pulse_type = 'arb_Q',
                        wave_name = 'w_zero_Q',
                        amplitude = 1,
                        pulse_length = n_points,
                        wfm_pars=wfm_pars,
                        output_order = '12'))
    ) 

    elif exp == 'ramsey':
        N = qb_pars['pi_len']
        amp =  qb_pars['pi_half_amp']
        
        sequence.waveforms[0] = (
        Wave(*make_wave(pulse_type = 'pi2',
                        wave_name = 'w_pi2_I',
                        amplitude = amp,
                        pulse_length =N,
                        output_order = '12')), 
        Wave(*make_wave(pulse_type = 'zero',
                        wave_name = 'w_pi2_Q',
                        amplitude = 0,
                        pulse_length = N,
                        output_order = '12'))
    )    

        sequence.waveforms[1] = (
        Wave(*make_wave(pulse_type = 'zero',
                        wave_name = 'w_zero_I',
                        amplitude = 0,
                        pulse_length = n_points,
                        output_order = '12')), 
        Wave(*make_wave(pulse_type = 'arb_Q',
                        wave_name = 'w_zero_Q',
                        amplitude = 1,
                        pulse_length = n_points,
                        wfm_pars = wfm_pars,
                        output_order = '12'))
    )    

    elif exp == 'echo':
        N = qb_pars['pi_len']
        amp_half =  qb_pars['pi_half_amp']
        amp = qb_pars['pi_amp']

        sequence.waveforms[0] = (
            Wave(*make_wave(pulse_type = 'pi2',
                        wave_name = 'w_pi2_I',
                        amplitude = amp_half,
                        pulse_length = N,
                        output_order = '12')), 
            Wave(*make_wave(pulse_type = 'zero',
                        wave_name = 'w_pi2_Q',
                        amplitude = 0,
                        pulse_length = N,
                        output_order = '12'))
    )    
    
        sequence.waveforms[1] = (
        Wave(*make_wave(pulse_type = 'pi',
                        wave_name = 'w_pi_I',
                        amplitude = amp,
                        pulse_length = N,
                        output_order = '12')), 
        Wave(*make_wave(pulse_type = 'zero',
                        wave_name = 'w_pi_Q',
                        amplitude = 0,
                        pulse_length = N,
                        output_order = '12'))
    ) 
    
        sequence.waveforms[2] = (
            Wave(*make_wave(pulse_type = 'zero',
                        wave_name = 'w_zero_I',
                        amplitude = 0,
                        pulse_length = n_points,
                        output_order = '12')), 
            Wave(*make_wave(pulse_type = 'zero',
                        wave_name = 'w_zero_Q',
                        amplitude = 0,
                        pulse_length = n_points,
                        output_order = '12'))
    )

    elif exp == 'z-gate':
        N = qb_pars['pi_len']
        amp = qb_pars['pi_half_amp']
        sequence.waveforms[0] = (
        Wave(*make_wave(pulse_type = 'pi2',
                        wave_name = 'w_pi2_I',
                        amplitude = amp_half,
                        pulse_length = N,
                        output_order = '12')), 
        Wave(*make_wave(pulse_type = 'zero',
                        wave_name = 'w_pi2_Q',
                        amplitude = 0,
                        pulse_length = N,
                        output_order = '12'))
    )    
        sequence.waveforms[1] = (
        Wave(*make_wave(pulse_type = 'zero',
                        wave_name= 'w_zero',
                        amplitude = 0,
                        pulse_length = n_points,
                        output_order = '12'))
    )

    elif exp == 'tomography':
        N = qb_pars['pi_len']
        amp_half = qb_pars['pi_half_amp']
        amp = qb_pars['pi_amp']

        sequence.waveforms[0] = (
        Wave(*make_wave(pulse_type = 'pi2',
                        wave_name = 'w_pi2x_I',
                        amplitude = amp_half,
                        pulse_length = N,
                        output_order = '12')), 
        Wave(*make_wave(pulse_type = 'zero',
                        wave_name = 'w_pi2x_Q',
                        amplitude = 0,
                        pulse_length = N,
                        output_order = '12'))
    )    

        sequence.waveforms[1] = (
        Wave(*make_wave(pulse_type = 'zero',
                        wave_name = 'w_zero_I',
                        amplitude = 0,
                        pulse_length = n_points,
                        output_order = '12')), 
        Wave(*make_wave(pulse_type = 'zero',
                        wave_name = 'w_zero_Q',
                        amplitude = 0,
                        pulse_length = n_points,
                        output_order = '12'))
    )

        sequence.waveforms[2] = (
        Wave(*make_wave(pulse_type = 'zero',
                        wave_name = 'w_pi2y_I',
                        amplitude = 0,
                        pulse_length = N,
                        output_order = '21')), 
        Wave(*make_wave(pulse_type = 'pi2',
                        wave_name = 'w_pi2y_Q',
                        amplitude = amp_half,
                        pulse_length = N,
                        output_order = '21'))
    )    

        sequence.waveforms[3] = (
        Wave(*make_wave(pulse_type = 'pi',
                        wave_name = 'w_pi_I',
                        amplitude = amp,
                        pulse_length = N,
                        output_order = '12')), 
        Wave(*make_wave(pulse_type = 'zero',
                        wave_name = 'w_pi_Q',
                        amplitude = 0,
                        pulse_length = N,
                        output_order = '12'))
        ) 
        
    elif exp == 'coherence-stabilization':
         N = qb_pars['pi_len']
         #amp_half = qb_pars['pi_half_amp']
         #pi_amp = qb_pars['pi_amp']
         # amp = wfm_pars['amp']
         
         sequence.waveforms[0] = (
         Wave(*make_wave(pulse_type = 'pi',
                         wave_name = 'w_prep_I',
                         amplitude = 1,
                         pulse_length = N,
                         output_order = '12')), 
         Wave(*make_wave(pulse_type = 'zero',
                         wave_name = 'w_prep_Q',
                         amplitude = 0,
                         pulse_length = N,
                         output_order = '12'))
         )    



         sequence.waveforms[1] = (
         Wave(*make_wave(pulse_type = 'arb_I',
                         wave_name = 'w_arb_I',
                         amplitude = 1,
                         pulse_length = n_points,
                         wfm_pars = wfm_pars,
                         exp_pars = exp_pars,
                         qb_pars = qb_pars,
                         output_order = '12')), 
         Wave(*make_wave(pulse_type = 'arb_Q',
                         wave_name = 'w_arb_Q',
                         amplitude = 1,
                         pulse_length = n_points,
                         wfm_pars=wfm_pars,
                         output_order = '12',
                         exp_pars = exp_pars,
                         qb_pars = qb_pars))
     )
         sequence.waveforms[2] = (
         Wave(*make_wave(pulse_type = 'pi',
                         wave_name = 'w_tom_I',
                         amplitude = 1,
                         pulse_length = N,
                         output_order = '12')), 
         Wave(*make_wave(pulse_type = 'zero',
                         wave_name = 'w_tom_Q',
                         amplitude = 0,
                         pulse_length = N,
                         output_order = '12')))
         
    elif exp == 'calibrate-rabi':
        N = qb_pars['pi_len']
        
        sequence.waveforms[0] = (
        Wave(*make_wave(pulse_type = 'pi',
                        wave_name = 'w_I',
                        amplitude = 1,
                        pulse_length = N,
                        output_order = '12')), 
        Wave(*make_wave(pulse_type = 'zero',
                        wave_name = 'w_Q',
                        amplitude = 0,
                        pulse_length = N,
                        output_order = '12')))
         
    elif exp == 'single-shot':
        N = qb_pars['pi_len']
        amp = qb_pars['pi_amp']
        sequence.waveforms[0] = (
        Wave(*make_wave(pulse_type = 'pi',
                        wave_name = 'w_pi_I',
                        amplitude = amp,
                        pulse_length = N,
                        output_order = '12')), 
        Wave(*make_wave(pulse_type = 'zero',
                        wave_name = 'w_pi_Q',
                        amplitude = 0,
                        pulse_length = N,
                        output_order = '12'))) 
        
    elif exp == 'mixer-calibration':
        N = 64
        amp = 1
        sequence.waveforms[0] = (
        Wave(*make_wave(pulse_type = 'constant',
                        wave_name = 'w_gauss1',
                        amplitude = amp,
                        pulse_length = N,
                        output_order = '12')), 
        Wave(*make_wave(pulse_type = 'constant',
                        wave_name = 'w_gauss2',
                        amplitude = amp,
                        pulse_length = N,
                        output_order = '12')))   

    if exp_pars['active_reset'] == True:
        N = qb_pars['pi_len']
        amp = qb_pars['pi_amp']
        sequence.waveforms[5] = (
        Wave(*make_wave(pulse_type = 'gaussian',
                        wave_name = 'w_reset_I',
                        amplitude = amp,
                        pulse_length = N,
                        output_order = '12')), 
        Wave(*make_wave(pulse_type = 'zero',
                        wave_name = 'w_reset_Q',
                        amplitude = 0,
                        pulse_length = N,
                        output_order = '12')))
    else:
        pass
        
    return sequence


def setup_seq_pars(sequence,exp,exp_pars={},qb_pars={},n_steps=100):
    
    # i = 0
    # for key,value in exp_pars.items():
    #     if (i > 1 and exp == 'spectroscopy') or i > 2:
    #         break
    #     else:
    #         if key == 'qubit_reset_time':
    #             value = roundTobase_theta(value*exp_pars['fsAWG']) # converts the reset time from us to num of samples
    #         else:
    #             pass
    #         sequence.constants[key] = value
    #     i += 1
    if exp_pars['exp'] != 'mixer-calibration':
        sequence.constants['n_avg'] = exp_pars['n_avg']
    if 'element' in exp_pars and exp_pars['element'] == 'rr':
        sequence.constants['qubit_reset_time'] = roundToBase(qb_pars['rr_reset_time']*1.17e6)
    else:
        if exp_pars['active_reset'] and exp_pars['exp'] != 'mixer-calibration':
                pass
        else:            
            sequence.constants['qubit_reset_time'] = roundToBase(qb_pars['qubit_reset_time']*1.17e6)
    if exp != 'spectroscopy':
        sequence.constants['n_steps'] = n_steps
        

#%% command_table_funcs
# Please refer to https://docs.zhinst.com/hdawg/commandtable/v2/schema for other settings
def make_ct(hdawg_core,exp_pars={},qb_pars={},wfm_pars={},x0=0,dx=16,n_steps=100):
   
    exp = exp_pars['exp']
    #base_theta = 180 + qb_pars['qb_mixer_imbalance'][1] 
    base_theta = -90
    
    if exp == 'p-rabi':
        sweep_var = 'amp'
    elif exp == 'z-gate':
        sweep_var = 'phase'
    elif exp == 'single-shot' or exp == 'coherence-stabilization' or exp == 'tomography-calibration' or exp == 'mixer-calibration':
        sweep_var = None
    elif exp == 'calibrate-rabi':
        sweep_var = None 
    elif exp != 'spectroscopy':
        sweep_var='time'
 
    
    
    ct = init_ct(hdawg_core)
    
    if exp == 't-rabi' or exp == 'echo':
        wfm_index = 2
    elif exp == 'p-rabi'  or exp == 'z-gate':
        wfm_index = 0
    else:
        wfm_index = 1
        
    if sweep_var == 'time':
        ct_sweep_length(ct,exp,wfm_index,qb_pars,x0,dx,0,n_steps)
    elif sweep_var == 'amp':
        ct_sweep_amp(ct,wfm_index,exp_pars,qb_pars,n_steps)
    elif sweep_var == 'phase':
        ct_sweep_phase(ct,exp_pars,wfm_index,n_steps)
    else:
        pass
        
    # if exp == 'coherence-stabilization':
    #     wfm_index = 0
    #     theta_prep,theta_tom,amp_prep,amp_tom = determine_axes_new(exp_pars,qb_pars)
    #     ct_sweep_length(ct,exp,1,qb_pars,x0,dx,90,n_steps) #Noise and coherent control wfms
    #     ct = arb_pulse(ct,n_steps,0,amp_prep,base_theta,theta_prep) # preparation pulse
    #     ct = arb_pulse(ct,n_steps+1,2,amp_tom,base_theta,theta_tom) # tomography pulse
    #     ct = arb_pulse(ct,n_steps+2,2,qb_pars['pi_amp'],base_theta,0) # pi pulse
    if exp == 'coherence-stabilization':
        wfm_index = 0
        theta_prep,amp_prep,base_amp = determine_axes_coherence(exp_pars,qb_pars)
        ct_sweep_length(ct,exp,1,qb_pars,x0,dx,90,n_steps) #Noise and coherent control wfms
        ct = arb_pulse(ct,n_steps,0,amp_prep,base_theta,theta_prep) # preparation pulse
        ct = arb_pulse(ct,n_steps+1,2,0.5*base_amp,base_theta,180) # tomography x pulse
        ct = arb_pulse(ct,n_steps+2,2,0.5*base_amp,base_theta,-90) # tomography y pulse
        ct = arb_pulse(ct,n_steps+3,2,0,base_theta,0) # tomography z pulse
        ct = arb_pulse(ct,n_steps+4,2,qb_pars['pi_amp'],base_theta,0) # pi pulse    
    elif exp == 't-rabi':
        wfm_index = 3
        theta_prep,theta_tom,amp_prep,amp_tom = determine_axes(exp_pars,qb_pars)
        ct = arb_pulse(ct,n_steps,wfm_index,amp_tom,base_theta,theta_tom) # tomography pulse
        wfm_index = 0 # gauss rise
        ct = arb_pulse(ct,n_steps+1,wfm_index,1,base_theta,0)
        wfm_index = 1 # gauss fall
        ct = arb_pulse(ct,n_steps+2,wfm_index,1,base_theta,0)
    elif exp == 'tomography-calibration':
        wfm_index = 0 # gauss rise
        ct = arb_pulse(ct,n_steps+1,wfm_index,1,base_theta,0)
        wfm_index = 1 # gauss fall
        ct = arb_pulse(ct,n_steps+2,wfm_index,1,base_theta,0)
        wfm_index = 2 # middle pulse
        ct_sweep_length(ct,exp,wfm_index,qb_pars,x0,dx,0,n_steps)
        wfm_index = 3
        exp_pars['tomographic-axis'] = 'X'
        theta_prep,theta_tom,amp_prep,amp_tom = determine_axes_new(exp_pars,qb_pars)
        ct = arb_pulse(ct,n_steps+3,wfm_index,amp_tom,base_theta,theta_tom) # X tomography pulse
        exp_pars['tomographic-axis'] = 'Y'
        theta_prep,theta_tom,amp_prep,amp_tom = determine_axes_new(exp_pars,qb_pars)
        ct = arb_pulse(ct,n_steps+4,wfm_index,amp_tom,base_theta,theta_tom) # Y tomography pulse
        ct = arb_pulse(ct,n_steps+5,wfm_index,0,base_theta,0) # Z tomography pulse
    elif exp == 'T1':
        wfm_index = 0
        ct = arb_pulse(ct,n_steps,wfm_index,1,base_theta,0) # pi_pulse
    elif exp == 'single-shot':
        wfm_index = 0
        ct = arb_pulse(ct,0,wfm_index,1,base_theta,0) # pi_pulse
    elif exp == 'ramsey':
        wfm_index = 0
        ct = arb_pulse(ct,n_steps+1,wfm_index,1,base_theta,0) # pi_pulse
    elif exp == 'mixer-calibration':
        wfm_index = 0
        ct = arb_pulse(ct,0,0,exp_pars['amp'],base_theta,90) 
        #ct = arb_pulse(ct,1,1,exp_pars['amp'],90,0) 
    elif exp == 'calibrate-rabi':
        wfm_index = 0
        ct = arb_pulse(ct,0,wfm_index,qb_pars['pi_half_amp']/4,base_theta,0) # pi pulse
    elif exp == 'echo':
        ct = arb_pulse(ct,n_steps+1,0,1,base_theta,0) # pi/2 pulse
        ct = arb_pulse(ct,n_steps+2,1,1,base_theta,0) #pi pulse
    if exp_pars['active_reset']:
        ct = arb_pulse(ct,1023,5,1,base_theta,0)
    
        
    else:
        pass
    
    return ct

    
def ct_sweep_length(ct,exp,wfm_index,qb_pars,x0,dx,theta,n_steps=100):
    
    # if exp == 't-rabi' or exp == 'tomography-calibration':
    #     n_steps += 1
    # elif exp == 'coherence-stabilization':
    #     n_steps += 2
    # else:
    #     pass
    
    for i in range(n_steps):
        wfm_length = x0 + i * dx
        # print(wfm_length)
        ct.table[i].waveform.index = wfm_index
        ct.table[i].waveform.length = wfm_length
        ct.table[i].amplitude0.value = 1.0
        ct.table[i].amplitude1.value = 1.0
        ct.table[i].phase0.value = theta
        #ct.table[i].phase1.value = 180 + qb_pars['qb_mixer_imbalance'][1] + theta
        ct.table[i].phase1.value = -90 + theta
        
def ct_sweep_amp(ct,wfm_index,exp_pars={},qb_pars={},n_steps=100):
    
    for i in range(n_steps):
        amp = exp_pars['x0'] + i * exp_pars['dx']
        ct.table[i].waveform.index = wfm_index
        ct.table[i].amplitude0.value = amp
        ct.table[i].amplitude1.value = amp
        ct.table[i].phase0.value = 90
        #ct.table[i].phase1.value =180 + qb_pars['qb_mixer_imbalance'][1]
        ct.table[i].phase1.value =0
# def ct_arb_amp(ct,wfm_index,exp_pars={},qb_pars9wfm_pars={},n_steps=100):
    
#     time_arr = np.arange(wfm_pars['t0'],wfm_pars['tmax']-wfm_pars['dt']/2,wfm_pars['dt'])
#     fun = lambda x : (1/np.sqrt(wfm_pars['tb'] - x))
#     amp = fun(time_arr)
#     plt.plot(time_arr,amp)
#     dx = int(len(time_arr)/n_steps)
#     for i in range(n_steps):
#         amplitude = amp[i * dx]
#         ct.table[i].waveform.index = wfm_index
#         ct.table[i].amplitude0.value = amplitude
#         ct.table[i].amplitude1.value = amplitude
#         ct.table[i].phase0.value = 0
#         ct.table[i].phase1.value = 90 + qb_pars['qb_mixer_imbalance'][1]

def ct_sweep_phase(ct,exp_pars,wfm_index,n_steps):
    
    for i in range(n_steps):
        phase = exp_pars['x0'] + i * exp_pars['dx']
        ct.table[i].waveform.index = wfm_index
        ct.table[i].phase0.value = phase
        ct.table[i].phase1.value = phase
        ct.table[i].waveform.samplingRateDivider = 0 # sets the AWG rate for the pi/2 pulses to 2.4 GHz
        
def init_ct(hdawg_core):
    
    ct_schema = hdawg_core.awgs[0].commandtable.load_validation_schema()
    ct = CommandTable(ct_schema)
    
    return ct

def arb_pulse(ct,table_index,wfm_index,amp,base_theta,theta):

    # print(amp,theta)
    ct.table[table_index].waveform.index = wfm_index
    ct.table[table_index].waveform.samplingRateDivider = 0
    ct.table[table_index].amplitude0.value = amp
    ct.table[table_index].amplitude1.value = amp
    ct.table[table_index].phase0.value = theta
    ct.table[table_index].phase1.value = base_theta + theta
    
    return ct

def determine_axes(exp_pars,qb_pars):
    print(f"Preparing {exp_pars['initial-state']} state")
    print(f"Measuring along {exp_pars['tomographic-axis']} axis")
    
    base_amp = qb_pars['pi_amp']
    
    if exp_pars['initial-state'] == 'X' and exp_pars['tomographic-axis'] == 'X':
        theta_prep = theta_tom = 0
        amp_prep = 0.5
        amp_tom = -0.5
    elif exp_pars['initial-state'] == 'Y' and exp_pars['tomographic-axis'] == 'X':
        theta_prep = 90
        theta_tom = 0
        amp_prep = 0.5
        amp_tom = -0.5
    elif exp_pars['initial-state'] == '0' and exp_pars['tomographic-axis'] == 'X':
        theta_prep = 0
        theta_tom = 180
        amp_prep = 0
        amp_tom = 0.5
    elif exp_pars['initial-state'] == '1' and exp_pars['tomographic-axis'] == 'X':
        theta_prep = 0
        theta_tom = 0
        amp_prep = 1
        amp_tom = -0.5
    elif exp_pars['initial-state'] == 'X' and exp_pars['tomographic-axis'] == 'Y':
        theta_prep = 0
        theta_tom = 90
        amp_prep = 0.5
        amp_tom = -0.5
    elif exp_pars['initial-state'] == 'Y' and exp_pars['tomographic-axis'] == 'Y':
        theta_prep = theta_tom = 90
        amp_prep = 0.5
        amp_tom = -0.5
    elif exp_pars['initial-state'] == '0' and exp_pars['tomographic-axis'] == 'Y':
        theta_prep = 0
        theta_tom = -90
        amp_prep = 0
        amp_tom = 0.5
    elif exp_pars['initial-state'] == '1' and exp_pars['tomographic-axis'] == 'Y':
        theta_prep = 0
        theta_tom = 90
        amp_prep = 1
        amp_tom = -0.5
    elif exp_pars['initial-state'] == 'X' and exp_pars['tomographic-axis'] == 'Z':
        theta_prep = 0
        theta_tom = 0
        amp_prep = 0.5
        amp_tom = 0
    elif exp_pars['initial-state'] == 'Y' and exp_pars['tomographic-axis'] == 'Z':
        theta_prep = 90
        theta_tom = 0
        amp_prep = 0.5
        amp_tom = 0
    elif exp_pars['initial-state'] == '0' and exp_pars['tomographic-axis'] == 'Z':
        theta_prep = 0
        theta_tom = 0
        amp_prep = 0
        amp_tom = 0
    elif exp_pars['initial-state'] == '1' and exp_pars['tomographic-axis'] == 'Z':
        theta_prep = 0
        theta_tom = 0
        amp_prep = 1
        amp_tom = 0
    elif exp_pars['initial-state']=='45' and exp_pars['tomographic-axis'] == 'X':
        theta_prep = 45
        theta_tom = 180
        amp_prep = 0.5
        amp_tom = 0.5
    elif exp_pars['initial-state']=='45' and exp_pars['tomographic-axis'] == 'Y':
        theta_prep = 45
        theta_tom = -90
        amp_prep = 0.5
        amp_tom = 0.5
    elif exp_pars['initial-state']=='45' and exp_pars['tomographic-axis'] == 'Z':
        theta_prep = 45
        theta_tom = 0
        amp_prep = 0.5
        amp_tom = 0
    elif exp_pars['tomographic-axis'] == 'X':
        theta_prep = 90
        theta_tom = 180
        amp_prep = 0.25
        amp_tom = 0.5
    elif exp_pars['tomographic-axis'] == 'Y':
        theta_prep = 90
        theta_tom = -90
        amp_prep = 0.25
        amp_tom = 0.5
    elif exp_pars['tomographic-axis'] == 'Z':
        theta_prep = 90
        theta_tom = 0
        amp_prep = 0.25
        amp_tom = 0
        # print('Invalid combination!!!')
    
    amp_prep = amp_prep*base_amp
    amp_tom = amp_tom*base_amp
    
    return theta_prep,theta_tom,amp_prep,amp_tom


def determine_axes_new(exp_pars,qb_pars):
    print(f"Preparing {exp_pars['initial-state']} state")
    print(f"Measuring along {exp_pars['tomographic-axis']} axis")
    
    base_amp = qb_pars['pi_amp']
    
    # In units of pi
    start_axis = exp_pars['initial-state'][0]
    rot_axis = exp_pars['initial-state'][1]
    
    # Convert starting axis to degrees to match tomographic axis:
    theta_prep = start_axis * 180

    
    if exp_pars['tomographic-axis'] == 'X':
        #theta_prep = 90
        theta_tom = 180
        #amp_prep = 0.25
        amp_tom = 0.5
    elif exp_pars['tomographic-axis'] == 'Y':
        #theta_prep = 90
        theta_tom = -90
        #amp_prep = 0.25
        amp_tom = 0.5
    elif exp_pars['tomographic-axis'] == 'Z':
        #theta_prep = 90
        theta_tom = 0
        #amp_prep = 0.25
        amp_tom = 0
    
    # Pulse amplitude for preparation and for tomography
    amp_prep = rot_axis*base_amp
    amp_tom = amp_tom*base_amp
    
    #print(theta_prep,theta_tom,amp_prep,amp_tom)
    return theta_prep,theta_tom,amp_prep,amp_tom



def determine_axes_coherence(exp_pars,qb_pars):
    '''
    Construct parameters for pulses
    
    OUTPUT:
    ------
    theta_prep (FLOAT): state preparation angle
    theta_tom (ARRAY): tomography axis angles (x,y,z)
    amp_prep (FLOAT): state preparation pulse amplitude
    amp_tom (ARRAY):  tomography pulse amplitude (x,y,z)
    '''
    print(f"Preparing {exp_pars['initial-state']} state")

    
    # base amplitude calibrated from Rabi
    base_amp = qb_pars['pi_amp']
    
    # In units of pi
    start_axis = exp_pars['initial-state'][0]
    rot_axis = exp_pars['initial-state'][1]
    
    # Convert starting axis to degrees to match tomographic axis:
    theta_prep = start_axis * 180
    #theta_tom = np.array((180,-90,0)) # x,y,z tomgraphy angle 
    #amp_tom = np.array((0.5,0.5,0)) # x,y,z tomography amplitude pulse
    
    # Pulse amplitude for preparation and for tomography
    amp_prep = rot_axis*base_amp
    #amp_tom = amp_tom*base_amp
    

    return theta_prep,amp_prep,base_amp