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
def gen_seq_code(exp,axis):

    if exp == 'spectroscopy':
        code = spec_sequence(on_off=False)
    elif exp == 't-rabi':
        code = time_rabi_sequence()
    elif exp == 'p-rabi':
        code = power_rabi_sequence()
    elif exp == 'T1':
        code = T1_sequence()
    elif exp == 'ramsey':
        code = ramsey_sequence()
    elif exp == 'echo':
        code = echo_sequence()
    elif exp == 'tomography':
        code = tomography_sequence()
    elif exp == 'coherence-stabilization':
        code = state_stabilization_sequence()
    elif exp == 'z-gate':
        code = z_gate_sequence()
    elif exp == 'single-shot':
        code = single_shot_sequence()
    
    code = finalize_sequence(code,axis)
    
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
    
    awg_program = '''playWave(1,w_marker);
            playZero(qubit_reset_time,AWG_RATE_1P2MHZ); '''
                     
    return awg_program

def tomographic_pulse_sequence(axis='Z'):
    
    if axis == 'X' or axis == 'Y':
        awg_program = '''executeTableEntry(n_steps);'''
    elif axis == 'Z':
        awg_program = ''''''
        
    return awg_program

def finalize_sequence(awg_program,axis='Z'):
    
    awg_program = awg_program.replace('_trigger_readout_',trigger_readout_sequence())
    awg_program = awg_program.replace('_tomography_pulse_',tomographic_pulse_sequence(axis))
    
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

def time_rabi_sequence(tomography=False):
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
                _tomography_pulse_
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
                executeTableEntry(n_steps);
                executeTableEntry(i);
                _tomography_pulse_
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
                    executeTableEntry(n_steps+1);
                    executeTableEntry(i);
                    executeTableEntry(n_steps+1);
                    _tomography_pulse_
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
                executeTableEntry(n_steps+1);
                executeTableEntry(i);
                executeTableEntry(n_steps+2);
                executeTableEntry(i);
                executeTableEntry(n_steps+1);
                _tomography_pulse_
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
                _trigger_readout_
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
        for (i=0; i<n_steps; i++) {
                executeTableEntry(n_steps);
                executeTableEntry(i);
                executeTableEntry(n_steps+1);
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
    repeat(512) {
        // OFF Measurement
        _trigger_readout_
        executeTableEntry(0);
        _trigger_readout_
        }'''

    return awg_program
 
def mixer_calib_sequence():
    '''Generate qubit spectroscopy sequence'''
   
    awg_program = '''
        resetOscPhase();
        wave w_gauss = amp*gauss(N, 1.0, N/2, N/2);
        wave w_zeros = zeros(N);
        
        
        playWave(1,2,w_gauss,1,2,w_zeros);
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

#%% waveform_funcs
def make_wave(pulse_type='gauss', wave_name='wave', wfm_pars = {}, amplitude = 1.0, pulse_length = 16 , output_order='1' ):
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
        pulse = amplitude*gen_arb_wfm('rising',wfm_pars,channel='I',normalize=True)
    elif pulse_type == 'arb_Q':
        pulse = amplitude*gen_arb_wfm('markov',wfm_pars,channel='Q',n_points=pulse_length)
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
                        output_order = '2'))
    )

    elif exp == 't-rabi':
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
                        output_order = '12'))
    )
    
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
        Wave(*make_wave(pulse_type = 'zero',
                        wave_name = 'w_zero_Q',
                        amplitude = 0,
                        pulse_length = n_points,
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
         amp_half = qb_pars['pi_half_amp']
         pi_amp = qb_pars['pi_amp']
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

        #   sequence.waveforms[1] = (
        #       Wave(*make_wave(pulse_type = 'pi',
        #                 wave_name = 'w_tom_I',
        #                 amplitude = 1,
        #                 pulse_length = N,
        #                 output_order = '12')), 
        #       Wave(*make_wave(pulse_type = 'zero',
        #                 wave_name = 'w_tom_Q',
        #                 amplitude = 0,
        #                 pulse_length = N,
        #                 output_order = '12'))
        # )  

         sequence.waveforms[1] = (
         Wave(*make_wave(pulse_type = 'arb_I',
                         wave_name = 'w_arb_I',
                         amplitude = 0,
                         pulse_length = n_points,
                         wfm_pars = wfm_pars,
                         output_order = '12')), 
         Wave(*make_wave(pulse_type = 'arb_Q',
                         wave_name = 'w_arb_Q',
                         amplitude = 0,
                         pulse_length = n_points,
                         wfm_pars=wfm_pars,
                         output_order = '12'))
     )

       

       
         

    elif exp == 'single-shot':
        N = qb_pars['pi_len']
        amp = qb_pars['pi_amp']
    #pi_pulse = qb_pars['pi_amp'] * gaussian (N,N/5)
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
    sequence.constants['n_avg'] = exp_pars['n_avg']
    if 'element' in exp_pars and exp_pars['element'] == 'rr':
        sequence.constants['qubit_reset_time'] = roundToBase(qb_pars['rr_reset_time']*1.17e6)
    else:
        sequence.constants['qubit_reset_time'] = roundToBase(qb_pars['qubit_reset_time']*1.17e6)
    if exp != 'spectroscopy':
        sequence.constants['n_steps'] = n_steps
        

#%% command_table_funcs
# Please refer to https://docs.zhinst.com/hdawg/commandtable/v2/schema for other settings
def make_ct(hdawg_core,exp_pars={},qb_pars={},wfm_pars={},x0=0,dx=16,n_steps=100):
   
    exp = exp_pars['exp']
    base_theta = 90 + qb_pars['qb_mixer_imbalance'][1] 
    
    if exp == 'p-rabi':
        sweep_var = 'amp'
    elif exp == 'z-gate':
        sweep_var = 'phase'
    elif exp == 'single-shot' or exp == 'coherence-stabilization':
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
        
    if exp == 'coherence-stabilization':
        wfm_index = 0
        theta_prep,theta_tom,amp_prep,amp_tom = determine_axes(exp_pars,qb_pars)
        ct_sweep_length(ct,exp,1,qb_pars,x0,dx,0,n_steps)
        ct = arb_pulse(ct,n_steps,wfm_index,amp_prep,base_theta,theta_prep) # preparation pulse
        ct = arb_pulse(ct,n_steps+1,wfm_index,amp_tom,base_theta,theta_tom) # tomography pulse
        
    elif exp == 't-rabi':
        wfm_index = 3
        theta_prep,theta_tom,amp_prep,amp_tom = determine_axes(exp_pars,qb_pars)
        ct = arb_pulse(ct,n_steps,wfm_index,amp_tom,base_theta,theta_tom) # tomography pulse
        wfm_index = 0 # gauss rise
        ct = arb_pulse(ct,n_steps+1,wfm_index,1,base_theta,0)
        wfm_index = 1 # gauss fall
        ct = arb_pulse(ct,n_steps+2,wfm_index,1,base_theta,0)
    elif exp == 'T1':
        wfm_index = 0
        ct = arb_pulse(ct,n_steps,wfm_index,1,base_theta,0) # pi_pulse
    elif exp == 'single-shot':
        wfm_index = 0
        ct = arb_pulse(ct,0,wfm_index,1,base_theta,0) # pi_pulse
    elif exp == 'ramsey':
        wfm_index = 0
        ct = arb_pulse(ct,n_steps+1,wfm_index,1,base_theta,0) # pi_pulse
        
    return ct
    
def ct_sweep_length(ct,exp,wfm_index,qb_pars,x0,dx,theta,n_steps=100):
    
    if exp == 't-rabi':
        stop = n_steps+1
    else:
        stop = n_steps
    
    for i in range(n_steps):
        wfm_length = x0 + i * dx
        # print(wfm_length)
        ct.table[i].waveform.index = wfm_index
        ct.table[i].waveform.length = wfm_length
        ct.table[i].phase0.value = 0
        ct.table[i].phase1.value = 90 + qb_pars['qb_mixer_imbalance'][1] + theta
        
def ct_sweep_amp(ct,wfm_index,exp_pars={},qb_pars={},n_steps=100):
    
    for i in range(n_steps):
        amp = exp_pars['x0'] + i * exp_pars['dx']
        ct.table[i].waveform.index = wfm_index
        ct.table[i].amplitude0.value = amp
        ct.table[i].amplitude1.value = amp
        ct.table[i].phase0.value = 0
        ct.table[i].phase1.value = 90 + qb_pars['qb_mixer_imbalance'][1]
        
def ct_arb_amp(ct,wfm_index,exp_pars={},qb_pars9wfm_pars={},n_steps=100):
    
    time_arr = np.arange(wfm_pars['t0'],wfm_pars['tmax']-wfm_pars['dt']/2,wfm_pars['dt'])
    fun = lambda x : (1/np.sqrt(wfm_pars['tb'] - x))
    amp = fun(time_arr)
    plt.plot(time_arr,amp)
    dx = int(len(time_arr)/n_steps)
    for i in range(n_steps):
        amplitude = amp[i * dx]
        ct.table[i].waveform.index = wfm_index
        ct.table[i].amplitude0.value = amplitude
        ct.table[i].amplitude1.value = amplitude
        ct.table[i].phase0.value = 0
        ct.table[i].phase1.value = 90 + qb_pars['qb_mixer_imbalance'][1]

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
        amp_prep = amp_tom = 0.5
    elif exp_pars['initial-state'] == 'Y' and exp_pars['tomographic-axis'] == 'X':
        theta_prep = 90
        theta_tom = 0
        amp_prep = amp_tom = 0.5
    elif exp_pars['initial-state'] == '0' and exp_pars['tomographic-axis'] == 'X':
        theta_prep = 0
        theta_tom = 0
        amp_prep = 0
        amp_tom = 0.5
    elif exp_pars['initial-state'] == '1' and exp_pars['tomographic-axis'] == 'X':
        theta_prep = 0
        theta_tom = 0
        amp_prep = 1
        amp_tom = 0.5
    elif exp_pars['initial-state'] == 'X' and exp_pars['tomographic-axis'] == 'Y':
        theta_prep = 0
        theta_tom = 90
        amp_prep = amp_tom = 0.5
    elif exp_pars['initial-state'] == 'Y' and exp_pars['tomographic-axis'] == 'Y':
        theta_prep = theta_tom = 90
        amp_prep = amp_tom = 0.5
    elif exp_pars['initial-state'] == '0' and exp_pars['tomographic-axis'] == 'Y':
        theta_prep = 0
        theta_tom = 90
        amp_prep = 0
        amp_tom = 0.5
    elif exp_pars['initial-state'] == '1' and exp_pars['tomographic-axis'] == 'Y':
        theta_prep = 0
        theta_tom = 90
        amp_prep = 1
        amp_tom = 0.5
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
    elif exp_pars['tomographic-axis'] == 'X':
        theta_prep = 90
        theta_tom = 0
        amp_prep = 0.25
        amp_tom = 0.5
    elif exp_pars['tomographic-axis'] == 'Y':
        theta_prep = 90
        theta_tom = 90
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