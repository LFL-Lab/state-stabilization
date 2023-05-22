# -*- coding: utf-8 -*-
"""
Created on Wed May 10 13:54:16 2023

@author: lfl
"""

freqs = np.arange(start=3.8,stop=3.9,step=25e-6) # frequencies are in GHz

qb.exp_pars = {
    'n_avg':                2048,
    'qubit_reset_time':     100e-6,
    'amp_q':                0.3,
    'satur_dur':            40e-6,
    'rr_attenuation':       15,
    }

p_data,I,Q = qb.spectroscopy(freqs,save_data=True)
qb.spec_plot(freq=freqs,I=I,Q=Q,element='qubit',find_peaks=True)

def spectroscopy(self,freqs,save_data=True):
    '''
    DESCRIPTION: Executes qubit spectroscopy.

    '''
    self.exp = 'spectroscopy'
    self.exp_pars['fsAWG'] = 600e6
    result_length=2*self.exp_pars['n_avg']
    instfuncs.set_attenuator(self.exp_pars['rr_attenuation'])

    readout_pulse_length = 2.3e-6 + self.qb_pars['cav_resp_time'] + 2e-6

    '''Setup HDAWG and UHFQA for pulsed spectroscopy'''
    self.awg.set('/dev8233/awgs/0/outputs/*/modulation/mode', 0)
    self.awg.setInt('/dev8233/awgs/0/time',2) # sets AWG sampling rate to 600 MHz
    self.setup_awg()
    
    # Initialize sequence class
    self.sequence = Sequence()
    # setup sequence code
    self.sequence.code = '''       

    // Beginning of the core sequencer program executed on the HDAWG at run time
    wave w_marker = 2*marker(512,1);

    repeat(n_avg) {
        // OFF Measurement
        playWave(1,w_marker);
        playZero(qubit_reset_time,AWG_RATE_1P2MHZ);
        // ON measurement
        playWave(1,w_const,2,w_zero);
        playWave(1,w_marker);
        playZero(qubit_reset_time,AWG_RATE_1P2MHZ);
                }'''
    
    # setup constants
    self.sequence.constants['n_avg'] = self.exp_pars['n_avg']
    self.sequence.constants['qubit_reset_time'] = self.roundToBase(self.exp_pars['qubit_reset_time']*1.17e6)
    # setup waveforms
    self.sequence.waveforms = Waveforms()
    qubit_drive_dur = self.roundToBase(self.exp_pars['satur_dur']*self.exp_pars['fsAWG'])
    const_pulse = self.exp_pars['amp_q'] * np.ones(qubit_drive_dur)
    self.sequence.waveforms[0] = (Wave(const_pulse, name="w_const", output=OutputType.OUT1),
        Wave(np.zeros(qubit_drive_dur), name="w_zero", output=OutputType.OUT2))
      # upload everything to awg
    with self.hdawg_core.set_transaction():
        self.hdawg_core.awgs[0].load_sequencer_program(self.sequence)
        self.hdawg_core.awgs[0].write_to_waveform_memory(self.sequence.waveforms)
      
  
    self.qa.setInt('/dev2528/awgs/0/time',2)
    
    self.qa_sequence = Sequence()
    self.qa_sequence.code = """
    const cav_resp_time = delay;

    // Readout pulse
    wave readoutPulse = 0.2*ones(N);

    while(true) {
        waitDigTrigger(1,1);
        startQA();
        playWave(1,readoutPulse);
    }
    """
    
    self.qa_sequence.constants['delay'] = self.roundToBase(self.qb_pars['cav_resp_time']*450e6)
    self.qa_sequence.constants['N'] = self.roundToBase(readout_pulse_length*450e6)
    
    with self.qa_awg_core.set_transaction():
        self.qa_awg_core.awgs[0].load_sequencer_program(self.qa_sequence)
        
    self.qa.setInt('/dev2528/awgs/0/auxtriggers/0/channel', 2) # sets the source of digital trigger 1 to be the signal at trigger input 3 (back panel)
    self.qa.setDouble('/dev2528/triggers/in/2/level', 0.1)
    self.qa.setInt('/dev2528/awgs/0/auxtriggers/0/slope', 1)
    # print('-------------Configuring QA-------------\n')
    base_rate=1.8e9
    bypass_crosstalk=0
    # set modulation frequency of QA AWG to some IF and adjust input range for better resolution
    self.qa.setDouble('/dev2528/sigins/0/range',0.5)
    # self.qa.setInt('/dev2528/oscs/0/freq'.format(device),int(rr_IF))
    # self.qa.setInt('/dev2528/awgs/0/outputs/0/mode', 1) # AWG set to modulation mode

    # QA setup Settings
    self.qa.setInt('/dev2528/qas/0/integration/sources/0', 0)
    self.qa.setInt('/dev2528/qas/0/delay',round(delay*base_rate))
    # if sequence =='spec':
    #     self.qa.setInt('/dev2528/qas/0/integration/mode'.format(device), 0) # 0 for standard (4096 samples max), 1 for spectroscopy
    # elif sequence =='pulse':
    #     self.qa.setInt('/dev2528/qas/0/integration/mode'.format(device), 0)
    self.qa.setInt('/dev2528/qas/0/integration/mode', 0) # 0 for standard (4096 samples max), 1 for spectroscopy
    self.qa.setInt('/dev2528/qas/0/integration/length', 4096)
    self.qa.setInt('/dev2528/qas/0/bypass/crosstalk', bypass_crosstalk)   #No crosstalk matrix
    self.qa.setInt('/dev2528/qas/0/bypass/deskew', 1)   #No crosstalk matrix
    # self.qa.setInt('/dev2528/qas/0/bypass/rotation'.format(device), 1)   #No rotation
    # x = np.linspace(0,integration_length,round(integration_length*base_rate))
    # weights_I = np.sin(2*np.pi*rr_IF*x)
    # weights_Q = np.zeros(round(integration_length*base_rate))
    weights_I = weights_Q = np.ones(4096)
    # weights_Q = np.zeros(round(integration_length*base_rate))
    self.qa.setVector('/dev2528/qas/0/integration/weights/0/real', weights_I)
    self.qa.setVector('/dev2528/qas/0/integration/weights/0/imag', weights_Q)
    self.qa.sync()
    self.qa.setInt('/dev2528/qas/0/integration/trigger/channel', 7); # 0 for trigger input ch 1, 7 for internal triggering (QA AWG -> QA)

    # QA Monitor Settings
    self.qa.setInt('/dev2528/qas/0/monitor/trigger/channel', 7)
    self.qa.setInt('/dev2528/qas/0/monitor/averages',1)
    self.qa.setInt('/dev2528/qas/0/result/averages', 1)
    self.qa.setInt('/dev2528/qas/0/monitor/length', 4096)
    # configure triggering (0=trigger input 1 7 for internal trigger)

    # QA Result Settings
    self.qa.setInt('/dev2528/qas/0/result/length', result_length)
    self.qa.setInt('/dev2528/qas/0/result/source', source) # 2 -> source = rotation | 7 = integration
    self.qa.setInt('/dev2528/qas/0/result/reset', 1)
    self.qa.setInt('/dev2528/qas/0/result/enable', 1)
    self.qa.setInt('/dev2528/qas/0/result/mode',0) # cyclic averaging
    self.qa.sync()
    
    print('Start measurement')
    sweep_data, paths = self.create_sweep_data_dict()
    data_ON = []
    data_OFF = []

    self.enable_awg(self.qa, 'dev2528') # start the readout sequence
      
    bt = time.time()
    
    with tqdm(total = len(freqs)) as progress_bar:
        for f in freqs:
            instfuncs.set_qb_LO(f*1e9,sweep=True)
            self.qa_result_reset()
            self.qa_result_enable()
            self.enable_awg(self.awg,'dev8233') #runs the drive sequence
            data = self.acquisition_poll(paths, result_length, timeout = 2*result_length*self.exp_pars['qubit_reset_time']) # transfers data from the QA result to the API for this frequency point
            # seperate OFF/ON data and average
            data_OFF = np.append(data_OFF, np.mean([data[paths[0]][k] for k in self.even(len(data[paths[0]]))]))
            data_ON =  np.append(data_ON, np.mean([data[paths[0]][k] for k in self.odd(len(data[paths[0]]))]))
            progress_bar.update(1)

    et = time.time()
    duration = et-bt
    print(f'Measurement time: {duration:.1f} seconds')

    data = (data_ON-data_OFF)/4096
    I_data= data.real
    Q_data = data.imag

    power_data = np.abs(I_data*I_data.conjugate()+Q_data*Q_data.conjugate())

    self.enable_awg(self.awg, 'dev8233',enable=0)
    self.stop_result_unit(paths)
    self.enable_awg(self.qa, 'dev2528', enable = 0)
    
    if save_data:
        self.save_data(project,device_name,data=np.array([[freqs],[I_data],[Q_data]]))
    
    return power_data,I_data,Q_data

