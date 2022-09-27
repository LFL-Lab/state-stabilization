# -*- coding: utf-8 -*-
"""
Created on Mon May 16 11:38:14 2022

@author: Evangelos
"""

from script import meas_device,options_ramsey

'''---------------------------------------Ramsey Parameter Sweep for Random Telegraph Noise-------------------------------------------'''

'''DESCRIPTION:
    1. Generate ~200 instances of telegraph noise for all parameter points and store them in separate folders
    2. Repeat ramsey measurement with different noise realizations
    3. Every x instances, do ramsey without RTN to get statistics for the bare T2*
    4. Repeat for every parameter point
'''

try:
    directory = 'E:\\generalized-markovian-noise\\%s\\sweep_data\\ramsey\\' %(meas_device)
    latest_sweep = max(glob.glob(os.path.join(directory, '*')), key=os.path.getmtime)
    sweep_count = int(latest_sweep[-3:].lstrip('0')) + 1
except:
    sweep_count = 1

numRealizations = 128
b_measurements = 128
numIterations = numRealizations + b_measurements
interval = int(numRealizations/b_measurements) # how often to take background data without generalized markovian noise
B0 = np.linspace(0.0005,0.05,10) #the relationship between B0 and frequency of oscillations is Omega_R = 25 MHz * A_q
# B0 = [0.0005 ,0.05]
# B0 = [0.02]
# nu = np.concatenate((np.linspace(0.1,2,5),np.linspace(3,20,4),np.linspace(25,100,5),[1000]))
# nu = np.concatenate((np.linspace(0.1,2,4),np.linspace(3,100,6),[1000]))
nu = [0.01,0.1,0.5,1,2,5,10,100,1e3,5e3]
# nu = [0.01,10]
# nu = [1e3]
# tau = np.concatenate((np.linspace(0.1,2,3),np.linspace(3,20,2),[1000]))
# tau = [0.1,0.5,1.5,3,5,10,20,100]
# tau = [0.1,10]
tau = [3]

# initialize arrays for data
stepSize_background = 100e-9
Tmax_background = 10e-6
stepSize = 40e-9
Tmax = 10e-6
nPoints, nStepsBackground, pulse_length_increment, pulse_length_start = expf.calc_nSteps(sequence='ramsey',fsAWG=1.2e9,stepSize=stepSize_background,Tmax=Tmax_background)
nPoints,nSteps, pulse_length_increment, pulse_length_start = expf.calc_nSteps(sequence='ramsey',fsAWG=1.2e9,stepSize=stepSize,Tmax=Tmax)
backData, measData = expf.init_arrays(numRealizations=numRealizations,interval=interval,nPointsBackground=nStepsBackground,nPoints=nSteps)

# make appropriate directories
sweep_name = 'sweep_%03d'%(sweep_count)
parent_dir = 'E:\\generalized-markovian-noise\\%s\\sweep_data\\ramsey\\'%(meas_device)
main_path = os.path.join(parent_dir,sweep_name)
data_path = os.path.join(main_path,'data')
noise_path = os.path.join(main_path,'noise_instances')
plot_path = os.path.join(main_path,'plot_images')
try:
    os.mkdir(main_path)
    os.mkdir(data_path)
    os.mkdir(noise_path)
    os.mkdir(plot_path)
except:
    print('Directories already exist')
# generate noise instances
str_gen_noise = time.time()
expf.gen_noise_realizations(par1_arr=tau,par2_arr=nu,numRealizations=numRealizations,nPoints=nPoints,T_max=Tmax,sweep_count=sweep_count,meas_device=meas_device)
end_gen_noise = time.time()
print('Generating noise realizations took: %.1f s' %(end_gen_noise-str_gen_noise))

plt.close('all')

# sweep ramsey
optionsRamsey_par_sweep = {
    'sampling_rate':    1.2e9,
    'active_reset':     True,
    'threshold':        options_ramsey['threshold'],
    'prePulseLength':   options_ramsey['prePulseLength'],
    'postPulseLength':  options_ramsey['postPulseLength'],
    'nAverages':        128,
    'Tmax':             3e-6,
    'stepSize':         100e-9,
    'pi2Width':         1/2*pi_pulse*1e-9,
    'amplitude_hd':     A_d,
    'cav_resp_time':    options_ramsey['cav_resp_time'],
    'measPeriod':       options_ramsey['measPeriod'],
    'qubitDriveFreq':   options_ramsey['qubitDriveFreq'],
    'mu':               0.309,
    'sigma':            0.015,
    'B0':               0,
    'nu':               0,
    'tauk':             0,
    'AC_freq':          acStarkLO.getFrequency(),
    'source':           2,
    'sweep':            1,
    }

# par_sweep_time = expf.calc_sweep_time(par1=nu, par2=tau,measTimeBackground=2.5,measTime=5.3,nMeasBackground=b_measurements,nMeas=numRealizations)
# print('Estimated Sweep Time: %.1f hours = %.1f days'%(par_sweep_time/3600,par_sweep_time/(3600*24)))
start_sweep = time.time()
# generate data
errors = 0

for i in range(0,2):
    for j in range(len(nu)):
        for k in range(len(tau)):

            # load noise instances
            noise_realizations = B0[i]*expf.pull_wfm(sweep_name,nu[j],tau[k])
            start_point_t = time.time()
            filename = 'B0_%d_uV_nu_%d_Hz_tau_%d_ns' %(round(B0[i]*1e6),round(nu[j]*1e3),round(tau[k]*1e3))
            b = 0 # keeps track of background measurements
            n = 0 # keeps track of noisy measurements
            m = 0 # keeps track of overall measurements

            #calibrate mixers
            expf.mixer_calib(sa,awg,'fine',mixer='qubit',threshold=-75,fc=optionsRamsey_par_sweep['qubitDriveFreq'],amp=options_rabi['amplitude_hd'])
            if optionsRamsey_par_sweep['mu'] != 0:
                expf.mixer_calib(sa,awg,'fine',mixer='ac',threshold=-75,fc=optionsRamsey_par_sweep['AC_freq'],amp=options_rabi['mu'])
            else:
                pass
            attn.setAttenuation(devID = 1, atten = 0)
            expf.mixer_calib(sa,daq,'fine',mixer='readout',threshold=-75,fc=7.2581e9,amp=0.7)
            attn.setAttenuation(devID = 1, atten = attenuation)

            #calibrate pi_pulse and set threshold
            t,data,nPoints = expf.pulse(daq,awg,setup=[0,1,0],sequence='rabi',**options_rabi)
            fitted_pars,error = pf.fit_data(x_vector=t,y_vector=data,sequence='rabi',dt=t[-1]/nPoints,**options_rabi)
            optionsRamsey_par_sweep['pi2Width'] = np.round(1/4*fitted_pars[1])*1e-9
            optionsRamsey_par_sweep['threshold'] = round(np.mean(data)*2**12)

            while m < numIterations:
                if m % (interval + 1) == 0 and b <= b_measurements: # executes a background measurement every n=interval noisy measurements
                    optionsRamsey_par_sweep['nAverages'] = 128
                    optionsRamsey_par_sweep['B0'] = 0
                    optionsRamsey_par_sweep['Tmax'] = Tmax_background
                    optionsRamsey_par_sweep['stepSize'] = stepSize_background
                    print('----------------------------------\nExecuting background Ramsey measurement')
                    try:
                        t2,data,nPoints = expf.pulse(daq,awg,setup=[0,1,0],sequence='ramsey',**optionsRamsey_par_sweep)
                        backData[b,:] = data
                        b += 1
                        m += 1
                    except:
                        print('Measurement Failed')
                        errors += 1
                else:
                    if m % (interval + 1) == 1:
                        setup = [0,1,0] # upload a new seqc file, and re-configure QA unit
                    else:
                        setup = [2,1,1] # replace noisy waveforms only; everything else stays the same
                    optionsRamsey_par_sweep['B0'] = B0[i]
                    optionsRamsey_par_sweep['nu'] = nu[j]
                    optionsRamsey_par_sweep['tauk'] = tau[k]
                    optionsRamsey_par_sweep['Tmax'] = Tmax
                    optionsRamsey_par_sweep['nAverages'] = 128
                    optionsRamsey_par_sweep['stepSize'] = stepSize
                    print('----------------------------------\nImplementing noise realization %d' %(n+1))
                    noise_instance = noise_realizations[n,:]
                    try:
                        t1,data,nPoints = expf.pulse(daq,awg,setup=setup,sequence='ramsey',noise_instance=noise_instance,**optionsRamsey_par_sweep)
                        measData[n,:] = data
                        n += 1
                        m += 1
                    except:
                        print('Measurement Failed')
                        errors += 1

            # save data after each parameter sweep point
            with open("E:\\generalized-markovian-noise\\%s\\sweep_data\\ramsey\\%s\\data\\data_%s.csv"%(meas_device,sweep_name,filename),"w",newline="") as datafile:
                writer = csv.writer(datafile)
                writer.writerow(optionsRamsey_par_sweep.keys())
                writer.writerow(optionsRamsey_par_sweep.values())
                writer.writerow(['Background Time Data'])
                writer.writerow(t2)
                writer.writerow(['Time Data'])
                writer.writerow(t1)
                writer.writerow(['Background Data'])
                writer.writerows(backData)
                writer.writerow(['Data'])
                writer.writerows(measData)

            end_point_t = time.time()
            print('Parameter point took %s\n%d errors'%(end_point_t-start_point_t,errors))

end_sweep = time.time()
print('Total Sweep Duration: %.1f s = %.1f hours = %.1f days' %(end_sweep-start_sweep,(end_sweep-start_sweep)/3600,(end_sweep-start_sweep)/(3600*24)))
sweep_count += 1
