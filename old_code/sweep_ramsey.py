# -*- coding: utf-8 -*-
from configuration import *

'''---------------------------------------Ramsey Parameter Sweep for Random Telegraph Noise-------------------------------------------'''

'''DESCRIPTION:
    1. Generate ~200 instances of telegraph noise for all parameter points and store them in separate folders
    2. Repeat ramsey measurement with different noise realizations
    3. Every 10 instances, do ramsey without RTN to get statistics for the bare T2*
    4. Repeat for every parameter point
'''
try:
    directory = 'E:\\generalized-markovian-noise\\%s\\sweep_data\\ramsey\\' %(meas_device)
    latest_sweep = max(glob.glob(os.path.join(directory, '*')), key=os.path.getmtime)
    sweep_count = int(latest_sweep[-3:].lstrip('0')) + 1
except:
    sweep_count = 1


numRealizations = 128
b_measurements = 32
numIterations = numRealizations + b_measurements
interval = int(numRealizations/b_measurements) # how often to take background data without generalized markovian noise
B0 = np.linspace(0.0005,0.05,11) #the relationship between B0 and frequency of oscillations is Omega_R = 25 MHz * A_q
# B0 = [0.0005,0.015]
# B0 = [0.05]
# nu = np.linspace(0.05,10,5)
# nu = [0.01,10]
nu = [0]
tau = np.concatenate((np.linspace(0.1,2,7),np.linspace(3,20,4)))
# tau = [0.1,1,2,3,10,10000]
# tau = [0.1,1000]
# tau = [10]

# initialize arrays for data
stepSize_background = 40e-9
Tmax_background = 8e-6
stepSize = 40e-9
Tmax = 15.5e-6
nPoints, nStepsBackground, pulse_length_increment, pulse_length_start = expf.calc_nSteps(sequence='ramsey',fsAWG=1.2e9,stepSize=stepSize_background,Tmax=Tmax_background)
nPoints,nSteps, pulse_length_increment, pulse_length_start = expf.calc_nSteps(sequence='ramsey',fsAWG=1.2e9,stepSize=stepSize,Tmax=Tmax)
bData_I, bData_Q, data_I, data_Q = expf.init_arrays(numRealizations=numRealizations,interval=interval,nPointsBackground=nStepsBackground,nPoints=nSteps)

# make appropriate directories
sweep_name = 'sweep_%03d'%(sweep_count)
parent_dir = 'E:\\generalized-markovian-noise\\%s\\sweep_data\\ramsey\\'%(meas_device)
main_path = os.path.join(parent_dir,sweep_name)
data_path = os.path.join(main_path,'data')
plot_path = os.path.join(main_path,'plot_images')
os.mkdir(main_path)
os.mkdir(data_path)
os.mkdir(plot_path)

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
    'nAverages':        256,
    'Tmax':             3e-6,
    'stepSize':         100e-9,
    'pi2Width':         1/2*pi_pulse*1e-9,
    'amplitude_hd':     A_d,
    'sequence':         'ramsey',
    'cav_resp_time': options_ramsey['cav_resp_time'],
    'sweep':            1,
    'measPeriod':       400e-6,
    'qubitDriveFreq':   3.3313e9,
    'AC_pars':          [0.315,0.015],
    'RT_pars':          [0,0,0],
    'AC_freq':          acStarkLO.getFrequency()
    }

par_sweep_time = expf.calc_sweep_time(par1=B0, par2=tau,measTimeBackground=5,measTime=5,nMeasBackground=b_measurements,nMeas=numRealizations)
print('Estimated Sweep Time: %.1f hours = %.1f days'%(par_sweep_time/3600,par_sweep_time/(3600*24)))
start_sweep = time.time()
# generate data
iteration_rabi = 1
for i in range(len(B0)):
    # optionsRamsey_par_sweep['RT_pars'][0] = B0[i]
    for j in range(len(tau)):
        # load noise instances
        noise_realizations = expf.pull_wfm(sweep_name=sweep_name, RT_pars=[B0[i],tau[j],nu[0]])
        start_point_t = time.time()
        filename = 'B0_%d_uV_nu_%d_kHz_tau_%d_ns' %(round(B0[i]*1e6),round(nu[0]*1e3),round(tau[j]*1e3))
        # optionsRamsey_par_sweep['RT_pars'][1] = tau[j]
        a = 0 # keeps track of background measurements
        b = 0 # keeps track of noisy measurements
        k = 0
        #calibrate pi_pulse
        t,I,Q,nPoints = expf.pulse(daq,awg,setup=[0,0,0],**options_rabi)
        fitted_pars,error = pf.fit_data(x_vector=t,y_vector=I,dt=t[-1]/nPoints,**options_rabi)
        optionsRamsey_par_sweep['pi2Width'] = np.round(1/4*fitted_pars[1])*1e-9
        optionsRamsey_par_sweep['threshold'] = round(np.mean([max(I),min(I)])*2**12)
        # fig = pf.plot_data(awg,x_vector=t,y_vector=I,fitted_pars=fitted_pars,**options_rabi,iteration=iteration_rabi)
        # plt.savefig(os.path.join(plot_path,filename+'rabi_fig_%03d.png' %(iteration_rabi)) , bbox_inches='tight')
        # plt.close(fig)
        iteration_rabi += 1
        # print('Next parameter point: B_0 = %.5f V and tau = %.3f microseconds' %(B0[i],tau[j]))
        while k < numIterations:
            if k % (interval + 1) == 0 and a != b_measurements:
                # get background T2* every 10 or so measurements
                optionsRamsey_par_sweep['nAverages'] = 256
                optionsRamsey_par_sweep['RT_pars'] = [0,0,0]
                optionsRamsey_par_sweep['Tmax'] = Tmax_background
                optionsRamsey_par_sweep['stepSize'] = stepSize_background
                print('----------------------------------\nExecuting background Ramsey measurement')
                t2,I,Q,nPoints = expf.pulse(daq,awg,setup=[0,1,0],sweep_name=sweep_name,**optionsRamsey_par_sweep)
                bData_I[a,:] = I
                bData_Q[a,:] = Q
                # fitted_pars,error = pf.fit_data(x_vector=t2,y_vector=I,dt=t2[-1]/nPoints,**optionsRamsey_par_sweep)
                # fig = pf.plot_data(awg,x_vector=t2,y_vector=I,fitted_pars=fitted_pars,dt=t2[-1]/nPoints,**optionsRamsey_par_sweep,iteration=a)
                # save and then clear plot
                # plt.savefig(os.path.join(plot_path,filename+'background_meas_fig_%03d.png'%(a+1)),bbox_inches='tight')
                # plt.close(fig)
                a += 1
                print('End measurement\n----------------------------------' )
            else:
                if k % (interval + 1) == 1:
                    setup = [0,1,0]
                else:
                    setup = [2,1,1]
                optionsRamsey_par_sweep['RT_pars'] = [B0[i],tau[j],nu[0]]
                # optionsRamsey_par_sweep['RT_pars'] = [B0,tau[j],nu[i]]
                optionsRamsey_par_sweep['Tmax'] = Tmax
                optionsRamsey_par_sweep['nAverages'] = 256
                optionsRamsey_par_sweep['stepSize'] = stepSize
                print('----------------------------------\nStart %s measurement' %("ramsey"))
                print('Implementing noise realization %d' %(b+1))
                noise_instance = noise_realizations[b,:]
                t1,I,Q,nPoints = expf.pulse(daq,awg,setup=setup,sweep_name=sweep_name,noise_instance=noise_instance,**optionsRamsey_par_sweep)
                data_I[b,0:nPoints] = I
                data_Q[b,0:nPoints] = Q
                # fitted_pars,error = pf.fit_data(x_vector=t1,y_vector=I,dt=t1[-1]/nPoints,**optionsRamsey_par_sweep)
                # fig = pf.plot_data(awg,x_vector=t1,y_vector=I,fitted_pars=fitted_pars,dt=t1[-1]/nPoints,**optionsRamsey_par_sweep,iteration=b)
                # save and then clear plot
                # plt.savefig(os.path.join(plot_path,filename+'_fig_%03d.png'%(b+1)),bbox_inches='tight')
                # plt.close(fig)
                b += 1
                print('End measurement\n----------------------------------' )
            k += 1


        # save data after each parameter sweep point
        with open("E:\\generalized-markovian-noise\\%s\\sweep_data\\ramsey\\%s\\data\\data_%s.csv"%(meas_device,sweep_name,filename),"w",newline="") as datafile:
            writer = csv.writer(datafile)
            writer.writerow(optionsRamsey_par_sweep.keys())
            writer.writerow(optionsRamsey_par_sweep.values())
            writer.writerow(['Background Time Data'])
            writer.writerow(t2)
            writer.writerow(['Time Data'])
            writer.writerow(t1)
            writer.writerow(['Background Data: Channel 1'])
            writer.writerows(bData_I)
            writer.writerow(['Background Data: Channel 2'])
            writer.writerows(bData_Q)
            writer.writerow(['Data: Channel 1'])
            writer.writerows(data_I)
            writer.writerow(['Data: Channel 2'])
            writer.writerows(data_Q)

        end_point_t = time.time()
        print('Parameter point took %s'%(end_point_t-start_point_t))
        plt.close('all')

end_sweep = time.time()
print('Total Sweep Duration: %.1f s = %.1f hours = %.1f days' %(end_sweep-start_sweep,(end_sweep-start_sweep)/3600,(end_sweep-start_sweep)/(3600*24)))
sweep_count += 1