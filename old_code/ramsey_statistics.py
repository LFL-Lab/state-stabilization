# -*- coding: utf-8 -*-
from configuration import *

'''------------------------------------------------------Ramsey Statistics---------------------------------------------------'''
'''DESCRIPTION: Repeat Ramsey Measurement for a couple of hours to determine timescale of environmental fluctuations'''
try:

    list_of_files = glob.glob('E:\generalized-markovian-noise\%s\Ramsey\ramsey_statistics\*.csv'%(meas_device))
    latest_file = max(list_of_files, key=os.path.getctime)
    iteration_ramsey_statistics = int(latest_file[-7:-4].lstrip('0')) + 1
except:
    iteration_ramsey_statistics = 1

exp = 'ramsey statistics'
options_ramsey_statistics = {
    'nAverages':        256,
    'Tmax':             options_ramsey['Tmax'],
    'stepSize':         options_ramsey['stepSize'],
    'integration_length': 2.3e-6,
    'cav_resp_time':    options_ramsey['cav_resp_time'],
    'amplitude_hd':     A_d,
    'sequence':         'ramsey',
    'nSteps':           nPoints,
    'measPeriod':       options_ramsey['measPeriod'],
    'qubitDriveFreq':   options_ramsey['qubitDriveFreq'],
    'sweep':            0,
    'pi2Width':         options_ramsey['pi2Width'],
    'AC_pars':          options_ramsey['AC_pars'],
    'RT_pars':          [0,0],
}

nReps = 150
rep = np.arange(nReps)
detun_arr = np.zeros(nReps)
T_phi_arr = np.zeros(nReps)
error_arr = np.zeros(nReps)
now = np.zeros(nReps)
plot = 1

start = time.time()
for i in range(nReps):
    t,I,Q,nPoints = expf.pulse(daq,awg,setup=[1,1,1],**options_ramsey_statistics)
    detun_arr[i],T_phi_arr[i],error = pf.pulse_plot1d(awg,x_vector=t, y_vector=I,plot=plot,dt=options_ramsey_statistics['Tmax']*1e6/nPoints,**options_ramsey_statistics, iteration=iteration_ramsey_statistics)
    error_arr[i] = max(error)
    now[i] = time.time()  - start

# # discard values with high errors
detun_arr_clean = np.zeros(nReps)
T_phi_arr_clean = np.zeros(nReps)
rep_clean = np.zeros(nReps)
now_clean = np.zeros(nReps)

def condition(x): return x > 5

bad_index_arr = [idx for idx, element in enumerate(error_arr) if condition(element)]

detun_arr_clean = np.delete(detun_arr,bad_index_arr)
T_phi_arr_clean = np.delete(T_phi_arr,bad_index_arr)
rep_clean = np.delete(rep,bad_index_arr)
# time_arr = 14.5*np.arange(len(rep_clean))
time_arr = 14.5*np.arange(len(now))

fig, (ax1,ax2) = plt.subplots(2,sharex=True)
ax1.plot(time_arr,np.abs(detun_arr), '-o', markersize = 3, c='C0')
# ax1.plot(rep_clean,np.abs(detun_arr_clean), '-o', markersize = 3, c='C0')
ax1.set_ylim((0,1))
ax1.set_ylabel('$\Delta$ (MHz)')
fig.suptitle('Ramsey Statistics %03d'%(iteration_ramsey_statistics))
ax2.plot(time_arr,T_phi_arr, '-o', markersize = 3, c='C0')

# ax2.plot(rep_clean,T_phi_arr_clean, '-o', markersize = 3, c='C0')
ax2.set_xlabel('time (sec)')
ax2.set_ylabel('$T_{\phi} (\mu s)$')
ax2.set_ylim((0,5))
textstr = "$\omega_d = 2\pi\\times%.4f$ GHz" %(options_ramsey_statistics['qubitDriveFreq']*1e-9)
plt.gcf().text(1, 0.25, textstr, fontsize=14)

# save data
with open("E:\\generalized-markovian-noise\\%s\\ramsey\\ramsey_statistics\\%s_data_%03d.csv"%(meas_device,'ramsey_statistics',iteration_ramsey_statistics),"w",newline="") as datafile:
    writer = csv.writer(datafile)
    writer.writerow(options_ramsey_statistics.keys())
    writer.writerow(options_ramsey_statistics.values())
    writer.writerow(detun_arr)
    writer.writerow(T_phi_arr)
    writer.writerow(now)
    writer.writerow(error_arr)

iteration_ramsey_statistics += 1