# -*- coding: utf-8 -*-

from configuration import *

'''------------------------------------------------------Sweep Rabi Amplitude----------------------------------------------'''
iteration_rabi_sweep = 2
exp = 'sweep_rabi_amp'
options_rabi_sweep = {
    'qubitDriveFreq':   3.8773e9,
    'nAverages':        256,
    'Tmax':             1.5e-6,
    'amplitude_hd':     1.0,
    'sequence':         'rabi',
    'channel':          0,
    'measPeriod':       200e-6,
    'AC_pars':          [0.4,0]
    }

amp = np.linspace(0.01,1,100)
Omega_Rabi = np.zeros(len(amp))

exp_pars = [options_rabi_sweep['qubitDriveFreq'],options_rabi_sweep['AC_pars']]
data = np.zeros((450,len(amp)))

for i in range(len(amp)):
    options_rabi_sweep['amplitude_hd'] = amp[i]
    if i == 0:
        setup = [0,1,0]
    else:
        setup = [0,1,1]
    # generate data
    print('$A_d$ = % V'%(amp))
    t,data[:,i],ch2Data,nPoints = expf.pulse(daq,awg,qubitLO,setup=[0,1,0],**options_rabi_sweep)
    pi_pulse,error = pf.pulse_plot1d(sequence='rabi',dt=options_rabi_sweep['Tmax']*1e6/nPoints,plot=1,qubitDriveFreq=options_rabi_sweep['qubitDriveFreq'],amplitude_hd=options_rabi_sweep['amplitude_hd'],x_vector=t, y_vector=data[:,i],fitting=1,AC_pars=options_rabi_sweep['AC_pars'],iteration=iteration_rabi_sweep)
    Omega_Rabi[i] = 1/(2*pi_pulse)*1e9

with open("E:\\generalized-markovian-noise\\rabi\\%s_data_%03d.csv"%(exp,iteration_rabi_sweep),"w",newline="") as datafile:
    writer = csv.writer(datafile)
    writer.writerow(exp_pars)
    writer.writerow(amp)
    writer.writerow(Omega_Rabi)
    writer.writerows(data)

def line(x,a,b):
    return a*x+b
best_vals, covar = scy.optimize.curve_fit(line, amp[0:25],Omega_Rabi[0:25]/1e6,xtol=1e-6,maxfev=3000)

fig, ax1 = plt.subplots()
plt.xticks(np.arange(-0.1,1.1,step=0.2))
plt.yticks(np.arange(0,11,step=2))
left,bottom,width,height = [0.5, 0.25, 0.3, 0.4]
ax2 = fig.add_axes([left,bottom,width,height])
ax1.plot(amp,Omega_Rabi/1e6, '-o', markersize = 3, c='C0')
ax1.set_xlabel('Qubit Drive Amplitude (Volts)')
ax1.set_ylabel('$\Omega_R$ (MHz)')
ax2.plot(amp[0:40],Omega_Rabi[0:40]/1e6,amp[0:40],line(amp[0:40],best_vals[0],best_vals[1]))
ax2.set_xlabel('Qubit Drive Amplitude (Volts)',fontsize=9)
ax2.set_ylabel('$\Omega_R$ (MHz) ',fontsize=10)
plt.xticks(np.arange(0.1,0.5,step=0.1),fontsize=9)
plt.yticks(np.arange(1,11,step=2),fontsize=9)
inset_box_txt = '$\Omega_R=$'+"{:.2e}".format(best_vals[0])+'$\\times A_d +$' +"{:.2e}".format(best_vals[1])
plt.gcf().text(0.5, 0.675, inset_box_txt, fontsize=10)
plt.show()
