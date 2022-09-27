#!/usr/bin/env python
# coding: utf-8

# ## Import some packages

# In[1]:


import textwrap
import zhinst.ziPython as zp
import zhinst.utils as zu
import numpy as np
import time
import import_ipynb 
import HDAWG as hd
import matplotlib.pyplot as plt 
get_ipython().run_line_magic('matplotlib', 'notebook')


# ## Connect instruments

# In[1]:


device_hd_id = 'dev8233'#'dev8248'#'dev8288'
host = 'localhost' #'10.42.0.225'#'10.42.34.48'#

daq = zp.ziDAQServer(host, 8004, 10)
daq.connectDevice(device_hd_id, 'usb') # connection: '1gbe' or 'usb'
device = device_hd_id


# In[ ]:


# daq.sync()
# time.sleep(1)
# data=hd.data_acquisation_cnts(daq,device,cnt=0,time_recording=1)
# t=data['/dev8288/cnts/0/sample']['timestamp']
# print(len(t))
# dt=np.diff(t)/1.8e9*1e9
# # plt.plot(time[0:130000])
# plt.hist(dt[0:130000],bins=1000)
# plt.yscale('log')


# ## Defined functions

# In[ ]:


def create_and_compile_awg(daq, device, awg_program, seqr_index= 0, timeout=1):
    awgModule = daq.awgModule()
    awgModule.set('device', device)
    awgModule.set('index', seqr_index)
    awgModule.execute()
    """Compile and upload awg_program as .elf file"""
    print("Starting compilation.")
    awgModule.set('compiler/sourcestring', awg_program)
    compilerStatus = -1
    while compilerStatus == -1:
        compilerStatus = awgModule.getInt('compiler/status')
        time.sleep(0.1)
    compilerStatusString = awgModule.getString('compiler/statusstring')
    print(f"Compiler messages:\n--------------\n{compilerStatusString}\n--------------")
    if compilerStatus == 1: # compilation failed
        raise Exception("Compilation failed.")
    if compilerStatus == 0:
        print("Compilation successful with no warnings.")
    if compilerStatus == 2:
        print("Compilation successful with warnings.")
    print("Waiting for the upload to the instrument.")
    elfProgress = 0
    elfStatus = 0
    lastElfProgressPrc = None
    while (elfProgress < 1.0) and (elfStatus != 1):
        elfProgress = awgModule.getDouble('progress')
        elfStatus = awgModule.getInt('elf/status')
        elfProgressPrc = round(elfProgress * 100);
        if elfProgressPrc != lastElfProgressPrc:
            print(f'Upload progress: {elfProgressPrc:2.0f}%')
            lastElfProgressPrc = elfProgressPrc
        time.sleep(0.1)
    if elfStatus == 0:
        print("Upload to the instrument successful.")
    if elfStatus == 1:
        raise Exception("Upload to the instrument failed.")    


# ## Simple play wave 

# In[ ]:


awg_program = textwrap.dedent(
'''
const n_sample  = 1024;
const amplitude = 0.5; 
const position  = n_sample/2;
const width     = n_sample/8;
wave w_gauss = gauss(n_sample, amplitude, position, width);

while(true){
playWave(w_gauss);
playZero(n_sample);
}
''')
daq.setInt(f'{device_hd_id}/system/awg/channelgrouping', 0)

for i in range(2):
    daq.setInt(f'/{device_hd_id}/sigouts/{i}/on', 1)


create_and_compile_awg(daq, device_hd_id, awg_program, 0, 5) # compile the seqc codes, timeout here is 10s. 
        
daq.setInt(f'{device_hd_id}/awgs/0/enable', 1)
print('enable hd_awg')
time.sleep(20)
daq.setInt(f'{device_hd_id}/awgs/0/enable', 0)
for i in range(2):
    daq.setInt(f'/{device_hd_id}/sigouts/{i}/on', 0)

print('disable hd_awg')


# ## Replace waveforms in 1x4 mode for HDAWG4

# In[ ]:


awg_program = textwrap.dedent(
'''
wave w1 = placeholder(20000, true);// true means using marker, default is false, available with labOne 20.07
wave w2 = placeholder(20000, true);
wave w3 = placeholder(20000, true);
wave w4 = placeholder(20000, true);

while(true){
playWave(1, w1, 2, w2, 3, w3, 4, w4);
waitWave();
wait(87482);
}

''')
daq.setInt(f'{device_hd_id}/system/awg/channelgrouping', 1)
for i in range(2):
    daq.setInt(f'/{device_hd_id}/sigouts/{i}/on', 1)

    
create_and_compile_awg(daq, device_hd_id, awg_program, 0, 5) # compile the seqc codes, timeout here is 10s.

# this is the configuration for grouping 1x8
N = 20000
w_1 = 0.5 * np.ones(N, dtype=float)   # wave 1 in physical output 1 of 1 AWG core
w_2 = 0.75 * np.ones(N, dtype=float)  # wave 2 in physical output 1 of 1 AWG core

marker = np.zeros(N)              
marker[3200: N-3200] = 1              # Marker setting for 1 AWG core

w_new1 = zu.convert_awg_waveform(wave1=w_1, wave2= w_2, markers=marker)

for i in range(2): 
    # here '.../waves/0' in fact for both channel 0 and 1, because the w_new1 includes 2 waveforms
    daq.setVector(f'/{device_hd_id:s}/awgs/{i}/waveform/waves/0', w_new1) 
    
        
daq.setInt(f'{device_hd_id}/awgs/0/enable', 1)
print('enable hd_awg')
time.sleep(10)
daq.setInt(f'{device_hd_id}/awgs/0/enable', 0)
print('disable hd_awg')
for i in range(2):
    daq.setInt(f'/{device_hd_id}/sigouts/{i}/on', 0)


# ## example_awg

# In[4]:


"""
Zurich Instruments LabOne Python API Example

Demonstrate how to connect to a Zurich Instruments HDAWG and upload and run an
AWG program.
"""

# Copyright 2018 Zurich Instruments AG

import os
import time
import textwrap
import numpy as np
import zhinst.utils


def run_example(device_id):
    """
    Run the example: Connect to a Zurich Instruments HDAWG upload and run a
    basic AWG sequence program. It also demonstrates how to upload (replace) a
    waveform without changing the sequencer program.

    Requirements:

       HDAWG Instrument.

    Arguments:

      device_id (str): The ID of the device to run the example with. For
        example, `dev8006` or `hdawg-dev8006`.

    Returns:

      No return value.

    Raises:

      Exception: If the device is not an HDAWG.

      RuntimeError: If the device is not "discoverable" from the API.

    See the "LabOne Programming Manual" for further help, available:
      - On Windows via the Start-Menu:
        Programs -> Zurich Instruments -> Documentation
      - On Linux in the LabOne .tar.gz archive in the "Documentation"
        sub-folder.
    """

    # Settings
    apilevel_example = 6  # The API level supported by this example.
    err_msg = "This example can only be ran on an HDAWG."
    # Call a zhinst utility function that returns:
    # - an API session `daq` in order to communicate with devices via the data server.
    # - the device ID string that specifies the device branch in the server's node hierarchy.
    # - the device's discovery properties.
#     (daq, device, _) = zhinst.utils.create_api_session(
#         device_id, apilevel_example, required_devtype="HDAWG", required_err_msg=err_msg
#     )
    zhinst.utils.api_server_version_check(daq)

    # Create a base configuration: Disable all available outputs, awgs, demods, scopes,...
    zhinst.utils.disable_everything(daq, device)

    # 'system/awg/channelgrouping' : Configure how many independent sequencers
    #   should run on the AWG and how the outputs are grouped by sequencer.
    #   0 : 4x2 with HDAWG8; 2x2 with HDAWG4.
    #   1 : 2x4 with HDAWG8; 1x4 with HDAWG4.
    #   2 : 1x8 with HDAWG8.
    # Configure the HDAWG to use one sequencer with the same waveform on all output channels.
    daq.setInt(f"/{device}/system/awg/channelgrouping", 0)

    # Some basic device configuration to output the generated wave.
    out_channel = 0
    awg_channel = 0
    amplitude = 1.0

    exp_setting = [
        ["/%s/sigouts/%d/on" % (device, out_channel), 1],
        ["/%s/sigouts/%d/range" % (device, out_channel), 1],
        ["/%s/awgs/0/outputs/%d/amplitude" % (device, awg_channel), amplitude],
        ["/%s/awgs/0/outputs/0/modulation/mode" % device, 0],
        ["/%s/awgs/0/time" % device, 0],
        ["/%s/awgs/0/userregs/0" % device, 0],
    ]
    daq.set(exp_setting)
    # Ensure that all settings have taken effect on the device before continuing.
    daq.sync()

    # Number of points in AWG waveform
    AWG_N = 2000

    # Define an AWG program as a string stored in the variable awg_program, equivalent to what would
    # be entered in the Sequence Editor window in the graphical UI.
    # This example demonstrates four methods of definig waveforms via the API
    # - (wave w0) loaded directly from programmatically generated CSV file wave0.csv.
    #             Waveform shape: Blackman window with negative amplitude.
    # - (wave w1) using the waveform generation functionalities available in the AWG Sequencer language.
    #             Waveform shape: Gaussian function with positive amplitude.
    # - (wave w2) using the vect() function and programmatic string replacement.
    #             Waveform shape: Single period of a sine wave.
    # - (wave w3) directly writing an array of numbers to the AWG waveform memory.
    #             Waveform shape: Sinc function. In the sequencer language, the waveform is initially
    #             defined as an array of zeros. This placeholder array is later overwritten with the
    #             sinc function.

    awg_program = textwrap.dedent(
        """\
        const AWG_N = _c1_;
        wave w0 = "wave0";
        wave w1 = gauss(AWG_N, AWG_N/2, AWG_N/20);
        wave w2 = vect(_w2_);
        wave w3 = zeros(AWG_N);
        while(getUserReg(0) == 0);
        setTrigger(1);
        setTrigger(0);
        playWave(w0);
        playWave(w1);
        playWave(w2);
        playWave(w3);
        """
    )

    # Define an array of values that are used to write values for wave w0 to a CSV file in the module's data directory
    waveform_0 = -1.0 * np.blackman(AWG_N)

    # Define an array of values that are used to generate wave w2
    waveform_2 = np.sin(np.linspace(0, 2 * np.pi, 96))

    # Fill the waveform values into the predefined program by inserting the array
    # as comma-separated floating-point numbers into awg_program.
    # Warning: Defining waveforms with the vect function can increase the code size
    #          considerably and should be used for short waveforms only.
    awg_program = awg_program.replace("_w2_", ",".join([str(x) for x in waveform_2]))

    # Fill in the integer constant AWG_N
    awg_program = awg_program.replace("_c1_", str(AWG_N))

    # Create an instance of the AWG Module
    awgModule = daq.awgModule()
    awgModule.set("device", device)
    awgModule.execute()

    # Get the modules data directory
    data_dir = awgModule.getString("directory")
    # All CSV files within the waves directory are automatically recognized by the AWG module
    wave_dir = os.path.join(data_dir, "awg", "waves")
    if not os.path.isdir(wave_dir):
        # The data directory is created by the AWG module and should always exist. If this exception is raised,
        # something might be wrong with the file system.
        raise Exception(f"AWG module wave directory {wave_dir} does not exist or is not a directory")
    # Save waveform data to CSV
    csv_file = os.path.join(wave_dir, "wave0.csv")
    np.savetxt(csv_file, waveform_0)

    # Transfer the AWG sequence program. Compilation starts automatically.
    awgModule.set("compiler/sourcestring", awg_program)
    # Note: when using an AWG program from a source file (and only then), the compiler needs to
    # be started explicitly with awgModule.set('compiler/start', 1)
    while awgModule.getInt("compiler/status") == -1:
        time.sleep(0.1)

    if awgModule.getInt("compiler/status") == 1:
        # compilation failed, raise an exception
        raise Exception(awgModule.getString("compiler/statusstring"))

    if awgModule.getInt("compiler/status") == 0:
        print("Compilation successful with no warnings, will upload the program to the instrument.")
    if awgModule.getInt("compiler/status") == 2:
        print("Compilation successful with warnings, will upload the program to the instrument.")
        print("Compiler warning: ", awgModule.getString("compiler/statusstring"))

    # Wait for the waveform upload to finish
    time.sleep(0.2)
    i = 0
    while (awgModule.getDouble("progress") < 1.0) and (awgModule.getInt("elf/status") != 1):
        print(f"{i} progress: {awgModule.getDouble('progress'):.2f}")
        time.sleep(0.2)
        i += 1
    print(f"{i} progress: {awgModule.getDouble('progress'):.2f}")
    if awgModule.getInt("elf/status") == 0:
        print("Upload to the instrument successful.")
    if awgModule.getInt("elf/status") == 1:
        raise Exception("Upload to the instrument failed.")

    # Replace the waveform w3 with a new one.
    waveform_3 = np.sinc(np.linspace(-6 * np.pi, 6 * np.pi, AWG_N))
    # Let N be the total number of waveforms and M>0 be the number of waveforms defined from CSV files. Then the index
    # of the waveform to be replaced is defined as following:
    # - 0,...,M-1 for all waveforms defined from CSV file alphabetically ordered by filename,
    # - M,...,N-1 in the order that the waveforms are defined in the sequencer program.
    # For the case of M=0, the index is defined as:
    # - 0,...,N-1 in the order that the waveforms are defined in the sequencer program.
    # Of course, for the trivial case of 1 waveform, use index=0 to replace it.
    # The list of waves given in the Waveform sub-tab of the AWG Sequencer tab can be used to help verify the index of
    # the waveform to be replaced.
    # Here we replace waveform w3, the 4th waveform defined in the sequencer program. Using 0-based indexing the
    # index of the waveform we want to replace (w3, a vector of zeros) is 3:
    # index of the waveform we want to replace (w3, a vector of zeros) is 3:
    index = 3
    waveform_native = zhinst.utils.convert_awg_waveform(waveform_3)
    path = f"/{device:s}/awgs/0/waveform/waves/{index:d}"
    daq.setVector(path, waveform_native)

    print(f"Enabling the AWG: Set /{device}/awgs/0/userregs/0 to 1 to trigger waveform playback.")
    # This is the preferred method of using the AWG: Run in single mode continuous waveform playback is best achieved by
    # using an infinite loop (e.g., while (true)) in the sequencer program.
    daq.setInt(f"/{device}/awgs/0/single", 1)
    daq.setInt(f"/{device}/awgs/0/enable", 1)
run_example(device_hd_id)


# ## example_precompensation_curve_fit

# In[ ]:


"""
Zurich Instruments LabOne Python API Example

Demonstrate how to connect to a Zurich Instruments HDAWG and
use the precompensation module to fit filter parameters for a
measured signal
"""

# Copyright 2019 Zurich Instruments AG

import time
import numpy as np
import zhinst.ziPython
import matplotlib.pyplot as plt


def get_precompensated_signal(module_handle, input_signal, amplitude, timeconstant):
    """
    Uploads the input_signal to the precompensationAdvisor module and returns the
    simulated forward transformed signal with an exponential filter(amplitude,timeconstant).
    """
    module_handle.set("exponentials/0/amplitude", amplitude)
    module_handle.set("exponentials/0/timeconstant", timeconstant)
    module_handle.set("wave/input/inputvector", input_signal)
    return np.array(module_handle.get("wave/output/forwardwave", True)["/wave/output/forwardwave"][0]["x"])


def run_example(device_id, do_plot=True):
    """
    Run the example: Connect to a Zurich Instruments HDAWG. The example uploads a signal to
    the precompensationAdvisor module and reads back the filtered signal. This functionality
    is used to feed a fitting algorithm for fitting filter parameters.

    Requirements:
      HDAWG


    Arguments:

      device_id (str): The ID of the device to run the example with. For
        example, `dev8050`.

      do_plot (bool, optional): Specify whether to plot the initial, target and fitted signals.


    See the "LabOne Programming Manual" for further help, available:
      - On Windows via the Start-Menu:
        Programs -> Zurich Instruments -> Documentation
      - On Linux in the LabOne .tar.gz archive in the "Documentation"
        sub-folder.
    """
    # Settings
    apilevel_example = 6  # The API level supported by this example.
    err_msg = "This example can only be ran on an HDAWG."
    # Call a zhinst utility function that returns:
    # - an API session `daq` in order to communicate with devices via the data server.
    # - the device ID string that specifies the device branch in the server's node hierarchy.
    # - the device's discovery properties.
#     (daq, device, _) = zhinst.utils.create_api_session(
#         device_id, apilevel_example, required_devtype="HDAWG", required_err_msg=err_msg
#     )

    zhinst.utils.api_server_version_check(daq)

    # Create a base configuration: Disable all available outputs, awgs, demods, scopes,...
    zhinst.utils.disable_everything(daq, device)

    pre = daq.precompensationAdvisor()

    sampling_rate = 2.4e9

    x, target_signal = generate_target_signal(sampling_rate=sampling_rate)
    actual_signal = generate_actual_signal(target_signal, sampling_rate=sampling_rate)

    # prepare the precompensationAdvisor module
    pre.set("exponentials/0/enable", 1)
    pre.set("wave/input/source", 3)
    pre.set("device", device_id)
    daq.setDouble("/" + device_id + "/system/clocks/sampleclock/freq", sampling_rate)
    # a short pause is needed for the precompensationAdvisor module to read
    # the updated the sampling rate from the device node
    time.sleep(0.05)
    sampling_rate = pre.getDouble("samplingfreq")

    # Fitting the parameters
    from lmfit import Model

    gmodel = Model(get_precompensated_signal, independent_vars=["module_handle", "input_signal"])
    result = gmodel.fit(
        target_signal,
        input_signal=actual_signal,
        module_handle=pre,
        amplitude=0.0,
        timeconstant=1e-4,
        fit_kws={"epsfcn": 1e-3},
    )  # 'epsfcn' is needed as filter parameters are discretized
    # in precompensationAdvisor module, otherwise fitting will
    # not converge

    print(result.fit_report())
    if do_plot:
        import matplotlib.pyplot as plt

        _, ax = plt.subplots()
        ax.plot(x, result.init_fit, "k", label="initial signal")
        ax.plot(x, result.best_fit, "r", label="fitted signal")
        ax.plot(x, target_signal, "b", label="target signal")
        ax.legend()
        ax.ticklabel_format(axis="both", style="sci", scilimits=(-2, 2))
        ax.set_xlabel("time [s]")
        ax.set_ylabel("Amplitude")
        plt.show()


def generate_target_signal(min_x=-96, max_x=5904, sampling_rate=2.4e9):
    """Returns a step function with given length and sampling interval."""
    x_values = np.array(range(min_x, max_x))
    x_values = [element / sampling_rate for element in x_values]
    signal2 = np.array(np.concatenate((np.zeros(-min_x), np.ones(max_x))))
    return x_values, signal2


def generate_actual_signal(initial_signal, amp=0.4, tau=100e-9, sampling_rate=2.4e9):
    """
    generate "actual signal" through filtering the initial signal with
    an exponential filter and add noise
    """
    from scipy import signal

    # calculate a and b from amplitude and tau
    alpha = 1 - np.exp(-1 / (sampling_rate * tau * (1 + amp)))
    if amp >= 0.0:
        k = amp / (1 + amp - alpha)
        a = [(1 - k + k * alpha), -(1 - k) * (1 - alpha)]
    else:
        k = -amp / (1 + amp) / (1 - alpha)
        a = [(1 + k - k * alpha), -(1 + k) * (1 - alpha)]
    b = [1, -(1 - alpha)]

    distorted_signal = np.array(
        signal.lfilter(b, a, initial_signal) + 0.01 * np.random.normal(size=initial_signal.size)
    )
    return distorted_signal
run_example(device_hd_id)


# ## example_awg_sourcefile

# In[ ]:


"""
Zurich Instruments LabOne Python API Example

Demonstrate how to connect to a Zurich Instruments Arbitrary Waveform Generator
and compile/upload an AWG program to the instrument.
"""

# Copyright 2018 Zurich Instruments AG

import time
import textwrap
import os
import zhinst.utils

# This is only used if this example is ran without the awg_sourcefile
# parameter: To ensure that we have a .seqc source file to use in this example,
# we write this to disk and then compile this file.
SOURCE = textwrap.dedent(
    """// Define an integer constant
    const N = 4096;
    // Create two Gaussian pulses with length N points,
    // amplitude +1.0 (-1.0), center at N/2, and a width of N/8
    wave gauss_pos = 1.0*gauss(N, N/2, N/8);
    wave gauss_neg = -1.0*gauss(N, N/2, N/8);
    // Continuous playback.
    while (true) {
      // Play pulse on AWG channel 1
      playWave(gauss_pos);
      // Wait until waveform playback has ended
      waitWave();
      // Play pulses simultaneously on both AWG channels
      playWave(gauss_pos, gauss_neg);
    }"""
)


def run_example(device_id, awg_sourcefile=None):
    """
    Connect to a Zurich Instruments HDAWG, compile, upload and run an AWG
    sequence program.

    Requirements:

      An HDAWG Instrument.

    Arguments:

      device_id (str): The ID of the device to run the example with. For
        example, `dev8006` or `hdawg-dev8006`.

      awg_sourcefile (str, optional): Specify an AWG sequencer file to compile
        and upload. This file must exist in the AWG source sub-folder of your
        LabOne data directory (this location is provided by the
        directory parameter). The source folder must not be included;
        specify the filename only with extension.

    Raises:

      Exception: AWG functionality is not available.

      RuntimeError: If the device is not "discoverable" from the API.

    See the "LabOne Programming Manual" for further help, available:
      - On Windows via the Start-Menu:
        Programs -> Zurich Instruments -> Documentation
      - On Linux in the LabOne .tar.gz archive in the "Documentation"
        sub-folder.
    """

    # Settings
    apilevel_example = 6  # The API level supported by this example.
    err_msg = "This example can only be ran on either an HDAWG with the AWG option enabled."
    # Call a zhinst utility function that returns:
    # - an API session `daq` in order to communicate with devices via the data server.
    # - the device ID string that specifies the device branch in the server's node hierarchy.
    # - the device's discovery properties.
#     (daq, device, _) = zhinst.utils.create_api_session(
#         device_id, apilevel_example, required_devtype="HDAWG", required_err_msg=err_msg
#     )
    
    zhinst.utils.api_server_version_check(daq)

    # Create a base configuration: Disable all available outputs, awgs, demods, scopes,...
    zhinst.utils.disable_everything(daq, device)

    # 'system/awg/channelgrouping' : Configure how many independent sequencers
    #   should run on the AWG and how the outputs are grouped by sequencer.
    #   0 : 4x2 with HDAWG8; 2x2 with HDAWG4.
    #   1 : 2x4 with HDAWG8; 1x4 with HDAWG4.
    #   2 : 1x8 with HDAWG8.
    # Configure the HDAWG to use one sequencer with the same waveform on all output channels.
    daq.setInt(f"/{device}/system/awg/channelgrouping", 0)

    # Create an instance of the AWG Module
    awgModule = daq.awgModule()
    awgModule.set("device", device)
    awgModule.execute()

    # Get the LabOne user data directory (this is read-only).
    data_dir = awgModule.getString("directory")
    # The AWG Tab in the LabOne UI also uses this directory for AWG seqc files.
    src_dir = os.path.join(data_dir, "awg", "src")
    if not os.path.isdir(src_dir):
        # The data directory is created by the AWG module and should always exist. If this exception is raised,
        # something might be wrong with the file system.
        raise Exception(f"AWG module wave directory {src_dir} does not exist or is not a directory")

    # Note, the AWG source file must be located in the AWG source directory of the user's LabOne data directory.
    if awg_sourcefile is None:
        # Write an AWG source file to disk that we can compile in this example.
        awg_sourcefile = "ziPython_example_awg_sourcefile.seqc"
        with open(os.path.join(src_dir, awg_sourcefile), "w") as f:
            f.write(SOURCE)
    else:
        if not os.path.exists(os.path.join(src_dir, awg_sourcefile)):
            raise Exception(
                f"The file {awg_sourcefile} does not exist, this must be specified via an absolute or relative path."
            )

    print("Will compile and load", awg_sourcefile, "from", src_dir)

    # Transfer the AWG sequence program. Compilation starts automatically.
    awgModule.set("compiler/sourcefile", awg_sourcefile)
    # Note: when using an AWG program from a source file (and only then), the compiler needs to
    # be started explicitly:
    awgModule.set("compiler/start", 1)
    timeout = 20
    t0 = time.time()
    while awgModule.getInt("compiler/status") == -1:
        time.sleep(0.1)
        if time.time() - t0 > timeout:
            Exception("Timeout")

    if awgModule.getInt("compiler/status") == 1:
        # compilation failed, raise an exception
        raise Exception(awgModule.getString("compiler/statusstring"))
    if awgModule.getInt("compiler/status") == 0:
        print("Compilation successful with no warnings, will upload the program to the instrument.")
    if awgModule.getInt("compiler/status") == 2:
        print("Compilation successful with warnings, will upload the program to the instrument.")
        print("Compiler warning: ", awgModule.getString("compiler/statusstring"))

    # Wait for the waveform upload to finish
    time.sleep(0.2)
    i = 0
    while (awgModule.getDouble("progress") < 1.0) and (awgModule.getInt("elf/status") != 1):
        print(f"{i} progress: {awgModule.getDouble('progress'):.2f}")
        time.sleep(0.5)
        i += 1
    print(f"{i} progress: {awgModule.getDouble('progress'):.2f}")
    if awgModule.getInt("elf/status") == 0:
        print("Upload to the instrument successful.")
    if awgModule.getInt("elf/status") == 1:
        raise Exception("Upload to the instrument failed.")

    print("Success. Enabling the AWG.")
    # This is the preferred method of using the AWG: Run in single mode continuous waveform playback is best achieved by
    # using an infinite loop (e.g., while (true)) in the sequencer program.
          
    daq.setInt(f'/{device_hd_id}/sigouts/0/on', 1)     
    daq.setInt(f"/{device}/awgs/0/single", 1)
    daq.setInt(f"/{device}/awgs/0/enable", 1)
run_example(device_hd_id)


# ## example_awg_commandtable

# In[ ]:


"""
Zurich Instruments LabOne Python API Example

Demonstrate how to connect to a Zurich Instruments HDAWG and upload and run an
AWG program using the command table.
"""

# Copyright 2020 Zurich Instruments AG

import os
import time
import json
import textwrap
import numpy as np
import zhinst.utils


def validate_json(json_str):
    """ Validate if string is a valid JSON """
    try:
        json.loads(json_str)
    except ValueError:
        return False

    return True


def run_example(device_id):
    """
    Run the example: Connect to a Zurich Instruments HDAWG upload and run a
    AWG sequence program using the command table, upload a command table,
    and replacing placeholder waveforms in the sequenc program.

    Requirements:

       HDAWG Instrument.

    Arguments:

      device_id (str): The ID of the device to run the example with. For
        example, `dev8006` or `hdawg-dev8006`.

    Returns:

      No return value.

    Raises:

      Exception: If the device is not an HDAWG.

      RuntimeError: If the device is not "discoverable" from the API.

    See the "LabOne Programming Manual" for further help, available:
      - On Windows via the Start-Menu:
        Programs -> Zurich Instruments -> Documentation
      - On Linux in the LabOne .tar.gz archive in the "Documentation"
        sub-folder.
    """

    # Settings
    apilevel_example = 6  # The API level supported by this example.
    err_msg = "This example can only be ran on an HDAWG."
    # Call a zhinst utility function that returns:
    # - an API session `daq` in order to communicate with devices via the data server.
    # - the device ID string that specifies the device branch in the server's node hierarchy.
    # - the device's discovery properties.
#     (daq, device, _) = zhinst.utils.create_api_session(
#         device_id, apilevel_example, required_devtype="HDAWG", required_err_msg=err_msg
#     )

#     daq = zp.ziDAQServer(host, 8004, 6)
#     daq.connectDevice(device_id, '1gbe')
#     zhinst.utils.api_server_version_check(daq)
    zhinst.utils.api_server_version_check(daq)

    # Create a base configuration: Disable all available outputs, awgs, demods, scopes,...
    zhinst.utils.disable_everything(daq, device)

    # 'system/awg/channelgrouping' : Configure how many independent sequencers
    #   should run on the AWG and how the outputs are grouped by sequencer.
    #   0 : 4x2 with HDAWG8; 2x2 with HDAWG4.
    #   1 : 2x4 with HDAWG8; 1x4 with HDAWG4.
    #   2 : 1x8 with HDAWG8.
    # Configure the HDAWG to use one sequencer with the same waveform on all output channels.
    daq.setInt(f"/{device}/system/awg/channelgrouping", 0)

    # Some basic device configuration to output the generated wave on Wave outputs 1 and 2.
    amplitude = 1.0

    exp_setting = [
        [f"/{device}/sigouts/0/on", 1],
        [f"/{device}/sigouts/1/on", 1],
        [f"/{device}/sigouts/0/range", 1],
        [f"/{device}/sigouts/1/range", 1],
        [f"/{device}/awgs/0/outputs/0/amplitude", amplitude],
        [f"/{device}/awgs/0/outputs/1/amplitude", amplitude],
        [f"/{device}/awgs/0/outputs/0/modulation/mode", 0],
        [f"/{device}/awgs/0/time", 0],
        [f"/{device}/awgs/*/enable", 0],
        [f"/{device}/awgs/0/userregs/0", 0],
    ]
    daq.set(exp_setting)
    # Ensure that all settings have taken effect on the device before continuing.
    daq.sync()

    # Number of points in AWG waveform
    AWG_N = 2000

    # Define an AWG program as a string stored in the variable awg_program, equivalent to what would
    # be entered in the Sequence Editor window in the graphical UI. Different to a self-contained program,
    # this example refers to a command table by the instruction "executeTableEntry", and to a placeholder
    # waveform p by the instruction "placeholder". Both the command table and the waveform data for the
    # waveform p need to be uploaded separately before this sequence program can be run.

    awg_program = textwrap.dedent(
        """\
        // Define placeholder with 1024 samples:
        wave p = placeholder(1024);

        // Assign placeholder to waveform index 10
        assignWaveIndex(p, p, 10);

        while(true) {
          executeTableEntry(0);
        }
        """
    )

    # JSON string specifying the command table to be uploadedCommand Table JSON
    json_str = textwrap.dedent(
        """
        {
          "$schema": "http://docs.zhinst.com/hdawg/commandtable/v2/schema",
          "header": {
            "version": "0.2",
            "partial": false
          },
          "table": [
            {
              "index": 0,
              "waveform": {
                "index": 10
              },
              "amplitude0": {
                "value": 1.0
              },
              "amplitude1": {
                "value": 1.0
              }
            }
          ]
        }
        """
    )

    # Ensure command table is valid (valid JSON and compliant with schema)
    assert validate_json(json_str)

    # Create an instance of the AWG Module
    awgModule = daq.awgModule()
    awgModule.set("device", device)
    awgModule.execute()

    # Get the modules data directory
    data_dir = awgModule.getString("directory")
    # All CSV files within the waves directory are automatically recognized by the AWG module
    wave_dir = os.path.join(data_dir, "awg", "waves")
    if not os.path.isdir(wave_dir):
        # The data directory is created by the AWG module and should always exist. If this exception is raised,
        # something might be wrong with the file system.
        raise Exception(
            f"AWG module wave directory {wave_dir} does not exist or is not a directory"
        )

    # Transfer the AWG sequence program. Compilation starts automatically.
    awgModule.set("compiler/sourcestring", awg_program)
    # Note: when using an AWG program from a source file (and only then), the compiler needs to
    # be started explicitly with awgModule.set('compiler/start', 1)
    while awgModule.getInt("compiler/status") == -1:
        time.sleep(0.1)

    if awgModule.getInt("compiler/status") == 1:
        # compilation failed, raise an exception
        raise Exception(awgModule.getString("compiler/statusstring"))

    if awgModule.getInt("compiler/status") == 0:
        print(
            "Compilation successful with no warnings, will upload the program to the instrument."
        )
    if awgModule.getInt("compiler/status") == 2:
        print(
            "Compilation successful with warnings, will upload the program to the instrument."
        )
        print("Compiler warning: ", awgModule.getString("compiler/statusstring"))

    # Wait for the waveform upload to finish
    time.sleep(0.2)
    i = 0
    while (awgModule.getDouble("progress") < 1.0) and (
        awgModule.getInt("elf/status") != 1
    ):
        print(f"{i} progress: {awgModule.getDouble('progress'):.2f}")
        time.sleep(0.2)
        i += 1
    print(f"{i} progress: {awgModule.getDouble('progress'):.2f}")
    if awgModule.getInt("elf/status") == 0:
        print("Upload to the instrument successful.")
    if awgModule.getInt("elf/status") == 1:
        raise Exception("Upload to the instrument failed.")

    # Upload command table
    daq.setVector(f"/{device}/awgs/0/commandtable/data", json_str)

    # Replace the placeholder waveform with a new one with a Gaussian shape
    x = np.linspace(0, 1024, 1024)
    x_center = 512
    sigma = 150
    w = np.array(
        np.exp(-np.power(x - x_center, 2.0) / (2 * np.power(sigma, 2.0))), dtype=float
    )
    waveform_native = zhinst.utils.convert_awg_waveform(w, -w)

    # Upload waveform data with the index 10 (this is the index assigned with the assignWaveIndex
    # sequencer instruction
    index = 10
    path = f"/{device:s}/awgs/0/waveform/waves/{index:d}"
    daq.setVector(path, waveform_native)

    print(f"Enabling the AWG.")
    # This is the preferred method of using the AWG: Run in single mode continuous waveform playback
    # is best achieved by using an infinite loop (e.g., while (true)) in the sequencer program.
    daq.setInt(f"/{device}/awgs/0/single", 1)
    daq.setInt(f"/{device}/awgs/0/enable", 1)
run_example(device_hd_id)


# ## Replace waveforms in 1x8 mode for HDAWG8

# In[ ]:


awg_program = textwrap.dedent(
'''
wave w1 = placeholder(20000, true);// true means using marker, default is false, available with labOne 20.07
wave w2 = placeholder(20000, true);
wave w3 = placeholder(20000, true);
wave w4 = placeholder(20000, true);
wave w5 = placeholder(20000, true);
wave w6 = placeholder(20000, true);
wave w7 = placeholder(20000, true);
wave w8 = placeholder(20000, true);

while(true){
playWave(1, w1, 2, w2, 3, w3, 4, w4, 5, w5, 6, w6, 7, w7, 8, w8);
waitWave();
wait(87482);
}

''')
daq.setInt(f'{device_hd_id}/system/awg/channelgrouping', 2)
for i in range(2):
    daq.setInt(f'/{device_hd_id}/sigouts/{i}/on', 1)
    daq.setInt(f'/{device_hd_id}/sigouts/{i}/on', 1)
    
create_and_compile_awg(daq, device_hd_id, awg_program, 0, 5) # compile the seqc codes, timeout here is 10s.

# this is the configuration for grouping 1x8
N = 20000
w_1 = 0.5 * np.ones(N, dtype=float)   # wave 1 in physical output 1 of 1 AWG core
w_2 = 0.75 * np.ones(N, dtype=float)  # wave 2 in physical output 1 of 1 AWG core

marker = np.zeros(N)              
marker[3200: N-3200] = 1              # Marker setting for 1 AWG core

w_new1 = zu.convert_awg_waveform(wave1=w_1, wave2= w_2, markers=marker)

for i in range(4): 
    # here '.../waves/0' in fact for both channel 0 and 1, because the w_new1 includes 2 waveforms
    daq.setVector(f'/{device_hd_id:s}/awgs/{i}/waveform/waves/0', w_new1) 
    
        
daq.setInt(f'{device_hd_id}/awgs/0/enable', 1)
print('enable hd_awg')
time.sleep(20)
daq.setInt(f'{device_hd_id}/awgs/0/enable', 0)
for i in range(2):
    daq.setInt(f'/{device_hd_id}/sigouts/{i}/on', 0)
    daq.setInt(f'/{device_hd_id}/sigouts/{i}/on', 0)
print('disable hd_awg')

