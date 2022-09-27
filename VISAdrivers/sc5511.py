# -*- coding: utf-8 -*-

import os
import ctypes
import logging
from typing import Any, Dict, Optional

import InstrumentDriver

class device_rf_params_t(ctypes.Structure):
    _fields_ = [("rf1_freq", ctypes.c_ulonglong),
                ("start_freq", ctypes.c_ulonglong),
                ("stop_freq", ctypes.c_ulonglong),
                ("step_freq", ctypes.c_ulonglong),
                ("sweep_dwell_time", ctypes.c_uint),
                ("sweep_cycles", ctypes.c_uint),
                ("buffer_time", ctypes.c_uint),
                ("rf_level", ctypes.c_float),
                ("rf2_freq", ctypes.c_short)
                ]


class Device_temperature_t(ctypes.Structure):
    _fields_ = [("device_temp", ctypes.c_float)]


class operate_status_t(ctypes.Structure):
    _fields_ = [("rf1_lock_mode", ctypes.c_ubyte),
                ("rf1_loop_gain", ctypes.c_ubyte),
                ("device_access", ctypes.c_ubyte),
                ("rf2_standby", ctypes.c_ubyte),
                ("rf1_standby", ctypes.c_ubyte),
                ("auto_pwr_disable", ctypes.c_ubyte),
                ("alc_mode", ctypes.c_ubyte),
                ("rf1_out_enable", ctypes.c_ubyte),
                ("ext_ref_lock_enable", ctypes.c_ubyte),
                ("ext_ref_detect", ctypes.c_ubyte),
                ("ref_out_select", ctypes.c_ubyte),
                ("list_mode_running", ctypes.c_ubyte),
                ("rf1_mode", ctypes.c_ubyte),
                ("harmonic_ss", ctypes.c_ubyte),
                ("over_temp", ctypes.c_ubyte)
                ]


class pll_status_t(ctypes.Structure):
    _fields_ = [("sum_pll_ld", ctypes.c_ubyte),
                ("crs_pll_ld", ctypes.c_ubyte),
                ("fine_pll_ld", ctypes.c_ubyte),
                ("crs_ref_pll_ld", ctypes.c_ubyte),
                ("crs_aux_pll_ld", ctypes.c_ubyte),
                ("ref_100_pll_ld", ctypes.c_ubyte),
                ("ref_10_pll_ld", ctypes.c_ubyte),
                ("rf2_pll_ld", ctypes.c_ubyte)]


class list_mode_t(ctypes.Structure):
    _fields_ = [("sss_mode", ctypes.c_ubyte),
                ("sweep_dir", ctypes.c_ubyte),
                ("tri_waveform", ctypes.c_ubyte),
                ("hw_trigger", ctypes.c_ubyte),
                ("step_on_hw_trig", ctypes.c_ubyte),
                ("return_to_start", ctypes.c_ubyte),
                ("trig_out_enable", ctypes.c_ubyte),
                ("trig_out_on_cycle", ctypes.c_ubyte)]


class device_status_t(ctypes.Structure):
    _fields_ = [("list_mode", list_mode_t),
                ("operate_status_t", operate_status_t),
                ("pll_status_t", pll_status_t)]

class device_info_t(ctypes.Structure):
    _fields_ = [("serial_number", ctypes.c_uint32),
                ("hardware_revision", ctypes.c_float),
                ("firmware_revision", ctypes.c_float),
                ("manufacture_date", ctypes.c_uint32)
                ]

# End of Structures------------------------------------------------------------  
'''             
class SignalCore_SC5511A(Instrument):
    def get_idn(self) -> Dict[str, Optional[str]]:
        logging.info(__name__ + " : Getting device info")
        self._handle = ctypes.c_void_p(self._dll.sc5511a_open_device(self._serial_number))
        self._dll.sc5511a_get_device_info(self._handle, ctypes.byref(self._device_info))
        device_info = self._device_info
        self._dll.sc5511a_close_device(self._handle)
        def date_decode(date_int:int):
            date_str = f"{date_int:032b}"
            yr = f"20{int(date_str[:8],2)}"
            month = f"{int(date_str[16:24],2)}"
            day = f"{int(date_str[8:16],2)}"
            return f"{month}/{day}/{yr}"
        IDN: Dict[str, Optional[str]] = {
            'vendor': "SignalCore",
            'model': "SC5511A",
            'serial_number': self._serial_number.value.decode("utf-8"),
            'firmware_revision': device_info.firmware_revision,
            'hardware_revision': device_info.hardware_revision,
            'manufacture_date': date_decode(device_info.manufacture_date)
            }
        return IDN
        '''

class Driver(InstrumentDriver.InstrumentWorker):
    '''
    def __init__(self, name: str, serial_number: str, dll = None, debug = False, **kwargs: Any):
        super().__init__(name, **kwargs)
        logging.info(__name__ + f' : Initializing instrument SignalCore generator {serial_number}')
        if dll is not None:
            self._dll = dll
        else:
            # this may not work
            self._dll = ctypes.CDLL('.//sc5511a.dll')

        if debug:
            print(self._dll)

        self._dll.sc5511a_open_device.restype = ctypes.c_uint64
        self._handle = ctypes.c_void_p(
            self._dll.sc5511a_open_device(ctypes.c_char_p(bytes(serial_number, 'utf-8'))))
        self._serial_number = ctypes.c_char_p(bytes(serial_number, 'utf-8'))
        self._open = False
        self._temperature = Device_temperature_t(0)

        if debug:
            print(serial_number, self._handle)
            self._dll.sc5511a_get_device_status(self._handle, ctypes.byref(self._device_status))
            status = self._device_status.operate_status_t.rf1_out_enable
            print('check status', status)

        self._dll.sc5511a_close_device(self._handle)
        self._device_info = device_info_t(0, 0, 0, 0)
        self.get_idn()
        self.do_set_auto_level_disable(0) # setting this to 1 will lead to unstable output power


        if self._device_status.operate_status_t.ext_ref_lock_enable == 0:
            self.do_set_reference_source(1)
        '''

    def checkIfOpen(self):
        if not self._open:
            self._dll.sc5511a_open_device.restype = ctypes.c_uint64
            self._handle_addr = ctypes.c_void_p(self._dll.sc5511a_open_device(self._serial_number))
            self._handle = ctypes.c_void_p(self._handle_addr.value)
            self._open = True  
            self._close = True

    def closeWhenOpened(self):
        if self._close:
            self._dll.sc5511a_close_device(self._handle)
            self._close = False
            self._open = False


    def performOpen(self, options={}):
        dirname = os.path.dirname(__file__)
        dll_path = os.path.join(dirname,'sc5511a.dll')
        os.chdir(dirname)
        self._dll = ctypes.CDLL(dll_path)

        self._serial_number = ctypes.c_char_p(bytes(self.getName(), 'utf-8'))
        self._rf_params = device_rf_params_t(0, 0, 0, 0, 0, 0, 0, 0, 0)
        self._status = operate_status_t(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        self._pll_status = pll_status_t(0, 0, 0, 0, 0, 0, 0, 0)
        self._list_mode = list_mode_t(0, 0, 0, 0, 0, 0, 0, 0)
        self._device_status = device_status_t(self._list_mode, self._status, self._pll_status)
        self._temperature = Device_temperature_t(0)
        self._dll.sc5511a_open_device.restype = ctypes.c_uint64
        self._handle_addr = ctypes.c_void_p(self._dll.sc5511a_open_device(self._serial_number))
        self._handle = ctypes.c_void_p(self._handle_addr.value)
        self._open = True
        self._close = False
        return self._handle

    def performClose(self, bError=False, options={}):
        self._dll.sc5511a_close_device(self._handle)
        self._open = False
        self._close = False
        return

    def performSetValue(self, quant, value, sweepRate=0.0, options={}):
        self.checkIfOpen()
        returnVal = value
        if quant.name == 'Frequency':
            """
            Sets RF1 frequency. Valid between 100MHz and 20GHz
                Args:
                    frequency (int) = frequency in Hz
            """
            frequency = value
            c_freq = ctypes.c_ulonglong(int(frequency))
            logging.info(__name__ + ' : Setting frequency to %s' % frequency)
            successVal = self._dll.sc5511a_set_freq(self._handle, c_freq)
        elif quant.name == 'Power':
            power = value
            logging.info(__name__ + ' : Setting power to %s' % power)
            c_power = ctypes.c_float(power)
            successVal = self._dll.sc5511a_set_level(self._handle, c_power)
        elif quant.name == 'Reference source':
            lock_to_external = value
            logging.info(__name__ + ' : Setting reference source to %s' % lock_to_external)
            high = ctypes.c_ubyte(0)
            lock = ctypes.c_ubyte(lock_to_external)
            returnVal = self._dll.sc5511a_set_clock_reference(self._handle, high, lock)
            successVal = True
        elif quant.name == 'Output status':
            """
            Turns the output of RF1 on or off.
                Input:
                    enable (int) = OFF = 0 ; ON = 1
            """
            enable = value
            logging.info(__name__ + ' : Setting output to %s' % enable)
            c_enable = ctypes.c_ubyte(enable)
            successVal = self._dll.sc5511a_set_output(self._handle, c_enable) 
        elif quant.name == 'Auto level disable':
            enable = value
            logging.info(__name__ + ' : Settingalc auto to %s' % enable)
            if enable:
                enable = 0
            else:
                enable = 1
            c_enable = ctypes.c_ubyte(enable)
            successVal = self._dll.sc5511a_set_auto_level_disable(self._handle, c_enable)
        self.closeWhenOpened()
        return returnVal
        

    def performGetValue(self, quant, options={}):
        self.checkIfOpen()
        returnVal = None
        self._dll.sc5511a_get_temperature(self._handle, ctypes.byref(self._temperature))
        self._dll.sc5511a_get_rf_parameters(self._handle, ctypes.byref(self._rf_params))
        self._dll.sc5511a_get_device_status(self._handle, ctypes.byref(self._device_status))
        self._status = self._device_status.operate_status_t
        # rf1 frequency
        if quant.name == 'Frequency':
            logging.info(__name__ + ' : Getting frequency')
            returnVal = self._rf_params.rf1_freq
        # rf1 power
        elif quant.name == 'Power':
            logging.info(__name__ + ' : Getting Power')
            returnVal = self._rf_params.rf_level
        elif quant.name == 'Output status':
            '''
            Reads the output status of RF1
                Output:
                    status (int) : OFF = 0 ; ON = 1
            '''
            logging.info(__name__ + ' : Getting output')
            returnVal = self._status.rf1_out_enable
        elif quant.name == 'Reference source':
            logging.info(__name__ + ' : Getting reference source')
            returnVal = self._status.ext_ref_lock_enable
        elif quant.name == 'Auto level disable':
            logging.info(__name__ + ' : Getting alc auto status')
            enabled = self._status.auto_pwr_disable
            returnVal = not enabled
        # temp of device
        elif quant.name == 'Temperature':
            logging.info(__name__ + " : Getting device temperature")
            returnVal = self._temperature.device_temp
        self.closeWhenOpened()
        return returnVal

