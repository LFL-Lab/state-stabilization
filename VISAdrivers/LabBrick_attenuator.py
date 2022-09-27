#!/usr/bin/env python
'''
    Python based USB CLI Application interface control for vaunix LDA devices on windows-64bit machine
    JA 05/05/2021    Initial Verison of USB CLI control Interface
    (c) 2021-2022 by Vaunix Technology Corporation, all rights reserved
'''
from ctypes import *
import ctypes
import os
import time
import sys, getopt
from time import sleep

from ctypes import c_int, c_uint, c_bool, c_float, POINTER, byref

class Error(Exception):
    pass

# open dll
# sPath = os.path.dirname(os.path.abspath(__file__))
# _lib = ctypes.CDLL(os.path.join(sPath, 'VNX_atten64'))

# # define data types used by this dll
# DEVID = c_uint
# LVSTATUS = c_uint
# STRING = ctypes.c_char_p
# ACTIVEDEVICES = DEVID*64

# def getDllObject(sName, argtypes=[DEVID,], restype=LVSTATUS):
#     """Create a dll ojbect with input and output types"""
#     obj = getattr(_lib, sName)
#     obj.restype = restype
#     obj.argypes = argtypes
#     return obj


# #VNX_FSYNSTH_API int fnLMS_InitDevice(DEVID deviceID);
# #VNX_FSYNSTH_API int fnLMS_CloseDevice(DEVID deviceID);
# fnLMS_InitDevice = getDllObject('fnLMS_InitDevice')
# fnLMS_CloseDevice = getDllObject('fnLMS_CloseDevice')

# #VNX_FSYNSTH_API void fnLMS_SetTestMode(bool testmode);
# #VNX_FSYNSTH_API int fnLMS_GetNumDevices();
# #VNX_FSYNSTH_API int fnLMS_GetDevInfo(DEVID *ActiveDevices);
# #VNX_FSYNSTH_API int fnLMS_GetModelName(DEVID deviceID, char *ModelName);
# #VNX_FSYNSTH_API int fnLMS_GetSerialNumber(DEVID deviceID);
# fnLMS_SetTestMode = getDllObject('fnLMS_SetTestMode', False, None)
# fnLMS_GetNumDevices = getDllObject('fnLMS_GetNumDevices', [], restype=c_int)
# fnLMS_GetDevInfo = getDllObject('fnLMS_GetDevInfo', [POINTER(ACTIVEDEVICES)], restype=c_int)
# fnLMS_GetModelName = getDllObject('fnLMS_GetModelNameA', [DEVID, STRING], restype=c_int)
# fnLMS_GetSerialNumber = getDllObject('fnLMS_GetSerialNumber', restype=c_int)

# #VNX_FSYNSTH_API LVSTATUS fnLMS_SetFrequency(DEVID deviceID, int frequency);
# #VNX_FSYNSTH_API LVSTATUS fnLMS_SetPowerLevel(DEVID deviceID, int powerlevel);
# #VNX_FSYNSTH_API LVSTATUS fnLMS_SetRFOn(DEVID deviceID, bool on);
# #VNX_FSYNSTH_API LVSTATUS fnLMS_SetUseInternalRef(DEVID deviceID, bool internal);
# #VNX_FSYNSTH_API LVSTATUS fnLMS_SetUseExternalPulseMod(DEVID deviceID, bool external);
# #VNX_FSYNSTH_API LVSTATUS fnLMS_SetFastPulsedOutput(DEVID deviceID, float pulseontime, float pulsereptime, bool on);
# fnLMS_SetFrequency = getDllObject('fnLMS_SetFrequency', [DEVID, c_int])
# fnLMS_SetPowerLevel = getDllObject('fnLMS_SetPowerLevel', [DEVID, c_int])
# fnLMS_SetRFOn = getDllObject('fnLMS_SetRFOn', [DEVID, c_bool])
# fnLMS_SetUseInternalRef = getDllObject('fnLMS_SetUseInternalRef', [DEVID, c_bool])
# fnLMS_SetUseExternalPulseMod = getDllObject('fnLMS_SetUseExternalPulseMod', [DEVID, c_bool])
# fnLMS_SetFastPulsedOutput = getDllObject('fnLMS_SetFastPulsedOutput', [DEVID, c_float, c_float, c_bool])


# # VNX_FSYNSTH_API int fnLMS_GetFrequency(DEVID deviceID);
# # VNX_FSYNSTH_API int fnLMS_GetPowerLevel(DEVID deviceID);
# # VNX_FSYNSTH_API int fnLMS_GetRF_On(DEVID deviceID);
# # VNX_FSYNSTH_API int fnLMS_GetUseInternalRef(DEVID deviceID);
# fnLMS_GetFrequency = getDllObject('fnLMS_GetFrequency', restype=c_int)
# fnLMS_GetPowerLevel = getDllObject('fnLMS_GetPowerLevel', restype=c_int)
# fnLMS_GetRF_On = getDllObject('fnLMS_GetRF_On', restype=c_int)
# fnLMS_GetUseInternalRef = getDllObject('fnLMS_GetUseInternalRef', restype=c_int)

# #VNX_FSYNSTH_API float fnLMS_GetPulseOnTime(DEVID deviceID);
# #VNX_FSYNSTH_API float fnLMS_GetPulseOffTime(DEVID deviceID);
# #VNX_FSYNSTH_API int fnLMS_GetPulseMode(DEVID deviceID);
# #VNX_FSYNSTH_API int fnLMS_GetUseInternalPulseMod(DEVID deviceID);
# fnLMS_GetPulseOnTime = getDllObject('fnLMS_GetPulseOnTime', restype=c_float)
# fnLMS_GetPulseOffTime = getDllObject('fnLMS_GetPulseOffTime', restype=c_float)
# fnLMS_GetPulseMode = getDllObject('fnLMS_GetPulseMode', restype=c_int)
# fnLMS_GetUseInternalPulseMod = getDllObject('fnLMS_GetUseInternalPulseMod', restype=c_int)

# #VNX_FSYNSTH_API int fnLMS_GetMaxPwr(DEVID deviceID);
# #VNX_FSYNSTH_API int fnLMS_GetMinPwr(DEVID deviceID);
# #VNX_FSYNSTH_API int fnLMS_GetMaxFreq(DEVID deviceID);
# #VNX_FSYNSTH_API int fnLMS_GetMinFreq(DEVID deviceID);
# fnLMS_GetMaxPwr = getDllObject('fnLMS_GetMaxPwr', restype=c_int)
# fnLMS_GetMinPwr = getDllObject('fnLMS_GetMinPwr', restype=c_int)
# fnLMS_GetMaxFreq = getDllObject('fnLMS_GetMaxFreq', restype=c_int)
# fnLMS_GetMinFreq = getDllObject('fnLMS_GetMinFreq', restype=c_int)

class LabBrick_attn():
    def __init__(self, name="ldaApi", port=''):
        self.vnx = cdll.VNX_atten64
        self.vnx.fnLDA_SetTestMode(False)
        self.devices_num = self.get_devices_number()
        # print(self.devices_num)
        DeviceIDArray = c_int * self.devices_num
        self.devices_list = DeviceIDArray()
        # print(self.devices_list)
        # fill the array with the ID's of connected attenuators
        self.vnx.fnLDA_GetDevInfo(self.devices_list)

    # def initDevice(self, serial):
    #     # get list of devices
    #     try:
    #         iSerial = int(serial)
    #     except:
    #         iSerial = 0
    #     lDev = self.get_devices_list()
    #     # lSerial = [d['serial'] for d in lDev]
    #     # if iSerial not in lSerial:
    #     #     # raise error if device not found
    #     #     sErr = ('Device with serial number "%d" cannot be found.\n\n' % iSerial) +\
    #     #             'Devices detected:\n'
    #     #     for dDev in lDev:
    #     #         sErr += ('Name: %s, Serial: %d\n' % (dDev['name'], dDev['serial']))
    #     #     raise Error(sErr)
    #     indx = lSerial.index(iSerial)
    #     self.device_id = DEVID(lDev[indx]['device_id'])
    #     status = fnLMS_InitDevice(self.device_id)
    #     self.maxPower = float(fnLMS_GetMaxPwr(self.device_id))*0.25
    #     self.minPower = float(fnLMS_GetMinPwr(self.device_id))*0.25
    #     self.maxFreq = float(fnLMS_GetMaxFreq(self.device_id))*10.
    #     self.minFreq = float(fnLMS_GetMinFreq(self.device_id))*10.
    #     self.check_error(status)



    def get_devices_list(self):
        print(self.devices_list)
        return self.devices_list

    def get_devices_number(self):
        return self.vnx.fnLDA_GetNumDevices()

    def get_serial_number(self, devid):
        id = int(devid)
        return self.vnx.fnLDA_GetSerialNumber(id)

    def check_deviceexists(self, devid):
        if((devid <= self.devices_num) and (self.devices_num >0)):
            return 1
        else:
            return 0

    # Get Model Name
    def get_modelname(self, devid):
        self.data = ctypes.create_string_buffer(32)
        self.vnx.fnLDA_GetModelNameA(devid, ctypes.byref(self.data))
        return self.data.value

    # Ger Serial Number
    def get_serial_number(self, devid):
        self.data= self.vnx.fnLDA_GetSerialNumber(devid)
        return self.data

    # Get Software Version
    def get_swversion(self):
        self.data = self.vnx.fnLDA_GetDLLVersion()
        return self.data

    # Get Min.Frequency
    def get_minfrequency(self, devid):
        self.data= self.vnx.fnLDA_GetMinWorkingFrequency(devid)
        return (self.data/10)

    # Get Max. Frequency
    def get_maxfrequency(self, devid):
        self.data= self.vnx.fnLDA_GetMaxWorkingFrequency(devid)
        return (self.data/10)

    # Get Min. Attenuation
    def get_minattenuation(self, devid):
        self.data= self.vnx.fnLDA_GetMinAttenuation(devid)
        return (self.data/4.0)

    # Get Max. Attenuation
    def get_maxattenuation(self, devid):
        self.data= self.vnx.fnLDA_GetMaxAttenuation(devid)
        return (self.data/4.0)

    # Get Current Frequency
    def get_currentfrequency(self, devid):
        self.data= self.vnx.fnLDA_GetWorkingFrequency(devid)
        return (self.data/10)

    # Get Current Attenuation
    def get_currentattenuation(self, devid):
        self.data= self.vnx.fnLDA_GetAttenuationHR(devid)
        return (self.data/20.0)

    # Get Ramp Start
    def get_rampstart(self, devid):
        self.data= self.vnx.fnLDA_GetRampStartHR(devid)
        return (self.data/20.0)

    # Get Ramp End
    def get_rampend(self, devid):
        self.data= self.vnx.fnLDA_GetRampEndHR(devid)
        return (self.data/20.0)

    # Get Dwell time
    def get_dwelltime(self, devid):
        self.data= self.vnx.fnLDA_GetDwellTime(devid)
        return self.data

    # Get bi-directional ramp dwelltime
    def get_bidirectional_dwelltime(self, devid):
        self.data= self.vnx.fnLDA_GetDwellTimeTwo(devid)
        return self.data

    # Get Idle time
    def get_idletime(self, devid):
        self.data= self.vnx.fnLDA_GetIdleTime(devid)
        return self.data

    # Get hold time
    def get_holdtime(self, devid):
        self.data= self.vnx.fnLDA_GetHoldTime(devid)
        return self.data

    # Get profile count
    def get_profilecount(self, devid):
        self.data= self.vnx.fnLDA_GetProfileCount(devid)
        return self.data

    # Get Profile Max length
    def get_profilemaxlength(self, devid):
        self.data= self.vnx.fnLDA_GetProfileMaxLength(devid)
        return self.data

    # Get Profile Dwell time
    def get_profiledwelltime(self, devid):
        self.data= self.vnx.fnLDA_GetProfileDwellTime(devid)
        return self.data

    # Get Profile Idle time
    def get_profileidletime(self, devid):
        self.data= self.vnx.fnLDA_GetProfileIdleTime(devid)
        return self.data

    # Get Profile Index
    def get_profileindex(self, devid):
        self.data= self.vnx.fnLDA_GetProfileIndex(devid)
        return self.data

    # Set Attenuation
    def set_attenuation(self, devid, attenuation):
        try:
            self.vnx.fnLDA_InitDevice(devid)
            self.vnx.fnLDA_SetAttenuationHR(devid, int(attenuation)*20)
            # self.vnx.fnLDA_CloseDevice(devid)
        except:
            print("An exception occurred")
        return True

    # Set Frequency
    def set_frequency(self, devid, frequency):
        try:
            self.vnx.fnLDA_InitDevice(devid)
            self.vnx.fnLDA_SetWorkingFrequency(devid, int(frequency)*10)
            self.vnx.fnLDA_CloseDevice(devid)
        except:
            print("An exception occurred")
        return True

    # Set Channel
    def set_channel(self, devid, channel):
        try:
            self.vnx.fnLDA_InitDevice(devid)
            self.vnx.fnLDA_SetChannel(devid, int(channel))
            self.vnx.fnLDA_CloseDevice(devid)
        except:
            print("An exception occurred")
        return True

    # Set Ramp Start
    def set_rampstartattenuation(self, devid, attenuation):
        try:
            self.vnx.fnLDA_InitDevice(devid)
            self.vnx.fnLDA_SetRampStartHR(devid, int(attenuation)*20)
            self.vnx.fnLDA_CloseDevice(devid)
        except:
            print("An exception occurred")
        return True

    # Set Ramp End
    def set_rampendattenuation(self, devid, attenuation):
        try:
            self.vnx.fnLDA_InitDevice(devid)
            self.vnx.fnLDA_SetRampEndHR(devid, int(attenuation)*20)
            self.vnx.fnLDA_CloseDevice(devid)
        except:
            print("An exception occurred")
        return True

    # Set dwell time
    def set_dwelltime(self, devid, dwelltime):
        try:
            self.vnx.fnLDA_InitDevice(devid)
            self.vnx.fnLDA_SetDwellTime(devid, int(dwelltime))
            self.vnx.fnLDA_CloseDevice(devid)
        except:
            print("An exception occurred")
        return True

    # Set idle time
    def set_idletime(self, devid, idletime):
        try:
            self.vnx.fnLDA_InitDevice(devid)
            self.vnx.fnLDA_SetIdleTime(devid, int(idletime))
            self.vnx.fnLDA_CloseDevice(devid)
        except:
            print("An exception occurred")
        return True

    # Set bidirectional dwell time
    def set_bidirectionaldwelltime(self, devid, dwelltime):
        try:
            self.vnx.fnLDA_InitDevice(devid)
            self.vnx.fnLDA_SetDwellTimeTwo(devid, int(dwelltime))
            self.vnx.fnLDA_CloseDevice(devid)
        except:
            print("An exception occurred")
        return True

    # Set hold time
    def set_holdtime(self, devid, holdtime):
        try:
            self.vnx.fnLDA_InitDevice(devid)
            self.vnx.fnLDA_SetHoldTime(devid, int(holdtime))
            self.vnx.fnLDA_CloseDevice(devid)
        except:
            print("An exception occurred")
        return True


    # Set Ramp direction  -- True: Upwards, False: Downwards
    def set_rampdirection(self, devid, rampdirection):
        try:
            self.vnx.fnLDA_InitDevice(devid)
            self.vnx.fnLDA_SetRampDirection(devid, int(rampdirection))
            self.vnx.fnLDA_CloseDevice(devid)
        except:
            print("An exception occurred")
        return True

    # Set bidirectional Ramp direction  -- True: Bi-directional, False: Uni-directional
    def set_rampbidirectional(self, devid, rampdirection):
        try:
            self.vnx.fnLDA_InitDevice(devid)
            self.vnx.fnLDA_SetRampBidirectional(devid, int(rampdirection))
            self.vnx.fnLDA_CloseDevice(devid)
        except:
            print("An exception occurred")
        return True

    # Set Rampmode -- True - Continuous, False - Once
    def set_rampmode(self, devid, rampmode):
        try:
            self.vnx.fnLDA_InitDevice(devid)
            self.vnx.fnLDA_SetRampMode(devid, int(rampmode))
            self.vnx.fnLDA_CloseDevice(devid)
        except:
            print("An exception occurred")
        return True

    # Set profile element data
    def set_profilelements(self, devid, profileindex, profiledata):
        try:
            self.vnx.fnLDA_InitDevice(devid)
#            print('profiledata:%r', int(profiledata)*10)
            self.vnx.fnLDA_SetProfileElement(devid, int(profileindex), int(profiledata)*20)
            self.vnx.fnLDA_CloseDevice(devid)
        except:
            print("An exception occurred")
        return True

    # Set profile count
    def set_profilecount(self, devid, profilelen):
        try:
            self.vnx.fnLDA_InitDevice(devid)
            self.vnx.fnLDA_SetProfileCount(devid, int(profilelen))
            self.vnx.fnLDA_CloseDevice(devid)
        except:
            print("An exception occurred")
        return True

    # Set profile Idletime
    def set_profileidletime(self, devid, idletime):
        try:
            self.vnx.fnLDA_InitDevice(devid)
            self.vnx.fnLDA_SetProfileIdleTime(devid, int(idletime))
            self.vnx.fnLDA_CloseDevice(devid)
        except:
            print("An exception occurred")
        return True

    # Set profile Dwell time
    def set_profiledwelltime(self, devid, dwelltime):
        try:
            self.vnx.fnLDA_InitDevice(devid)
            self.vnx.fnLDA_SetProfileDwellTime(devid, int(dwelltime))
            self.vnx.fnLDA_CloseDevice(devid)
        except:
            print("An exception occurred")
        return True

    # Set profile mode  0 - Off, 1 - Profile Once, 2 - Repeat
    def set_profilemode(self, devid, profilemode):
        try:
            self.vnx.fnLDA_InitDevice(devid)
            self.vnx.fnLDA_StartProfile(devid, int(profilemode))
            self.vnx.fnLDA_CloseDevice(devid)
        except:
            print("An exception occurred")
        return True

    def set_savesettings(self, devid):
        try:
            self.vnx.fnLDA_InitDevice(devid)
            self.vnx.fnLDA_SaveSettings(devid)
            self.vnx.fnLDA_CloseDevice(devid)
        except:
            print("An exception occurred")
        return True

    # Main Function call
    def main(self, argv):
        devid = '26551'
        attenuation = '10'
        frequency = ''
        attobj = LabBrick_attn(devid)
        try:
            opts, args = getopt.getopt(argv,"hi:a:f:c:s:e:w:d:o:b:D:M:B:C:I:W:O:S:F:r",["idevid=","aattn=","ffreq=", "cchnl=", "srmst=", "ermed=", "wdwtm=", "didtm=", "ohdtm=", "bdwtm=", "Drmdi=", "Bbimd=", "Mrmmd=", "Cprct=", "Iprit=", "Wprdt=", "Oprmd=", "Svst=", "Fprel="
            ])
        except getopt.GetoptError as err:
            print(err)
            print ('ldausbcli.py argument error')
            sys.exit(2)
        for opt, arg in opts:
            if opt == '-h':
                print ('ldausbcli.py -i <deviceid> -a <attenuation> -f <frequency> -c <channel> -s <rampstart> -e <rampend> -w <dwelltime> -d < idletime> -o <holdtime> -b <bidirectional-dwelltime> -D <rampdirection[0-Up,1-Down]> -M <rampmode[0-Once,1-Continuous,2-Off]> -B <rampbidirectional[0-Unidirectional,1-Bidirectional]> -C <profilecount> -I <profileidletime> -W <profiledwelltime>  -O <profilemode[0-Off,1-Once,2-Repeat]> -S <savesetting> -F <profilefile> -r <read>')
                sys.exit()
            elif opt in ("-i", "--idevid"):
                if(int(arg) > 0):
                    devid = int(arg)-1
                    print ('Device ID:', devid+1)
                else:
                    print ('Wrong device id..Please Check !')
                    exit()

            # Set Attenuation
            elif opt in ("-a", "--aattn"):
                attenuation = arg
                print ('Attenuation:', attenuation)
                if(devid !=''):
                    if attobj.check_deviceexists(devid):
                        attobj.set_attenuation(self.devices_list[devid], attenuation)
                    else:
                        print ("Device not exists.. Please Check!")
                else:
                    if(len(self.devices_list) > 0):
                        for index in self.devices_list:
                            attobj.set_attenuation(index, attenuation)
                    else:
                        print("No Devices Connected.. Please Check!")

            # Set Frequency
            elif opt in ("-f", "--ffreq"):
                frequency = arg
                print ('frequency:', frequency)
                if(devid !=''):
                    if attobj.check_deviceexists(devid):
                        attobj.set_frequency(self.devices_list[devid], frequency)
                    else:
                        print ("Device not exists.. Please Check!")
                else:
                    if(len(self.devices_list) > 0):
                        for index in self.devices_list:
                            attobj.set_frequency(index, frequency)
                    else:
                        print("No Devices Connected.. Please Check!")

            # Set Channel
            elif opt in ("-c", "--cchnl"):
                channel = arg
                print ('channel:', channel)
                if(devid !=''):
                    if attobj.check_deviceexists(devid):
                        attobj.set_channel(self.devices_list[devid], channel)
                    else:
                        print ("Device not exists.. Please Check!")
                else:
                    if(len(self.devices_list) > 0):
                        for index in self.devices_list:
                            attobj.set_channel(index, channel)
                    else:
                        print("No Devices Connected.. Please Check!")

            # Set Rampstart
            elif opt in ("-s", "--srmst"):
                rampstart = arg
                print ('rampstart:', rampstart)
                if(devid !=''):
                    if attobj.check_deviceexists(devid):
                        attobj.set_rampstartattenuation(self.devices_list[devid], rampstart)
                    else:
                        print ("Device not exists.. Please Check!")
                else:
                    if(len(self.devices_list) > 0):
                        for index in self.devices_list:
                            attobj.set_rampstartattenuation(index, rampstart)
                    else:
                        print("No Devices Connected.. Please Check!")

            # Set RampEnd
            elif opt in ("-e", "--ermed"):
                rampend = arg
                print ('rampend:', rampend)
                if(devid !=''):
                    if attobj.check_deviceexists(devid):
                        attobj.set_rampendattenuation(self.devices_list[devid], rampend)
                    else:
                        print ("Device not exists.. Please Check!")
                else:
                    if(len(self.devices_list) > 0):
                        for index in self.devices_list:
                            attobj.set_rampendattenuation(index, rampend)
                    else:
                        print("No Devices Connected.. Please Check!")

            # Set Dwell time
            elif opt in ("-w", "--wdwtm"):
                dwelltime = arg
                print ('dwelltime:', dwelltime)
                if(devid !=''):
                    if attobj.check_deviceexists(devid):
                        attobj.set_dwelltime(self.devices_list[devid], dwelltime)
                    else:
                        print ("Device not exists.. Please Check!")
                else:
                    if(len(self.devices_list) > 0):
                        for index in self.devices_list:
                            attobj.set_dwelltime(index, dwelltime)
                    else:
                        print("No Devices Connected.. Please Check!")

            # Set Idle time
            elif opt in ("-d", "--didtm"):
                idletime = arg
                print ('idletime:', idletime)
                if(devid !=''):
                    if attobj.check_deviceexists(devid):
                        attobj.set_idletime(self.devices_list[devid], idletime)
                    else:
                        print ("Device not exists.. Please Check!")
                else:
                    if(len(self.devices_list) > 0):
                        for index in self.devices_list:
                            attobj.set_idletime(index, idletime)
                    else:
                        print("No Devices Connected.. Please Check!")

            # Set hold time
            elif opt in ("-o", "--ohdtm"):
                holdtime = arg
                print ('holdtime:', holdtime)
                if(devid !=''):
                    if attobj.check_deviceexists(devid):
                        attobj.set_holdtime(self.devices_list[devid], holdtime)
                    else:
                        print ("Device not exists.. Please Check!")
                else:
                    if(len(self.devices_list) > 0):
                        for index in self.devices_list:
                            attobj.set_holdtime(index, holdtime)
                    else:
                        print("No Devices Connected.. Please Check!")

            # Set bi-directional dwell time
            elif opt in ("-b", "--bdwtm"):
                bddwelltime = arg
                print ('bddwelltime:', bddwelltime)
                if(devid !=''):
                    if attobj.check_deviceexists(devid):
                        attobj.set_bidirectionaldwelltime(self.devices_list[devid], bddwelltime)
                    else:
                        print ("Device not exists.. Please Check!")
                else:
                    if(len(self.devices_list) > 0):
                        for index in self.devices_list:
                            attobj.set_bidirectionaldwelltime(index, bddwelltime)
                    else:
                        print("No Devices Connected.. Please Check!")

            # Set ramp direction
            elif opt in ("-D", "--Drmdi"):
                rmpdir = arg
                print ('ramp direction:', rmpdir)
                if(devid !=''):
                    if attobj.check_deviceexists(devid):
                        attobj.set_rampdirection(self.devices_list[devid], rmpdir)
                    else:
                        print ("Device not exists.. Please Check!")
                else:
                    if(len(self.devices_list) > 0):
                        for index in self.devices_list:
                            attobj.set_rampdirection(index, rmpdir)
                    else:
                        print("No Devices Connected.. Please Check!")

            # Set ramp bi-direction
            elif opt in ("-B", "--Bbimd"):
                rmpdir = arg
                print ('ramp bi-direction:', rmpdir)
                if(devid !=''):
                    if attobj.check_deviceexists(devid):
                        attobj.set_rampbidirectional(self.devices_list[devid], rmpdir)
                    else:
                        print ("Device not exists.. Please Check!")
                else:
                    if(len(self.devices_list) > 0):
                        for index in self.devices_list:
                            attobj.set_rampbidirectional(index, rmpdir)
                    else:
                        print("No Devices Connected.. Please Check!")

            # Set ramp mode
            elif opt in ("-M", "--Mrmmd"):
                rmpmode = arg
                print ('ramp mode:', rmpmode)
                if(devid !=''):
                    if attobj.check_deviceexists(devid):
                        attobj.set_rampmode(self.devices_list[devid], rmpmode)
                    else:
                        print ("Device not exists.. Please Check!")
                else:
                    if(len(self.devices_list) > 0):
                        for index in self.devices_list:
                            attobj.set_rampmode(index, rmpmode)
                    else:
                        print("No Devices Connected.. Please Check!")

            # Set Profile count
            elif opt in ("-C", "--Cprct"):
                profilecount = arg
                print ('profile count:', profilecount)
                if(devid !=''):
                    if attobj.check_deviceexists(devid):
                        attobj.set_profilecount(self.devices_list[devid], profilecount)
                    else:
                        print ("Device not exists.. Please Check!")
                else:
                    if(len(self.devices_list) > 0):
                        for index in self.devices_list:
                            attobj.set_profilecount(index, profilecount)
                    else:
                        print("No Devices Connected.. Please Check!")

            # Set Profile idletime
            elif opt in ("-I", "--Iprit"):
                profileidletime = arg
                print ('profile idle-time:', profileidletime)
                if(devid !=''):
                    if attobj.check_deviceexists(devid):
                        attobj.set_profileidletime(self.devices_list[devid], profileidletime)
                    else:
                        print ("Device not exists.. Please Check!")
                else:
                    if(len(self.devices_list) > 0):
                        for index in self.devices_list:
                            attobj.set_profileidletime(index, profileidletime)
                    else:
                        print("No Devices Connected.. Please Check!")


            # Set Profile dwelltime
            elif opt in ("-W", "--Wprdt"):
                profiledwelltime = arg
                print ('profile dwell-time:', profiledwelltime)
                if(devid !=''):
                    if attobj.check_deviceexists(devid):
                        attobj.set_profiledwelltime(self.devices_list[devid], profiledwelltime)
                    else:
                        print ("Device not exists.. Please Check!")
                else:
                    if(len(self.devices_list) > 0):
                        for index in self.devices_list:
                            attobj.set_profiledwelltime(index, profiledwelltime)
                    else:
                        print("No Devices Connected.. Please Check!")

            # Set Profile mode
            elif opt in ("-O", "--Oprmd"):
                profilemode = arg
                print ('profile mode:', profilemode)
                if(devid !=''):
                    if attobj.check_deviceexists(devid):
                        attobj.set_profilemode(self.devices_list[devid], profilemode)
                    else:
                        print ("Device not exists.. Please Check!")
                else:
                    if(len(self.devices_list) > 0):
                        for index in self.devices_list:
                            attobj.set_profilemode(index, profilemode)
                    else:
                        print("No Devices Connected.. Please Check!")

            # Set savesettings
            elif opt in ("-S", "--Svst"):
                if(devid !=''):
                    if attobj.check_deviceexists(devid):
                        attobj.set_savesettings(self.devices_list[devid])
                    else:
                        print ("Device not exists.. Please Check!")
                else:
                    if(len(self.devices_list) > 0):
                        for index in self.devices_list:
                            attobj.set_savesettings(index)
                    else:
                        print("No Devices Connected.. Please Check!")

            # Set Profile File
            elif opt in ("-F", "--Fprel"):
                profilefilename = arg
                print ('profile Filename:', profilefilename)
                if(devid !=''):
                    if attobj.check_deviceexists(devid):
                        profilefile = open(profilefilename, 'r')
                        count = 0
                        dwelltime = 0
                        idletime = 0
                        profilelength = 0
                        profileindex=0
                        for linedata in profilefile.readlines():
                            if "dwell=" in linedata:
                                dwelltime = linedata.split("dwell=",1)[1]
                                attobj.set_profiledwelltime(self.devices_list[devid], int(float(dwelltime)*1000))
                            elif "idle=" in linedata:
                                idletime = linedata.split("idle=",1)[1]
                                attobj.set_profileidletime(self.devices_list[devid], int(float(idletime)*1000))
                            elif "length=" in linedata:
                                profilelength = linedata.split("length=",1)[1]
                                attobj.set_profilecount(self.devices_list[devid], int(profilelength))
                            else:
#                                print("Line{}: {}".format(profileindex, linedata.strip()))
#                                print('profiledata-1:%r',int(float(linedata.strip())))
                                attobj.set_profilelements(self.devices_list[devid])
                                profileindex = profileindex + 1
                                sleep(0.05) # delay 50mec

                        print("Reading File Done")
                    else:
                        print ("Device not exists.. Please Check!")
                else:
                    print("Enter the device id..Please Check!")

            elif opt in ("-r", "--rdata"):
                if(devid !=''):
                    if attobj.check_deviceexists(devid):
                        print("*****************Current Information of the device**********")
                        self.vnx.fnLDA_InitDevice(self.devices_list[devid])
                        print("Model Name:",attobj.get_modelname(self.devices_list[devid]))
                        print("Serial #:",attobj.get_serial_number(self.devices_list[devid]))
                        print("SW Version:",attobj.get_swversion())
                        print("Min Frequency(MHz):",attobj.get_minfrequency(self.devices_list[devid]))
                        print("Max Frequency(MHz):",attobj.get_maxfrequency(self.devices_list[devid]))
                        print("Min Attenuation(dB):",attobj.get_minattenuation(self.devices_list[devid]))
                        print("Max Attenuation(dB):",attobj.get_maxattenuation(self.devices_list[devid]))
                        print("Frequency(MHz):",attobj.get_currentfrequency(self.devices_list[devid]))
                        print("Attenuation(dB):",attobj.get_currentattenuation(self.devices_list[devid]))
                        print("Ramp Start Attn(dB):",attobj.get_rampstart(self.devices_list[devid]))
                        print("Ramp End Attn(dB):",attobj.get_rampend(self.devices_list[devid]))
                        print("Dwell Time(ms):",attobj.get_dwelltime(self.devices_list[devid]))
                        print("Idle Time(ms):",attobj.get_idletime(self.devices_list[devid]))
                        print("Hold Time:",attobj.get_holdtime(self.devices_list[devid]))
                        print("Bi-directional Dwell Time:",attobj.get_bidirectional_dwelltime(self.devices_list[devid]))
                        print("Profile Count:",attobj.get_profilecount(self.devices_list[devid]))
                        print("Profile Maxlength:",attobj.get_profilemaxlength(self.devices_list[devid]))
                        print("Profile dwelltime(ms):",attobj.get_profiledwelltime(self.devices_list[devid]))
                        print("Profile idletime(ms):",attobj.get_profileidletime(self.devices_list[devid]))
                        print("Profile index:",attobj.get_profileindex(self.devices_list[devid]))
                        self.vnx.fnLDA_CloseDevice(self.devices_list[devid])
                    else:
                        print ("Device not exists.. Please Check!")
                else:
                    if(len(self.devices_list) > 0):
                        for index in self.devices_list:
                            print("*****************Current Information of the devices**********")
                            self.vnx.fnLDA_InitDevice(index)
                            print("Device#:",index)
                            print("Model Name:",attobj.get_modelname(index))
                            print("Serial #:",attobj.get_serial_number(index))
                            print("SW Version:",attobj.get_swversion())
                            print("Min Frequency(MHz):",attobj.get_minfrequency(index))
                            print("Max Frequency(MHz):",attobj.get_maxfrequency(index))
                            print("Min Attenuation(dB):",attobj.get_minattenuation(index))
                            print("Max Attenuation(dB):",attobj.get_maxattenuation(index))
                            print("Frequency(MHz):",attobj.get_currentfrequency(index))
                            print("Attenuation(dB):",attobj.get_currentattenuation(index))
                            print("Ramp Start Attn(dB):",attobj.get_rampstart(index))
                            print("Ramp End Attn(dB):",attobj.get_rampend(index))
                            print("Dwell Time(ms):",attobj.get_dwelltime(index))
                            print("Idle Time(ms):",attobj.get_idletime(index))
                            print("Hold Time:",attobj.get_holdtime(index))
                            print("Bi-directional Dwell Time:",attobj.get_bidirectional_dwelltime(index))
                            print("Profile Count:",attobj.get_profilecount(index))
                            print("Profile Maxlength:",attobj.get_profilemaxlength(index))
                            print("Profile dwelltime(ms):",attobj.get_profiledwelltime(index))
                            print("Profile idletime(ms):",attobj.get_profileidletime(index))
                            print("Profile index:",attobj.get_profileindex(index))
                            self.vnx.fnLDA_CloseDevice(index)
                    else:
                        print("No Devices Connected.. Please Check!")

if __name__ == '__main__':
    pass
    # LabBrick_attn().main(sys.argv[1:])
