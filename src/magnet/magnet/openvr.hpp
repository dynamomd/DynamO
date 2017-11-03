/*
   (C) copyright 2017, Marcus Bannerman, All Rights Reserved

   LICENSING INFORMATION:
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <openvr.h>

namespace magnet {
  class OpenVRTracker {
  public:
    OpenVRTracker(std::function<void(std::string)> log = [](std::string){}):
      _vr(nullptr),
      _log(log)
    {}

    ~OpenVRTracker() { shutdown(); }

    void setLog(std::function<void(std::string)> log = [](std::string){}) {
      _log = log;
    }
    
    void init() {
      if (_vr)
	_log("OpenVR already initialized!");
      
      if (!vr::VR_IsRuntimeInstalled())
	{
	  _log("Error: No OpenVR runtime library detected, have you installed SteamVR?");
	  return;
	}
      else
	_log(std::string("Using OpenVR runtime at ")+vr::VR_RuntimePath());
	
      if (!vr::VR_IsHmdPresent())
	{
	  _log("Error: No HMD detected, have you started SteamVR?");
	  return;
	}
	    
      vr::EVRInitError eError = vr::VRInitError_None;
      _vr = vr::VR_Init(&eError, vr::VRApplication_Scene);
	
      if (eError != vr::VRInitError_None)
	{
	  _vr = nullptr;
	  _log("Error: "+to_string(eError));
	  return;
	}

      std::string val;
      bool success;

      for (vr::TrackedDeviceIndex_t i(0); i < vr::k_unMaxTrackedDeviceCount; ++i)
	if (_vr->IsTrackedDeviceConnected(i)) {
	  std::tie(success, val) = GetTrackedDeviceString(i, vr::Prop_TrackingSystemName_String);
	  
	  if (!success) {
	    _log("Error: While fetching Device#"+std::to_string(i)+" name, "+val);
	    shutdown();
	    return;
	  } else
	    _log("Device#"+std::to_string(i)+" name: "+val);
	  
	  std::tie(success, val) = GetTrackedDeviceString(i, vr::Prop_SerialNumber_String);
	  if (!success) {
	    _log("Error: While fetching Device#"+std::to_string(i)+" serial, "+val);
	    shutdown();
	    return;
	  } else
	    _log("Device#"+std::to_string(i)+" serial: "+val);
	  
	  _log("Device#"+std::to_string(i)+" class: "+to_string(_vr->GetTrackedDeviceClass(i)));
	}
    }

    bool initialised() const { return _vr != nullptr; }

    void shutdown() {
      if (_vr != nullptr) {
	vr::VR_Shutdown();
	_vr = nullptr;
      }
    }

  protected:
    static std::string to_string(vr::EVRInitError e) {
      return vr::VR_GetVRInitErrorAsEnglishDescription(e);
    }

    static std::string to_string(vr::TrackedPropertyError e) {
      switch(e) {
      case vr::TrackedProp_Success:
	return "Success!";
      case vr::TrackedProp_WrongDataType:
	return "The property was requested with the wrong typed function.";
      case vr::TrackedProp_NotYetAvailable:
	return "The property is not yet available.";
      case vr::TrackedProp_PermissionDenied:
	return "Permission denied.";
      case vr::TrackedProp_InvalidOperation:
	return "Invalid operation";
      case vr::TrackedProp_WrongDeviceClass:
	return "The property was requested on a tracked device with the wrong class.";
      case vr::TrackedProp_BufferTooSmall:
	return "The string property will not fit in the provided buffer. The buffer size needed is returned.";
      case vr::TrackedProp_UnknownProperty:
	return "The property enum value is unknown.";
      case vr::TrackedProp_InvalidDevice:
	return "The tracked device index was invalid.";
      case vr::TrackedProp_CouldNotContactServer:
	return "OpenVR could not contact vrserver to query the device for this property.";
      case vr::TrackedProp_ValueNotProvidedByDevice:
	return "The driver for this device returned that it does not provide this specific property for this device.";
      case vr::TrackedProp_StringExceedsMaximumLength:
	return "The string property value returned by a driver exceeded the maximum property length of 32k.";
      default:
	return "Unhandled tracked property error "+std::to_string(e);
      }
    }

    static std::string to_string(vr::TrackedDeviceClass c) {
      switch(c) {
      case vr::TrackedDeviceClass_Invalid:
	return "Invalid - no device";
      case vr::TrackedDeviceClass_HMD:
	return "HMD device";
      case vr::TrackedDeviceClass_Controller:
	return "Controller device";
      case vr::TrackedDeviceClass_GenericTracker:
	return "Tracking generic device";
      case vr::TrackedDeviceClass_TrackingReference:
	return "Tracking reference device";
      case vr::TrackedDeviceClass_DisplayRedirect:
	return "Display redirect device";
      default:
	return "Unhandled device class "+std::to_string(c);
      }
    }

    std::tuple<bool, std::string> GetTrackedDeviceString(vr::TrackedDeviceIndex_t unDevice, vr::TrackedDeviceProperty prop)
    {
      vr::TrackedPropertyError error;
      uint32_t unRequiredBufferLen = _vr->GetStringTrackedDeviceProperty(unDevice, prop, NULL, 0, &error);

      if ((error != vr::TrackedProp_Success) && (error != vr::TrackedProp_BufferTooSmall)) {
	return std::make_tuple(false, to_string(error));
      }
      
      if(unRequiredBufferLen == 0)
	return std::make_tuple(true, "");
      
      char *pchBuffer = new char[ unRequiredBufferLen ];
      unRequiredBufferLen = _vr->GetStringTrackedDeviceProperty( unDevice, prop, pchBuffer, unRequiredBufferLen, &error);

      if (error != vr::TrackedProp_Success) {
	return std::make_tuple(false, to_string(error));
      }
      
      std::string sResult = pchBuffer;
      delete [] pchBuffer;
      return std::make_tuple(true, sResult);
    }
    
    vr::IVRSystem* _vr;
    std::function<void(std::string)> _log;
  };
}
