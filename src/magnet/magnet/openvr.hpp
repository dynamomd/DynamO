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
      _HMD(nullptr),
      _log(log)
    {}

    ~OpenVRTracker() { shutdown(); }

    void setLog(std::function<void(std::string)> log = [](std::string){}) {
      _log = log;
    }
    
    void init() {
      if (_HMD)
	_log("OpenVR already initialized!\n");
      
      if (!vr::VR_IsRuntimeInstalled())
	{
	  _log("Error: No OpenVR runtime library detected, have you installed SteamVR?\n");
	  return;
	}
      else
	_log(std::string("Using OpenVR runtime at ")+vr::VR_RuntimePath()+std::string("\n"));
	
      if (!vr::VR_IsHmdPresent())
	{
	  _log("Error: No HMD detected, have you started SteamVR?\n");
	  return;
	}
	    
      vr::EVRInitError eError = vr::VRInitError_None;
      _HMD = vr::VR_Init(&eError, vr::VRApplication_Scene);
	
      if (eError != vr::VRInitError_None)
	{
	  _HMD = nullptr;
	  std::string error = vr::VR_GetVRInitErrorAsEnglishDescription(eError);
	  _log("Error: "+error+"\n");
	  return;
	}
    }

    bool initialised() const { return _HMD != nullptr; }

    void shutdown() {
      if (_HMD != nullptr) {
	vr::VR_Shutdown();
	_HMD = nullptr;
      }
    }
    
  protected:
    vr::IVRSystem* _HMD;
    std::function<void(std::string)> _log;
  };
}
