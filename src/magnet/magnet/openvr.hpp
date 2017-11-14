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

#include <magnet/GL/camera.hpp>
#include <openvr.h>

namespace magnet {
  class OpenVRTracker : public magnet::GL::Camera{
  public:
    
    OpenVRTracker(std::function<void(std::string)> log = [](std::string){}):
      Camera(0.1f, 30.0f),
      _vr(nullptr),
      _log(log),
      _eye(vr::Eye_Left)
    {
      _hmd_pose = magnet::GL::GLMatrix::identity();
    }

    void setEye(vr::EVREye x) {
      _eye = x;
    }
    
    virtual magnet::GL::GLMatrix getViewMatrix() const {
      switch(_eye){
      case vr::Eye_Left:
	return _eyePosLeft * _hmd_pose * magnet::GL::translate(0,1.5,0);
      case vr::Eye_Right:
	return _eyePosRight * _hmd_pose * magnet::GL::translate(0,1.5,0);
      default:
	M_throw() << "Bad enumeration in OpenVRTracker::getViewMatrix";
      }
    }
    
    virtual magnet::GL::GLMatrix getProjectionMatrix() const {
      switch(_eye){
      case vr::Eye_Left:
	return _projectionLeft;
      case vr::Eye_Right:
	return _projectionRight;
      default:
	M_throw() << "Bad enumeration in OpenVRTracker::getProjectionMatrix";
      }
    }
    
    //Not implemented!
    virtual void setUp(math::Vector newup, math::Vector axis = math::Vector{0,0,0}) {}

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

      //Initialise Compositor
      if (vr::VRCompositor())
	_log("Compositor initialised.");
      else {
	_log("Error: Compositor initialisation failed.");
	shutdown();
	return;
      }

      std::array<uint32_t, 2> dims = getRenderDims();
      resize(dims[0],dims[1]);
      
    }

    magnet::GL::FBO r_renderTarget;


    virtual magnet::GL::FBO& getResolveBuffer() {
      switch(_eye){
      case vr::Eye_Left:
	return _renderTarget;
      case vr::Eye_Right:
	return r_renderTarget;
      default:
	M_throw() << "Bad enumeration in OpenVRTracker::getProjectionMatrix";
      }
      return r_renderTarget;
    }
    
    virtual void deinit() {
      magnet::GL::Camera::deinit();
      r_renderTarget.deinit();
    }

    void submit() {
      vr::Texture_t tex = {reinterpret_cast<void*>(intptr_t(getResolveBuffer().getColorTexture()->getGLHandle())),
			   vr::TextureType_OpenGL, vr::ColorSpace_Gamma};
      vr::EVRCompositorError err = vr::VRCompositor()->Submit(_eye, &tex);
      if (err != vr::VRCompositorError_None)
	_log("Error: "+to_string(err));
    }

    
    virtual void resize(size_t width, size_t height) {
      magnet::GL::Camera::resize(width, height);

      std::shared_ptr<magnet::GL::Texture2D> r_colorTexture(new magnet::GL::Texture2D);
      r_colorTexture->init(width, height, GL_RGBA8);
      r_colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      r_colorTexture->parameter(GL_TEXTURE_MAX_LEVEL, 0);

      std::shared_ptr<magnet::GL::Texture2D> r_depthTexture(new magnet::GL::Texture2D);
      r_depthTexture->init(width, height, GL_DEPTH_COMPONENT);
      r_depthTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      r_depthTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      r_depthTexture->parameter(GL_TEXTURE_COMPARE_MODE, GL_NONE);
      
      r_renderTarget.init();
      r_renderTarget.attachTexture(r_colorTexture, 0);
      r_renderTarget.attachTexture(r_depthTexture);      
    }
    
    
    bool initialised() const { return _vr != nullptr; }
    
    std::array<uint32_t, 2> getRenderDims() const {
      uint32_t renderWidth, renderHeight;
      _vr->GetRecommendedRenderTargetSize(&renderWidth, &renderHeight);
      return std::array<uint32_t, 2>{renderWidth, renderHeight};
    }
    
    void shutdown() {
      if (_vr != nullptr) {
	deinit();
	vr::VR_Shutdown();
	_vr = nullptr;
	_log("Shutdown of VR complete.");
      }
    }
    
    void handleEvents() {
      	vr::VREvent_t event;
	while( _vr->PollNextEvent( &event, sizeof( event ) ) ) {
	  switch( event.eventType )
	    {
	    case vr::VREvent_TrackedDeviceActivated:
	      //SetupRenderModelForTrackedDevice( event.trackedDeviceIndex );
	      _log( "Device "+std::to_string(event.trackedDeviceIndex)+" attached." );
	      break;
	    case vr::VREvent_TrackedDeviceDeactivated:
	      _log( "Device "+std::to_string(event.trackedDeviceIndex)+" detached.");
	      break;
	    case vr::VREvent_TrackedDeviceUpdated:
	      _log( "Device "+std::to_string(event.trackedDeviceIndex)+" updated.");
	      break;
	    }
	}
	
	//// Process SteamVR controller state
	//for( vr::TrackedDeviceIndex_t unDevice = 0; unDevice < vr::k_unMaxTrackedDeviceCount; unDevice++ )
	//  {
	//    vr::VRControllerState_t state;
	//    if(_vr->GetControllerState( unDevice, &state, sizeof(state) ) )
	//      _tracked_devices[ unDevice ] = state.ulButtonPressed == 0;
	//  }
    }

    void getPosesAndSync() {
      vr::EVRCompositorError err = vr::VRCompositor()->WaitGetPoses(_tracked_devices.data(), vr::k_unMaxTrackedDeviceCount, nullptr, 0);
      
      if (err != vr::VRCompositorError_None)
	_log("Error: "+to_string(err));
      
      if (_tracked_devices[vr::k_unTrackedDeviceIndex_Hmd].bPoseIsValid) {
	_hmd_pose = magnet::math::inverse(convert(_tracked_devices[vr::k_unTrackedDeviceIndex_Hmd].mDeviceToAbsoluteTracking));
	_eyePosLeft = magnet::math::inverse(convert(_vr->GetEyeToHeadTransform(vr::Eye_Left)));
	_eyePosRight = magnet::math::inverse(convert(_vr->GetEyeToHeadTransform(vr::Eye_Right)));
	
	_projectionLeft = convert(_vr->GetProjectionMatrix(vr::Eye_Left, _zNearDist, _zFarDist));
	_projectionRight = convert(_vr->GetProjectionMatrix(vr::Eye_Right, _zNearDist, _zFarDist));
      }
      
    }
    
    void PostPresentHandoff() {
      vr::VRCompositor()->PostPresentHandoff();
    }
    
  protected:
    static magnet::GL::GLMatrix convert(vr::HmdMatrix44_t m) {
      magnet::GL::GLMatrix retval;
      for (size_t i(0); i < 4; ++i)
	for (size_t j(0); j < 4; ++j)
	  retval(i,j) = m.m[i][j];
      return retval;
    }

    static magnet::GL::GLMatrix convert(vr::HmdMatrix34_t m) {
      magnet::GL::GLMatrix retval;
      for (size_t i(0); i < 3; ++i)
	for (size_t j(0); j < 4; ++j)
	  retval(i,j) = m.m[i][j];

      retval(3,0) = 0;
      retval(3,1) = 0;
      retval(3,2) = 0;
      retval(3,3) = 1.0f;
      return retval;
    }
    
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

    static std::string to_string(vr::EVRCompositorError err) {
      switch (err) {
      case vr::VRCompositorError_None:
	return "None";
      case vr::VRCompositorError_RequestFailed:
	return "Request failed";
      case vr::VRCompositorError_IncompatibleVersion:
	return "Incompatible version";
      case vr::VRCompositorError_DoNotHaveFocus:
	return "Do not have focus";
      case vr::VRCompositorError_InvalidTexture:
	return "Invalid texture";
      case vr::VRCompositorError_IsNotSceneApplication:
	return "Is not a scene application";
      case vr::VRCompositorError_TextureIsOnWrongDevice:
	return "Texture is on wrong device";
      case vr::VRCompositorError_TextureUsesUnsupportedFormat:
	return "Texture uses unsupported format";
      case vr::VRCompositorError_SharedTexturesNotSupported:
	return "Shared textures are not supported";
      case vr::VRCompositorError_IndexOutOfRange:
	return "Index out of range";
      case vr::VRCompositorError_AlreadySubmitted:
	return "Texture already submitted";
      case vr::VRCompositorError_InvalidBounds:
	return "Invalid bounds";
      default:
	return "Unhandled VR Compositor Error "+std::to_string(err);
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

    magnet::GL::GLMatrix _projectionLeft, _projectionRight, _eyePosLeft, _eyePosRight, _hmd_pose;
    vr::EVREye _eye;
    std::array<vr::TrackedDevicePose_t, vr::k_unMaxTrackedDeviceCount> _tracked_devices;
  };
}
