/*    dynamo:- Event driven molecular dynamics simulator 
 *    http://www.marcusbannerman.co.uk/dynamo
 *    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>
 *
 *    This program is free software: you can redistribute it and/or
 *    modify it under the terms of the GNU General Public License
 *    version 3 as published by the Free Software Foundation.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include <magnet/CL/GLBuffer.hpp>
#include <magnet/exception.hpp>

namespace magnet {
  namespace CL {
    class CLGLState {
    public:
      CLGLState():_initialised(false) {}

      //Only call this when you're in a valid OpenGL context
      inline void init()
      {
	if (_initialised) M_throw() << "Initialising twice!";
	initContext();
	_commandQ =cl::CommandQueue(_context, _device);
	_initialised = true;
      }
 
      inline cl::Platform& getPlatform() 
      { 
#ifdef MAGNET_DEBUG
	if (!_initialised) M_throw() << "Not initialised()!";
#endif
	return _platform; 
      }

      inline cl::Context& getContext() 
      { 
#ifdef MAGNET_DEBUG
	if (!_initialised) M_throw() << "Not initialised()!";
#endif
	return _context; 
      }

      inline cl::Device& getDevice() 
      { 
#ifdef MAGNET_DEBUG
	if (!_initialised) M_throw() << "Not initialised()!";
#endif
	return _device; 
      }

      inline cl::CommandQueue& getCommandQueue() 
      { 
#ifdef MAGNET_DEBUG
	if (!_initialised) M_throw() << "Not initialised()!";
#endif
	return _commandQ; 
      }
      
    private:

     cl::Platform _platform;
     cl::Context _context;
     cl::Device _device;
     cl::CommandQueue _commandQ;

      bool _initialised;

      inline void initContext()
      {
	std::vector<cl::Platform> platforms;
       cl::Platform::get(&platforms);
      
        //Now cycle through the platforms trying to get a context with a GPU
        for (std::vector<cl::Platform>::const_iterator iPtr = platforms.begin();
             iPtr != platforms.end(); ++iPtr)
	  {
	    std::vector<cl::Device> devices;
	    try {
	      iPtr->getDevices(CL_DEVICE_TYPE_GPU, &devices);
	    } catch (...)
	      { continue; }

	    for (std::vector<cl::Device>::const_iterator devPtr = devices.begin();
		 devPtr != devices.end(); ++devPtr)
	      {
		if (!getCLGLContext(*iPtr, *devPtr)) continue;
		
		//Success! now set the platform+device and return!
		_platform = *iPtr;
		cl::GLBuffer::hostTransfers() = false;
		_device = *devPtr;
		return;
	      }
	  }

        //Try and see if there's another OpenCL-OpenGL platform available
        for (std::vector<cl::Platform>::const_iterator iPtr = platforms.begin();
             iPtr != platforms.end(); ++iPtr)
	  {
	    std::vector<cl::Device> devices;
	    try {
	      iPtr->getDevices(CL_DEVICE_TYPE_ALL, &devices);
	    } catch (...)
	      { continue; }

	    for (std::vector<cl::Device>::const_iterator devPtr = devices.begin();
		 devPtr != devices.end(); ++devPtr)
	      {
		if (!getCLGLContext(*iPtr, *devPtr)) continue;
		
		//Success! now set the platform+device and return!
		_platform = *iPtr;
		cl::GLBuffer::hostTransfers() = false;
		_device = *devPtr;
		return;
	      }
	  }
	
	//No CLGL platform was found so just give the first platform
	_platform = platforms.front();
	//Make sure to set host transfers on!
       cl::GLBuffer::hostTransfers() = true;
	
	cl_context_properties cpsFallBack[] = {CL_CONTEXT_PLATFORM, 
					       (cl_context_properties)_platform(),
					       0};
	try {
	  _context = cl::Context(CL_DEVICE_TYPE_ALL, cpsFallBack, NULL, NULL);
	  _device = _context.getInfo<CL_CONTEXT_DEVICES>().front();
	} catch (...)
	  {
	    throw std::runtime_error("Failed to create a normal OpenCL context!");
	  }
      }
      
      inline bool getCLGLContext(::cl::Platform clplatform, ::cl::Device dev)
      {
	if (glXGetCurrentContext() == NULL)
	  throw std::runtime_error("Failed to obtain the GL context");
	
	cl_context_properties cpsGL[] = { CL_CONTEXT_PLATFORM, (cl_context_properties) clplatform(),
					  CL_GLX_DISPLAY_KHR, (cl_context_properties) glXGetCurrentDisplay(),
					  CL_GL_CONTEXT_KHR, (cl_context_properties) glXGetCurrentContext(), 0};
	
	std::vector< ::cl::Device> devlist;
	devlist.push_back(dev);
	
	//create the OpenCL context 
	try {
	  _context =cl::Context(devlist, cpsGL);
	} catch(...)
	  //If we fail, it's probably because it's not the correct
	  //device for the GL context
	  { return false; }
	
	return true;
      }
    };
  }
}
