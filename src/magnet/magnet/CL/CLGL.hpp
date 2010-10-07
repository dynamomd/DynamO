/*    DYNAMO:- Event driven molecular dynamics simulator 
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
	initDevice();
	_commandQ =cl::CommandQueue(_context, _device);
	_initialised = true;
      }
 
      inline cl::Platform getPlatform() 
      { 
	if (!_initialised) M_throw() << "Not initialised()!";
	return _platform; 
      }

      inline cl::Context getContext() 
      { 
	if (!_initialised) M_throw() << "Not initialised()!";
	return _context; 
      }

      inline cl::Device getDevice() 
      { 
	if (!_initialised) M_throw() << "Not initialised()!";
	return _device; 
      }

      inline cl::CommandQueue getCommandQueue() 
      { 
	if (!_initialised) M_throw() << "Not initialised()!";
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
      
        //Now cycle through the platforms trying to get a context
        for (std::vector<cl::Platform>::const_iterator iPtr = platforms.begin();
             iPtr != platforms.end(); ++iPtr)
	  {
	    try {
	      _context = getCLGLContext(*iPtr);
	      //Success! now set the platform and return!
	      _platform = *iPtr;
	     cl::GLBuffer::hostTransfers() = false;
	      return;
	    } catch (...)
	      {/*Failed so we just carry on*/}
	  }
	
	//No CLGL platform was found so just give the first platform
	_platform = platforms.front();
	//Make sure to set host transfers on!
       cl::GLBuffer::hostTransfers() = true;
	
	cl_context_properties cpsFallBack[] = {CL_CONTEXT_PLATFORM, 
					       (cl_context_properties)getPlatform()(),
					       0};
	try {
	  _context =cl::Context(CL_DEVICE_TYPE_ALL, cpsFallBack, NULL, NULL);
	} catch (...)
	  {
	    throw std::runtime_error("Failed to create a normal OpenCL context!");
	  }
      }
      
      inline cl::Context getCLGLContext(::cl::Platform clplatform)
      {
	if (glXGetCurrentContext() == NULL)
	  throw std::runtime_error("Failed to obtain the GL context");
	
	cl_context_properties cpsGL[] = { CL_CONTEXT_PLATFORM, (cl_context_properties) clplatform(),
					  CL_GLX_DISPLAY_KHR, (cl_context_properties) glXGetCurrentDisplay(),
					  CL_GL_CONTEXT_KHR, (cl_context_properties) glXGetCurrentContext(), 0};
	
	::cl::Context clcontext;
	
	//create the OpenCL context 
	try {
	  clcontext =cl::Context(CL_DEVICE_TYPE_ALL, cpsGL, NULL, NULL);
	} catch(...)
	  {
	    throw std::runtime_error("Could not generate CLGL context");
	  }
	
	return clcontext;
      }

      inline void initDevice()
      {
	//Grab the first device
	std::vector<cl::Device> devices = getContext().getInfo<CL_CONTEXT_DEVICES>();
	
	//Default to the first device
	_device = devices.front();
	
	//But check if there is a GPU to use
	for (std::vector<cl::Device>::const_iterator iPtr = devices.begin();
	     iPtr != devices.end(); ++iPtr)
	  if (iPtr->getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_GPU)
	    {
	      //Take the first GPU
	      _device = *iPtr;
	      break;
	    }
      }
    };
  }
}
