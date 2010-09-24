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

namespace magnet {
  namespace cl {
    class GLInterop {
    private:
      static inline ::cl::Platform& setPlatform() 
      {
	static ::cl::Platform _platform;
	return _platform;
      }

      static inline ::cl::Context& setContext() 
      {
	static ::cl::Context _context;
	return _context;
      }

    public:

      static inline const ::cl::Platform& getPlatform() 
      {
	return setPlatform();
      }

      static inline const ::cl::Context& getContext() 
      {
	return setContext();
      }

      //Setup the best OpenCL/GL context we can
      static bool init()
      {
	std::vector< ::cl::Platform> platforms;
        ::cl::Platform::get(&platforms);
      
        //Now cycle through the platforms trying to get a context
        for (std::vector< ::cl::Platform>::const_iterator iPtr = platforms.begin();
             iPtr != platforms.end(); ++iPtr)
	  {
	    try {
	      setContext() = getCLGLContext(*iPtr);
	      //Success! now set the platform and return!
	      setPlatform() = *iPtr;
	      ::cl::GLBuffer::hostTransfers() = false;
	      return true;
	    } catch (...)
	      {/*Failed so we just carry on*/}
	  }

	//No CLGL platform was found so just give the first platform
	setPlatform() = platforms.front();
	//Make sure to set host transfers on!
        ::cl::GLBuffer::hostTransfers() = true;
	
	cl_context_properties cpsFallBack[] = {CL_CONTEXT_PLATFORM, 
					       (cl_context_properties)getPlatform()(),
					       0};
	try {
	  setContext() = ::cl::Context(CL_DEVICE_TYPE_ALL, cpsFallBack, NULL, NULL);
	} catch (::cl::Error& err)
	  {
	    throw std::runtime_error("Failed to create a normal OpenCL context!");
	  }

	return true;
      }

      static ::cl::Context getCLGLContext(::cl::Platform clplatform)
      {
	GLXContext GLContext = glXGetCurrentContext();
	
	if (GLContext == NULL)
	  throw std::runtime_error("Failed to obtain the GL context");
	
	cl_context_properties cpsGL[] = { CL_CONTEXT_PLATFORM, 
					  (cl_context_properties)clplatform(),
					  CL_GLX_DISPLAY_KHR, (intptr_t) glXGetCurrentDisplay(),
					  CL_GL_CONTEXT_KHR, (intptr_t) GLContext, 0};
	
	::cl::Context clcontext;
	
	//create the OpenCL context 
	try {
	  clcontext = ::cl::Context(CL_DEVICE_TYPE_ALL, cpsGL, NULL, NULL);
	} catch(::cl::Error& err)
	  {
	    throw std::runtime_error("Could not generate CLGL context");
	  }
	
	return clcontext;
      }
    };
  }
}
