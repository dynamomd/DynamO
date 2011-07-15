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

//Here we have the correct order of GL includes
#define GL_GLEXT_PROTOTYPES
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>
#include <GL/glext.h>
#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>

#include <magnet/GL/detail/traits.hpp>
#include <magnet/GL/detail/enums.hpp>
#include <magnet/exception.hpp>
#include <map>

namespace magnet {
  namespace GL {
    template<class T> class Buffer;

    /*! \brief Class representing an OpenGL context (and its
        associated OpenCL context if required).
     *
     * The purpose of this class is to track the state of an OpenGL
     * context, allowing queries as to the currently bound shader,
     * textures and so on.
     *
     * This class also establishes the corresponding CL context for
     * the GL context.
     *
     */
    class Context
    {
    public:
      /*! \brief Method to fetch the current OpenGL context.
       *
       * This function is used to make sure that whenever the context
       * is requested, the same copy is always returned.
       */
      inline static Context& getContext()
      {
	static std::map<ContextKey, Context> contexts;

	ContextKey key = getCurrentContextKey();

	std::map<ContextKey, Context>::iterator ctxPtr = contexts.find(key);
	if (ctxPtr == contexts.end())
	  contexts[key]._context = getCurrentContextKey();
	  
	return contexts[key];
      }

      /** @name The OpenGL state functions. */
      /**@{*/
      /*! \brief Draw the elements described in the passed buffer
       *
       */
      template<class T>
      inline void drawElements(Buffer<T>& elementBuffer, element_type::Enum etype)
      {
	elementBuffer.bind(buffer_targets::ELEMENT_ARRAY);
	glDrawElements(etype, elementBuffer.size(), detail::c_type_to_gl_enum<T>::val, 0);
      }
      /**@}*/

      /** @name The OpenCL-OpenGL interface. */
      /**@{*/
      //! \brief Fetch the OpenCL platform for this OpenGL context.
      inline const cl::Platform& getCLPlatform() 
      { initCL(); return _clplatform; }

      //! \brief Fetch the OpenCL context for this OpenGL context.
      inline const cl::Context& getCLContext() 
      { initCL(); return _clcontext; }

      //! \brief Fetch the OpenCL device for this OpenGL context.
      inline const cl::Device& getCLDevice() 
      { initCL(); return _cldevice; }

      //! \brief Fetch the OpenCL command queue for this OpenGL context.
      inline const cl::CommandQueue& getCLCommandQueue() 
      { initCL(); return _clcommandQ; }

    protected:
      //! \brief The OpenCL platform for this GL context.
      cl::Platform _clplatform;
      //! \brief The OpenCL context for this GL context.
      cl::Context _clcontext;
      //! \brief The OpenCL device for this GL context.
      cl::Device _cldevice;
      //! \brief The OpenCL command queue for this GL context.
      cl::CommandQueue _clcommandQ;

      //! \brief Flag set if the OpenCL state has been initialised.
      bool _clInitialised;

      /*! \brief If a matching OpenCL context does not exist, it will
       * create one from the current OpenGL context.
       */
      inline void initCL()
      {
	if (_clInitialised) return;
	_clInitialised = true;

	initOpenCLContext();
	_clcommandQ = cl::CommandQueue(_clcontext, _cldevice);
      }

      /*! \brief Initializes an OpenCL context, platform and device
       * from the current OpenGL context.
       */
      inline void initOpenCLContext()
      {
	if (initOpenCLContextType(CL_DEVICE_TYPE_GPU)) return;
	if (initOpenCLContextType(CL_DEVICE_TYPE_ALL)) return;

	M_throw() << "Failed to create an OpenCL context from the OpenGL context!";
      }

      /*! \brief Attempts to Initializes an OpenCL context given a
       * specific device type.
       *
       * \returns True if a matching device and context could be
       * constructed.
       */
      inline bool initOpenCLContextType(cl_device_type devType)
      {
	std::vector<cl::Platform> platforms;
        cl::Platform::get(&platforms);
	
	//Now cycle through the platforms trying to get a context with a GPU
	for (std::vector<cl::Platform>::const_iterator iPtr = platforms.begin();
	     iPtr != platforms.end(); ++iPtr)
	  {
	    std::vector<cl::Device> devices;
	    try {
	      iPtr->getDevices(devType, &devices);
	    } catch (...)
	      { continue; }
    
	    for (std::vector<cl::Device>::const_iterator devPtr = devices.begin();
		 devPtr != devices.end(); ++devPtr)
	      {
		if (!getCLGLContext(*iPtr, *devPtr)) continue;
		
		//Success! now set the platform+device and return!
		_clplatform = *iPtr;
		_cldevice = *devPtr;
		return true;
	      }
	  }

	return false;
      }

      /*! \fn bool getCLGLContext(cl::Platform clplatform, cl::Device dev)
       * 
       * \brief This is a system specific command to build a OpenCL
       * context from the current OpenGL context.
       *
       * \returns True if a context could be created for the passed
       * device and platform.
       */

      /*! \fn ContextKey getCurrentContextKey()
       * 
       * \brief This is a system-specific command to fetch the system-specific
       * OpenGL-context-handle of the current GL context.
       */

      ////////////////X11 specific bindings////////////////
      inline bool getCLGLContext(cl::Platform clplatform, cl::Device dev)
      {
	cl_context_properties cpsGL[] = { CL_CONTEXT_PLATFORM, (cl_context_properties) clplatform(),
					  CL_GLX_DISPLAY_KHR, (cl_context_properties) glXGetCurrentDisplay(),
					  CL_GL_CONTEXT_KHR, (cl_context_properties) _context, 0};
	
	std::vector<cl::Device> devlist;
	devlist.push_back(dev);
	
	//create the OpenCL context 
	try {
	  _clcontext = cl::Context(devlist, cpsGL);
	} catch(...)
	  //If we fail, it's probably because it's not the correct
	  //device for the GL context
	  { return false; }
	
	return true;
      }

      /**@}*/

      typedef GLXContext ContextKey;

      inline static ContextKey getCurrentContextKey()
      { 
	ContextKey key = glXGetCurrentContext();
	if (!key) M_throw() << "Not in a valid GLX context";
	return key; 
      }
      /////////////////////////////////////////////////////

      //! \brief Only allow \ref getContext() to instantiate Context objects.
      Context(): _clInitialised(false) {}

      //! \brief The system-dependent handle to the GL context.
      ContextKey _context;

      //! \brief Friend statement needed for \ref getContext().
      friend class std::map<ContextKey, Context>;
    };
  }
}
