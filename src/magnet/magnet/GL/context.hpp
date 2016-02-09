/*    dynamo:- Event driven molecular dynamics simulator 
 *    http://www.dynamomd.org
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
#include <GL/glew.h>
#include <GL/glx.h>
#include <GL/freeglut.h>

#ifdef MAGNET_CLGL
# define __CL_ENABLE_EXCEPTIONS
# include <CL/cl.hpp>
#endif

#include <magnet/GL/detail/traits.hpp>
#include <magnet/GL/detail/enums.hpp>
#include <magnet/GL/detail/typesafe_get.hpp>
#include <magnet/GL/detail/error_check.hpp>
#include <magnet/thread/taskQueue.hpp>
#include <magnet/exception.hpp>
#include <magnet/GL/matrix.hpp>
#include <memory>
#include <set>
#include <map>
#include <iostream>

namespace magnet {
  namespace GL {
    namespace shader { namespace detail { class Shader; } }
    template<class T> class Buffer;

    /*! \brief Class representing an OpenGL context (and its
      associated OpenCL context if required).
	
      The purpose of this class is to track the state of an OpenGL
      context, allowing queries as to the currently bound shader,
      textures and so on. It also tracks the GL state to minimise the
      number of GL state changes and redundant state changes are
      ignored.
     
      This class also establishes the corresponding OpenCL context for
      the GL context on access.
     */
    class Context
    {
    protected:
      /** @name Platform specific code. */
      /**@{*/
      typedef GLXContext ContextKey;

      inline static ContextKey getCurrentContextKey()
      { 
	ContextKey key = glXGetCurrentContext();
	detail::errorCheck();
	if (!key) M_throw() << "Not in a valid GLX context";
	return key; 
      }
      /**@}*/

    public:
      /*!\brief The reference counting type to use for a holding a
         reference to a context.
       */
      typedef std::shared_ptr<Context> ContextPtr;

      /*! \brief Method to fetch the current OpenGL context.
       
        This function is used to make sure that whenever the context
        is requested, the same copy is always returned.
       */
      inline static ContextPtr getContext()
      {
	typedef std::map<ContextKey, ContextPtr> ContextStoreType; 
	static ContextStoreType contexts;

	ContextKey key = getCurrentContextKey();

	ContextStoreType::iterator ctxPtr = contexts.find(key);
	if (ctxPtr == contexts.end())
	  {
	    contexts[key].reset(new Context);
	    contexts[key]->init();
	  }
	  
	return contexts[key];
      }
      
      /** @name Vertex attribute array interface. */
      /**@{*/
      /*! \brief Enables a vertex attribute array index. */
      inline void enableAttributeArray(GLuint attrnum)
      {
	if (attrnum >= _vertexAttributeState.size())
	  M_throw() << "Attribute index out of range";
	if (_vertexAttributeState[attrnum].active) return;
	glEnableVertexAttribArray(attrnum);
	detail::errorCheck();
	_vertexAttributeState[attrnum].active = true;
      }

      /*! \brief Disables a vertex attribute array index. */
      inline void disableAttributeArray(GLuint attrnum)
      {
	if (attrnum >= _vertexAttributeState.size())
	  M_throw() << "Attribute index out of range";
	if (!(_vertexAttributeState[attrnum].active)) return;
	glDisableVertexAttribArray(attrnum);
	detail::errorCheck();
	_vertexAttributeState[attrnum].active = false;
      }

      /*! \brief Disable all active vertex attribute arrays. */
      inline void cleanupAttributeArrays()
      {
	resetInstanceTransform();
	for (GLuint i(0); i < _vertexAttributeState.size(); ++i)
	  {
	    disableAttributeArray(i);
	    detail::errorCheck();
	  }
      }

      /*! \brief Sets the value of a vertex attribute, if no attribute
        array is bound.
       
        This function only sets the state if it has been updated.
       */
      inline void setAttribute(GLuint idx, GLfloat x = 0, GLfloat y = 0, GLfloat z = 0, GLfloat w = 0)
      {
	if (idx >= _vertexAttributeState.size())
	  M_throw() << "Attribute index out of range";

	if (idx == 0)
	  M_throw() << "Cannot set the value of the 0th vertex attribute.";

	std::array<GLfloat, 4> newval = {{x,y,z,w}};

#ifdef MAGNET_DEBUG
	std::array<GLfloat, 4> oldval;
	glGetVertexAttribfv(idx, GL_CURRENT_VERTEX_ATTRIB, &oldval[0]);
	detail::errorCheck();

	if (oldval != _vertexAttributeState[idx].current_value)
	  M_throw() << "Vertex attribute state changed without using the GL context!";
#endif

	if (newval == _vertexAttributeState[idx].current_value) return;
	_vertexAttributeState[idx].current_value = newval;

	glVertexAttrib4f(idx, newval[0], newval[1], newval[2], newval[3]);
	detail::errorCheck();
      }
      
      /*! \brief Sets the divisor of a vertex attribute.
       
        The divisor is used in instancing to set the rate at which
        vertex attributes are incremented.
       */
      inline void setAttributeDivisor(GLuint idx, GLuint divisor)
      {
	if (idx >= _vertexAttributeState.size())
	  M_throw() << "Attribute index out of range";

	if (divisor == _vertexAttributeState[idx].divisor) return;
	_vertexAttributeState[idx].divisor = divisor;
	glVertexAttribDivisorARB(idx, divisor);
	detail::errorCheck();
      }

      /*! \brief The index of the automatically-indexed position
        vertex attribute.
       
        This index for the vertex position is set in the OpenGL
        standard.

        \sa detail::Shader::build()
       */
      static const GLuint vertexPositionAttrIndex = 0;

      /*! \brief The index of the automatically-indexed position color
        attribute.

        \sa detail::Shader::build()
       */
      static const GLuint vertexColorAttrIndex = 1;

      /*! \brief The index of the automatically-indexed normal vertex
        attribute.
	
	\sa detail::Shader::build()
       */
      static const GLuint vertexNormalAttrIndex = 2;

      /*! \brief The index of the automatically-indexed instance
        origin vertex attribute.

        \sa detail::Shader::build()
       */
      static const GLuint instanceOriginAttrIndex = 3;

      /*! \brief The index of the automatically-indexed instance
        orientation vertex attribute.
       
	\sa detail::Shader::build()
       */
      static const GLuint instanceOrientationAttrIndex = 4;

      /*! \brief The index of the automatically-indexed instance scale
        vertex attribute.
       
	\sa detail::Shader::build()
       */
      static const GLuint instanceScaleAttrIndex = 5;

      /*! \brief The index of the automatically-indexed 
       * texture coordinate vertex attribute.
       * \sa detail::Shader::build()
       */
      static const GLuint vertexTexCoordAttrIndex = 6;

      /*! \brief Convenience function to set the vertex attribute
        representing the color in a shader.
       
        This uses the \ref vertexColorAttrIndex value for the index of
        the color attribute.
       */
      inline void color(GLfloat r = 0, GLfloat g = 0, GLfloat b = 0, GLfloat a = 1) 
      { setAttribute(vertexColorAttrIndex, r, g, b, a); }

      /*! \brief Convenience function to set the instance rotation.
       
        This uses the \ref vertexColorAttrIndex value for the index of
        the color attribute.
       
        \param angle The rotation angle in radians.
        \param axis The rotation axis.
       */
      inline void rotation(GLfloat angle, math::Vector axis) 
      {
	GLfloat s = std::sin(angle / 2);
	GLfloat c = std::cos(angle / 2);
	setAttribute(instanceOrientationAttrIndex, axis[0] * s, axis[1] * s, axis[2] * s, c); 
      }

      /*! \brief Resets the vertex attributes used in instancing to
        avoid unintended transformations of the instanced object.
       */
      inline void resetInstanceTransform()
      {
	setAttribute(instanceOriginAttrIndex, 0, 0, 0, 0);
	setAttribute(instanceOrientationAttrIndex, 0, 0, 0, 1);
	setAttribute(instanceScaleAttrIndex, 1, 1, 1, 1);
      }
      /**@}*/

#ifdef MAGNET_CLGL
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
      /**@}*/
#endif

      /*! \brief Returns the currently attached shader program.
       
        The currently attached shader program is managed by the \ref
        Shader classes \ref Shader::attach() and \ref Shader::detach()
        functions.
       */
      inline shader::detail::Shader& getAttachedShader()
      {
	if (_shaderStack.empty())
	  M_throw() << "No shader attached to the GL context!";

	return *_shaderStack.back();
      }
      
      /*! \brief Sets the current viewport.
       
        \param x The coordinate of the leftmost pixel in the viewport.
        \param y The coordinate of the lowest pixel in the viewport. 
        \param width The width of the viewport in pixels. 
        \param height The height of the viewport in pixels. 
       */
      inline void setViewport(GLint x, GLint y, GLsizei width, GLsizei height)
      {
	std::array<GLint, 4> val = {{x, y, width, height}};
	setViewport(val);
      }

      /*! \brief Sets the viewport using the passed viewport state.
       */
      inline void setViewport(const std::array<GLint, 4>& val)
      {
	if (val == _viewPortState) return;

	_viewPortState = val;
	glViewport(_viewPortState[0], _viewPortState[1], _viewPortState[2], _viewPortState[3]);
	detail::errorCheck();
      }
      
      /*! \brief Returns the current viewport state.
       
        The returned array contains, in order, the leftmost pixel, the
        lowest pixel, the width and the height of the viewport.
       */
      inline const std::array<GLint, 4>& getViewport() const
      { return _viewPortState; }

      /*! \brief Swaps the front and back buffers.
       
        This command performs a glutSwapBuffers() and then executes
        any tasks left in the OpenGL task list. These tasks might have
        arisen from host program communication or some other
        asynchronous communication.
       */
      inline void swapBuffers()
      {
	glutSwapBuffers();
	detail::errorCheck();
	_glTasks.drainQueue();
	++_frameCounter;
      }

      /*! \brief Add a task to be performed after the next \ref swapBuffers.
       
        This function is used to allow other threads to instruct the
        OpenGL render thread to perform some task. This is usually
        used when a simulation thread wishes to update some data used
        for rendering.
       */
      inline void queueTask(std::function<void()> threadfunc)
      { _glTasks.queueTask(threadfunc); }

      /*! \brief The total number of \ref swapBuffers() calls.
	
	This function should give a the count of the number of frames
	drawn to the screen, assuming \ref swapBuffers() is used to
	paint the back buffer to the screen.
       */
      inline size_t getFrameCount() const 
      { return _frameCounter; }

      bool testExtension(std::string extension)
      {
	if (_extensions.empty())
	  {
	    detail::errorCheck();
	    size_t numExtensions = detail::glGet<GL_NUM_EXTENSIONS>();
	    for (size_t i(0); i < numExtensions; ++i)
	      _extensions.insert(std::string(reinterpret_cast<const char*>(glGetStringi(GL_EXTENSIONS, i))));
	  }

	return _extensions.find(extension) != _extensions.end();
      }

      inline void setBlend(bool newstate) { testAndSetState(newstate, _blend, GL_BLEND); }
      inline void setAlphaTest(bool newstate) { testAndSetState(newstate, _alphaTest, GL_ALPHA_TEST); }
      inline void setDepthTest(bool newstate) { testAndSetState(newstate, _depthTest, GL_DEPTH_TEST); }
      inline void setCullFace(bool newstate) { testAndSetState(newstate, _cullFace, GL_CULL_FACE); }

#ifndef GL_SAMPLE_SHADING
# define GL_SAMPLE_SHADING GL_SAMPLE_SHADING_ARB
#endif
      inline void setSampleShading(bool newstate) { testAndSetState(newstate, _sampleShading, GL_SAMPLE_SHADING); }

      void bindDefaultVAO() {
	glBindVertexArray(_dummyVAO);
	detail::errorCheck();
      }

    protected:
      inline void testAndSetState(const bool newstate, bool& oldstate, const GLenum cap)
      {
#ifdef MAGNET_DEBUG
	if (glIsEnabled(cap) != oldstate)
	  M_throw() << "Something is altering the GL state outside of Magnet!";
#endif

	if (newstate==oldstate) return;

	if (newstate) 
	  glEnable(cap);
	else 
	  glDisable(cap);
	oldstate = newstate;
      }

      /*! \brief State tracking variable for blending */
      bool _blend;
      bool _alphaTest;
      bool _depthTest;
      bool _cullFace;
      bool _sampleShading;


      /*! \brief A dummy vertex array object (VAO).
	
	The OpenGL standard requires from 3.2 (core), a VAO must be
	bound to the state as the default VAO is not available.
       */
      GLuint _dummyVAO;

      /*! \brief A counter of the number of calls to \ref
          swapBuffers(). 
      */
      size_t _frameCounter;

#ifdef MAGNET_CLGL
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
#endif

      std::set<std::string> _extensions;

      friend class shader::detail::Shader;

      /*! \brief The stack of bound shaders.
       */
      std::vector<shader::detail::Shader*> _shaderStack;

      /*! \brief A cache of the current OpenGL viewport state.
       */
      std::array<GLint, 4> _viewPortState;

      /*! \brief A queue of tasks to complete in the GL thread.
	
        These tasks are issued after the next \ref swapBuffers
        function call.
       */
      magnet::thread::TaskQueue _glTasks;

#ifdef MAGNET_CLGL
      /*! \brief If a matching OpenCL context does not exist, it will
        create one from the current OpenGL context.
       */
      inline void initCL()
      {
	if (_clInitialised) return;
	_clInitialised = true;

	initOpenCLContext();
	_clcommandQ = cl::CommandQueue(_clcontext, _cldevice);
      }

      /*! \brief Initializes an OpenCL context, platform and device
        from the current OpenGL context.
       */
      inline void initOpenCLContext()
      {
	std::cout << "GL-Context " << _context << ": Creating an OpenCL GPU context" << std::endl;	
	if (initOpenCLContextType(CL_DEVICE_TYPE_GPU)) return;
	std::cout << "GL-Context " << _context << ": Failed to create OpenCL GPU context, trying all devices for OpenCL context" << std::endl;
	if (initOpenCLContextType(CL_DEVICE_TYPE_ALL)) return;

	M_throw() << "Failed to create an OpenCL context from the OpenGL context!";
      }

      /*! \brief Attempts to Initializes an OpenCL context given a
        specific device type.
       
        \returns True if a matching device and context could be
        constructed.
       */
      inline bool initOpenCLContextType(cl_device_type devType)
      {
	std::vector<cl::Platform> platforms;
        cl::Platform::get(&platforms);
	
	//Now cycle through the platforms trying to get a context with a GPU
	for (auto& platform : platforms)
	  {
	    std::cout << "GL-Context " << _context << ":   Trying OpenCL platform - " 
		      << platform.getInfo<CL_PLATFORM_VENDOR>() 
		      << " - " << platform.getInfo<CL_PLATFORM_NAME>()
		      << " - " << platform.getInfo<CL_PLATFORM_VERSION>()
		      << std::endl;	

	    std::vector<cl::Device> devices;
	    try {
	      platform.getDevices(devType, &devices);
	    } catch (...)
	      { continue; }
    
	    for (auto& device : devices)
	      {
		std::cout << "GL-Context " << _context << ":     Trying  Device - " 
			  << device.getInfo<CL_DEVICE_NAME>() 
			  << " - " << device.getInfo<CL_DRIVER_VERSION>() 
			  << std::endl;	
		if (!getCLGLContext(platform, device)) continue;

		std::cout << "GL-Context " << _context << ": Success" << std::endl;
		
		//Success! now set the platform+device and return!
		_clplatform = platform;
		_cldevice = device;
		return true;
	      }
	  }

	return false;
      }

      /*! \fn bool getCLGLContext(cl::Platform clplatform, cl::Device dev)
        
        \brief This is a system specific command to build a OpenCL
        context from the current OpenGL context.
       
        \returns True if a context could be created for the passed
        device and platform.
       */

      /*! \fn ContextKey getCurrentContextKey()
        
        \brief This is a system-specific command to fetch the
        system-specific OpenGL-context-handle of the current GL
        context.
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
#endif

#ifdef MAGNET_DEBUG
      static void DebugCallback(unsigned int source, unsigned int type,
				unsigned int id, unsigned int severity,
				int length, const char* message, void* userParam)
      {
	std::cerr << "GLError: Source=";
	switch (source)
	  {
	  case GL_DEBUG_SOURCE_API_ARB:
	    std::cerr << "OpenGL"; break;
	  case GL_DEBUG_SOURCE_WINDOW_SYSTEM_ARB:
	    std::cerr << "Window System"; break;
	  case GL_DEBUG_SOURCE_SHADER_COMPILER_ARB:
	    std::cerr << "Shader Compiler"; break;
	  case GL_DEBUG_SOURCE_THIRD_PARTY_ARB:
	    std::cerr << "Third Party"; break;
	  case GL_DEBUG_SOURCE_APPLICATION_ARB:
	    std::cerr << "Application"; break;
	  case GL_DEBUG_SOURCE_OTHER_ARB:
	    std::cerr << "Other"; break;
	  default:
	    std::cerr << "Unknown(" << source << ")"; break;
	  }
	
	std::cerr << ", Type=";

	switch (type)
	  {
	  case GL_DEBUG_TYPE_ERROR_ARB:
	    std::cerr << "Error"; break;
	  case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR_ARB:
	    std::cerr << "Deprecated behavior"; break;
	  case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR_ARB:
	    std::cerr << "Undefined behavior"; break;
	  case GL_DEBUG_TYPE_PORTABILITY_ARB:
	    std::cerr << "Portability"; break;
	  case GL_DEBUG_TYPE_PERFORMANCE_ARB:
	    std::cerr << "Performance"; break;
	  case GL_DEBUG_TYPE_OTHER_ARB:
	    std::cerr << "Other"; break;
	  default:
	    std::cerr << "Unknown(" << type << ")"; break;
	  }

	std::cerr << ", Severity=";

	switch (severity)
	  {
	  case GL_DEBUG_SEVERITY_HIGH_ARB:
	    std::cerr << "High"; break;
	  case GL_DEBUG_SEVERITY_MEDIUM_ARB:
	    std::cerr << "Medium"; break;
	  case GL_DEBUG_SEVERITY_LOW_ARB:
	    std::cerr << "Low"; break;
	  default:
	    std::cerr << "Unknown(" << severity; break;
	  }

	std::cerr << ", ID=" << id << ", Message=\"";
	if (length == -1)
	  std::cerr << message;
	else
	  std::cerr << std::string(message, message + length);
	  
	std::cerr << "\"\n";
	std::cerr << "Stack trace:" << magnet::stacktrace() << "\n";
      }
#endif

      /*! \brief Initializes the OpenGL context and state tracking.
       */
      inline void init()
      {
	_frameCounter = 0;
	_context = getCurrentContextKey();
	std::cout << "GL-Context " << _context << ": Created a new OpenGL context" << std::endl;	
	//////////Capability testing /////////////////////////////

	//Check for errors before running glew
	detail::errorCheck();
	//We must set this true as glew uses deprecated code
	//(glGetString) for testing for extensions
	glewExperimental = GL_TRUE;
	if (glewInit() != GLEW_OK)
	  M_throw() << "GL-Context " << _context << "Failed to initialise GLEW!";
	//Suppress GL errors generated by glew, by fetching and
	//discarding the last GL error;
	glGetError();

	std::cout << "GL-Context " << _context 
		  << ": OpenGL version " 
		  << detail::glGet<GL_MAJOR_VERSION>() << "."
		  << detail::glGet<GL_MINOR_VERSION>() << std::endl;	
#ifdef MAGNET_DEBUG
	if (testExtension("GL_ARB_debug_output"))
	  {
	    glDebugMessageCallbackARB(&Context::DebugCallback, NULL);
	    glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS_ARB);
	  }
#endif
	///////Variable initialisation ///////////////////////////
	_viewPortState = detail::glGet<GL_VIEWPORT>();

	_vertexAttributeState.resize(detail::glGet<GL_MAX_VERTEX_ATTRIBS>());
	for (GLuint i(1); i < _vertexAttributeState.size(); ++i)
	  {
	    glVertexAttrib4f(i, 0,0,0,1);
	    detail::errorCheck();
	  }	  

	color(0,1,1,1);
	resetInstanceTransform();

	//Bind the dummy VAO object (so hacky)
	glGenVertexArrays(1, &_dummyVAO);
	detail::errorCheck();
	glBindVertexArray(_dummyVAO);
	detail::errorCheck();
      }

      /////////////////////////////////////////////////////

      //! \brief Only allow \ref getContext() to instantiate Context objects.
      Context(): 
	_blend(false),
	_alphaTest(false),
	_depthTest(false),
	_cullFace(false),
	_sampleShading(false)
#ifdef MAGNET_CLGL
	,
	_clInitialised(false)
#endif
      {}

      //! \brief The system-dependent handle to the GL context.
      ContextKey _context;

      /*! \brief Class used to track the state of a vertex attribute
       * array.
       */
      struct VertexAttrState
      {
	VertexAttrState(): 
	  active(false),
	  divisor(0)
	{
	  current_value[0] = current_value[1] = current_value[2] = 0; 
	  current_value[3] = 1; 
	}

	bool active;
	std::array<GLfloat, 4> current_value;
	GLuint divisor;
      };

      /*! \brief The state of the vertex attributes */
      std::vector<VertexAttrState> _vertexAttributeState;
    };
  }
}
