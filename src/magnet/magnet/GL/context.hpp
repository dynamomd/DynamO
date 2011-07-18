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
#include <magnet/GL/detail/typesafe_get.hpp>
#include <magnet/exception.hpp>
#include <magnet/GL/matrix.hpp>
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
	  contexts[key].init();
	  
	return contexts[key];
      }
      
      /** @name Vertex attribute array interface. */
      /**@{*/
      /*! \brief Enables a vertex attribute array index. */
      inline void enableAttributeArray(GLuint attrnum)
      {
	if (attrnum >= _vertexAttributeState.size())
	  M_throw() << "Attribute index out of range";
	glEnableVertexAttribArray(attrnum);
	_vertexAttributeState[attrnum].active = true;
      }

      /*! \brief Disable all active vertex attribute arrays. */
      inline void cleanupAttributeArrays()
      {
	resetInstanceTransform();

	for (GLuint i(0); i < _vertexAttributeState.size(); ++i)
	  if (_vertexAttributeState[i].active) 
	    { glDisableVertexAttribArray(i); _vertexAttributeState[i].active = false;}
      }

      /*! \brief Sets the value of a vertex attribute, if no attribute
       * array is bound.
       *
       * This function only sets the state if it has been updated.
       */
      inline void setAttribute(GLuint idx, GLfloat x = 0, GLfloat y = 0, GLfloat z = 0, GLfloat w = 0)
      {
	if (idx >= _vertexAttributeState.size())
	  M_throw() << "Attribute index out of range";

	std::tr1::array<GLfloat, 4> newval = {{x,y,z,w}};

#ifdef MAGNET_DEBUG
	std::tr1::array<GLfloat, 4> oldval;
	glGetVertexAttribfv(idx, GL_CURRENT_VERTEX_ATTRIB, &oldval[0]);
	if (oldval != _vertexAttributeState[idx].current_value)
	  M_throw() << "Vertex attribute state changed without using the GL context!";
#endif

	if (newval == _vertexAttributeState[idx].current_value) return;
	_vertexAttributeState[idx].current_value = newval;

	glVertexAttrib4f(idx, newval[0], newval[1], newval[2], newval[3]);
      }
      
      /*! \brief Sets the divisor of a vertex attribute.
       *
       * The divisor is used in instancing to set the rate at which
       * vertex attributes are incremented.
       */
      inline void setAttributeDivisor(GLuint idx, GLuint divisor)
      {
	if (idx >= _vertexAttributeState.size())
	  M_throw() << "Attribute index out of range";

	if (divisor && !GL_ARB_instanced_arrays)
	  M_throw() << "Cannot perform instanced vertex attributes, GL_ARB_instanced_arrays is not supported";

	if (divisor == _vertexAttributeState[idx].divisor) return;
	_vertexAttributeState[idx].divisor = divisor;
	glVertexAttribDivisorARB(idx, divisor);
      }


      /*! \brief The index of the automatically-indexed position
       * vertex attribute.  
       *
       * This index for the vertex position is used in the OpenGL
       * standard.
       * \sa detail::Shader::build()
       */
      static const GLuint vertexPositionAttrIndex = 0;

      /*! \brief The index of the automatically-indexed position
       * color attribute.  
       * \sa detail::Shader::build()
       */
      static const GLuint vertexColorAttrIndex = 1;

      /*! \brief The index of the automatically-indexed normal vertex
       * attribute.  
       * \sa detail::Shader::build()
       */
      static const GLuint vertexNormalAttrIndex = 2;

      /*! \brief The index of the automatically-indexed instance origin
       * vertex attribute.  
       * \sa detail::Shader::build()
       */
      static const GLuint instanceOriginAttrIndex = 3;

      /*! \brief The index of the automatically-indexed instance
       * orientation vertex attribute.
       * \sa detail::Shader::build()
       */
      static const GLuint instanceOrientationAttrIndex = 4;

      /*! \brief The index of the automatically-indexed instance
       * scale vertex attribute.
       * \sa detail::Shader::build()
       */
      static const GLuint instanceScaleAttrIndex = 5;

      /*! \brief Convenience function to set the vertex attribute
       * representing the color in a shader.
       *
       * This uses the \ref vertexColorAttrIndex value for the index
       * of the color attribute.
       */
      inline void color(GLfloat r = 0, GLfloat g = 0, GLfloat b = 0, GLfloat a = 1) 
      { setAttribute(vertexColorAttrIndex, r, g, b, a); }

      /*! \brief Resets the vertex attributes used in instancing to
       * avoid unintended transformations of the instanced object.
       */
      inline void resetInstanceTransform()
      {
	setAttribute(instanceOriginAttrIndex, 0, 0, 0, 0);
	setAttribute(instanceOrientationAttrIndex, 0, 0, 0, 1);
	setAttribute(instanceScaleAttrIndex, 1, 1, 1, 0);
      }
      /**@}*/

      /** @name The OpenGL matrix interface. */
      /**@{*/
      void setProjectionMatrix(const GLMatrix& mat)
      { 
	_projectionMatrix = mat;
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf(_projectionMatrix);
      }

      void setViewMatrix(const GLMatrix& mat)
      { 
	_modelViewMatrix = mat;
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixf(_modelViewMatrix);
      }

      const GLMatrix& getModelViewMatrix() { return _modelViewMatrix; }
      const GLMatrix& getProjectionMatrix() { return _projectionMatrix; }

      void setTextureMatrix(GLuint textureUnit, const GLMatrix& mat)
      { 
	glActiveTextureARB(GL_TEXTURE0 + textureUnit);
	glMatrixMode(GL_TEXTURE);
	glLoadIdentity();
	glMultMatrixf(mat);
      }
      /**@}*/

      /*! \brief Draws a fullscreen quad, with texture coordinates
       * set from the bottom left (0,0) to the top right (1,1).
       */
      void drawScreenQuad()
      {
	setViewMatrix(GLMatrix::identity());
	setProjectionMatrix(GLMatrix::identity());
	
	glBegin(GL_QUADS);
	glTexCoord2f(0.0f, 0.0f);
	glVertex2d(-1, -1);
	glTexCoord2f(1.0f, 0.0f);
	glVertex2d(1, -1);
	glTexCoord2f( 1.0f, 1.0f);
	glVertex2d(1, 1);
	glTexCoord2f(0.0f, 1.0f);
	glVertex2d(-1, 1);
	glEnd();
      }


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

      /*! \brief Initializes the OpenGL context and state tracking.
       */
      inline void init()
      {
	_context = getCurrentContextKey();
	_modelViewMatrix = GLMatrix::identity();
	_projectionMatrix = GLMatrix::identity();

	if (glewInit() != GLEW_OK)
	  M_throw() << "Failed to initialise GLEW!";

	_vertexAttributeState.resize(detail::glGet<GL_MAX_VERTEX_ATTRIBS>());

	for (GLuint i(0); i < _vertexAttributeState.size(); ++i)
	  glVertexAttrib4f(i, 0,0,0,1);
	
	color(0,1,1,1);
	resetInstanceTransform();
      }

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
	std::tr1::array<GLfloat, 4> current_value;
	GLuint divisor;
      };

      /*! \brief The state of the vertex attributes */
      std::vector<VertexAttrState> _vertexAttributeState;

      GLMatrix _modelViewMatrix;
      GLMatrix _projectionMatrix;
    };
  }
}
