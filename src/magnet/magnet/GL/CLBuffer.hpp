/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

    This program is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    version 3 as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#pragma once

#define GL_GLEXT_PROTOTYPES
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>
#include <GL/glext.h>
#include <CL/cl.hpp>
#include <stdexcept>

namespace cl {
/*! \class Buffer
 * \brief Memory buffer interface.
 */

#if defined(__CL_ENABLE_EXCEPTIONS)
#define __ERR_STR(x) #x
#else
#define __ERR_STR(x) NULL
#endif // __CL_ENABLE_EXCEPTIONS


class GLBuffer : public Buffer
{
  bool _hostTransfer;
  ::GLuint _bufobj;

  ::GLenum _bufType;

public:
  GLBuffer(const Context& context,
	   cl_mem_flags flags,
	   ::GLuint bufobj,
	   ::GLenum bufType,
	   bool hostTransfer,
	   cl_int* err = NULL
	   ):
    _hostTransfer(hostTransfer),
    _bufobj(bufobj),
    _bufType(bufType)
  {
    if (_hostTransfer)
      {
	if ((flags & CL_MEM_COPY_HOST_PTR) || (flags & CL_MEM_USE_HOST_PTR))
	  throw std::runtime_error("Cannot use CL_MEM_COPY_HOST_PTR/CL_MEM_USE_HOST_PTR on a host transfer GLBuffer");

	glBindBuffer(_bufType, _bufobj);
        ::GLint size;
	glGetBufferParameteriv(_bufType, GL_BUFFER_SIZE, &size);

        cl_int error;
        object_ = ::clCreateBuffer(context(), flags, size, NULL, &error);

        detail::errHandler(error, __CREATE_GL_BUFFER_ERR);
        if (err != NULL) *err = error;
      }
    else
      {
	cl_int error;
	object_ = ::clCreateFromGLBuffer(context(), flags, bufobj, &error);
	
	detail::errHandler(error, __CREATE_GL_BUFFER_ERR);
	if (err != NULL) *err = error;
      }
  }

  Event acquire(CommandQueue & cmdq, cl_int* err = NULL)
  {
    Event retEvent;

    if (_hostTransfer)
      {
	glBindBuffer(_bufType, _bufobj);
	const void* glBufPointer = glMapBuffer(_bufType, GL_READ_ONLY);

	::GLint size;
	glGetBufferParameteriv(_bufType, GL_BUFFER_SIZE, &size);
	
	void* clBufPointer = cmdq.enqueueMapBuffer(*this, true, CL_MAP_WRITE, 0, size);
	
	memcpy(clBufPointer, glBufPointer, size);

	glUnmapBuffer(_bufType);
	cmdq.enqueueUnmapMemObject(*this, (void*)clBufPointer);
      }
    else
      {
	cl_int error = ::clEnqueueAcquireGLObjects(cmdq(),
						   1,
						   &(*this)(),
						   0,
						   0,
						   &retEvent());
	
	detail::errHandler(error, __ENQUEUE_ACQUIRE_GL_ERR);
	
	if (err != NULL) {
	  *err = error;
	}
      }
    return retEvent;
  }

  Event release(CommandQueue & cmdq, cl_int* err = NULL)
  {
    Event retEvent;

    if (_hostTransfer)
      {
	glBindBuffer(_bufType, _bufobj);
	void* glBufPointer = glMapBuffer(_bufType, GL_WRITE_ONLY);

	::GLint size;
	glGetBufferParameteriv(_bufType, GL_BUFFER_SIZE, &size);
	
	void* clBufPointer = cmdq.enqueueMapBuffer(*this, true, CL_MAP_READ, 0, size);

	
	memcpy(glBufPointer, clBufPointer, size);

	glUnmapBuffer(_bufType);
	cmdq.enqueueUnmapMemObject(*this, (void*)clBufPointer);	
      }
    else
      {
	cl_int error = ::clEnqueueReleaseGLObjects(cmdq(),
						   1,
						   &(*this)(),
						   0,
						   NULL,
						   &retEvent());
	
	detail::errHandler(error, __ENQUEUE_RELEASE_GL_ERR);
	
	if (err != NULL) {
	  *err = error;
	}
      }
    return retEvent;
  }

  //! Default constructor; buffer is not valid at this point.
  GLBuffer():_hostTransfer(false)  {}

  operator Buffer() {return *this;}
  
};
};
