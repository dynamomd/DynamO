/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>

#include <magnet/GL/buffer.hpp>
#include <stdexcept>

namespace cl {
#if defined(__CL_ENABLE_EXCEPTIONS)
#define __ERR_STR(x) #x
#else
#define __ERR_STR(x) NULL
#endif // __CL_ENABLE_EXCEPTIONS

class GLBuffer : public Buffer
{
  magnet::GL::Buffer* _bufobj;

public:
  inline static bool& hostTransfers()
  {
    static bool _hostTransfer = false;
    return _hostTransfer;
  }

  GLBuffer(const Context& context,
	   cl_mem_flags flags,
	   magnet::GL::Buffer& buff,
	   cl_int* err = NULL
	   ):
    _bufobj(&buff)
  {
    if (hostTransfers())
      {
	if ((flags & CL_MEM_COPY_HOST_PTR) || (flags & CL_MEM_USE_HOST_PTR))
	  throw std::runtime_error("Cannot use CL_MEM_COPY_HOST_PTR/CL_MEM_USE_HOST_PTR on a host transfer GLBuffer");
	
        ::GLint size = _bufobj->byte_size();
        cl_int error;
        object_ = ::clCreateBuffer(context(), flags, size, NULL, &error);

        detail::errHandler(error, __CREATE_GL_BUFFER_ERR);
        if (err != NULL) *err = error;
      }
    else
      {
	cl_int error;
	object_ = ::clCreateFromGLBuffer(context(), flags, _bufobj->getGLObject(), &error);
	
	detail::errHandler(error, __CREATE_GL_BUFFER_ERR);
	if (err != NULL) *err = error;
      }
  }

  void acquire(CommandQueue cmdq, cl_int* err = NULL)
  {
    if (hostTransfers())
      {
	const void* glBufPointer = _bufobj->map<void>();

	::GLint size = _bufobj->byte_size();

	void* clBufPointer = cmdq.enqueueMapBuffer(*this, true, CL_MAP_WRITE, 0, size);
	
	memcpy(clBufPointer, glBufPointer, size);

	_bufobj->unmap();
	cmdq.enqueueUnmapMemObject(*this, (void*)clBufPointer);
      }
    else
      {
	cl_int error = ::clEnqueueAcquireGLObjects(cmdq(),
						   1,
						   &((*this)()),
						   0,
						   NULL,
						   NULL);
	
	detail::errHandler(error, __ENQUEUE_ACQUIRE_GL_ERR);
	
	if (err != NULL) {
	  *err = error;
	}
      }
  }

  void release(CommandQueue cmdq, cl_int* err = NULL)
  {
    if (hostTransfers())
      {
	void* glBufPointer = _bufobj->map<void>();

	::GLint size = _bufobj->byte_size();

	void* clBufPointer = cmdq.enqueueMapBuffer(*this, true, CL_MAP_READ, 0, size);
	
	memcpy(glBufPointer, clBufPointer, size);

	_bufobj->unmap();
	cmdq.enqueueUnmapMemObject(*this, (void*)clBufPointer);	
      }
    else
      {
	cl_int error = ::clEnqueueReleaseGLObjects(cmdq(),
						   1,
						   &(*this)(),
						   0,
						   NULL,
						   NULL);
	
	detail::errHandler(error, __ENQUEUE_RELEASE_GL_ERR);
	
	if (err != NULL) {
	  *err = error;
	}
      }
  }

  //! Default constructor; buffer is not valid at this point.
  GLBuffer(): _bufobj(NULL) {}

};
};
