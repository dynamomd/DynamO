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

#include <magnet/exception.hpp>

namespace magnet {
  //! \brief Namespace for all magnet OpenGL classes and functions.
  namespace GL {
    /*! \brief An OpenGL buffer object.
     *
     * This class is used to represent vertex/element/normal buffer
     * objects and provides some automatic memory handling for them.
     */
    class Buffer
    {
    public:
      //! \brief The available GL targets to which this Buffer may be
      //! bound.
      enum BufferBindTargets
	{ ARRAY = GL_ARRAY_BUFFER,
	  ELEMENT_ARRAY = GL_ELEMENT_ARRAY_BUFFER, 
	  PIXEL_PACK_BUFFER = GL_PIXEL_PACK_BUFFER,
	  PIXEL_UNPACK_BUFFER = GL_PIXEL_UNPACK_BUFFER
	};

      /*! \brief The possible host access patterns, if the host
       * accesses the data.
       */
      enum BufferUsage
	{ 
	  STREAM_DRAW = GL_STREAM_DRAW, 
	  STREAM_READ = GL_STREAM_READ, 
	  STREAM_COPY = GL_STREAM_COPY, 
	  STATIC_DRAW = GL_STATIC_DRAW, 
	  STATIC_READ = GL_STATIC_READ, 
	  STATIC_COPY = GL_STATIC_COPY, 
	  DYNAMIC_DRAW = GL_DYNAMIC_DRAW, 
	  DYNAMIC_READ = GL_DYNAMIC_READ, 
	  DYNAMIC_COPY = GL_DYNAMIC_COPY
	};

      inline Buffer(): _size(0) {}
      inline ~Buffer() { deinit(); }

      /*! \brief Initialises the Buffer object with some data.
       *
       *  This will create the underlying OpenGL buffer, and load it
       *  with the contents of data.
       *
       * \param data A vector containing the data to be loaded to the
       * Buffer.
       *
       * \param usage The expected host memory access pattern, used to
       * optimise performance.
       */
      template<class T>
      inline void init(const std::vector<T>& data, BufferUsage usage = STATIC_DRAW)
      {
	if (data.empty())
	  M_throw() << "Cannot initialise GL::Buffer with no data!";

	deinit();
	_size = data.size();

	glGenBuffersARB(1, &_buffer);	
	bind(ARRAY);
	glBufferData(ARRAY, _size * sizeof(T), &data[0], usage);
      }
      
      //! \brief Attach the Buffer to a OpenGL target
      inline void bind(BufferBindTargets target) const 
      {
	glBindBufferARB(target, _buffer);	
      }

      //! \brief Map a buffer onto the host device memory space;
      template<class T>
      inline T* map()
      {
	bind(ARRAY);
	return static_cast<T*>(glMapBuffer(GL_ARRAY_BUFFER, GL_READ_WRITE));
      }

      //! \brief Map a buffer onto the host device memory space;
      template<class T>
      inline const T* map() const
      {
	bind(ARRAY);
	return static_cast<const T*>(glMapBuffer(GL_ARRAY_BUFFER, GL_READ_ONLY));
      }

      inline void unmap() const
      {
	bind(ARRAY);
	glUnmapBuffer(GL_ARRAY_BUFFER);
      }

      /*! \brief Destroys any OpenGL resources associated with this
       * object.
       */
      inline void deinit() 
      {
	if (_size)
	  glDeleteBuffersARB(1, &_buffer);

	_size = 0;
      }

      /*!\brief Returns the size of the allocated buffer, or 0 if not
       * allocated.  
       *
       * This is the number of elements stored in the buffer, NOT the
       * byte size of the buffer. As type information is lost across
       * the GL interface, we must request the size in bytes using
       * \ref size_byte().
       */
      inline size_t size() const { return _size; }

      /*!\brief Returns the size in bytes of the allocated buffer, or
       * 0 if not allocated.
       */
      inline size_t byte_size() const 
      { 
	if (size())
	  {
	    bind(ARRAY);
	    ::GLint size;
	    glGetBufferParameteriv(GL_ARRAY_BUFFER, GL_BUFFER_SIZE, &size);
	    return size;
	  }
	else 
	  return 0;
      }

      /*! \brief Returns the underlying OpenGL handle for the
       * buffer */
      GLuint getGLObject() { return _buffer; }

    protected:
      size_t _size;
      GLuint _buffer;
    };
  }
}
