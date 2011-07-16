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
#include <tr1/array>

namespace magnet {
  namespace GL {
    namespace detail {
#define GL_GET_ENUM_TYPE_TRAIT_FACTORY(F)		\
      F(GL_ACCUM_ALPHA_BITS, GLint, 1)			\
      F(GL_ACCUM_BLUE_BITS, GLint, 1)			\
      F(GL_ACCUM_GREEN_BITS, GLint, 1)			\
      F(GL_ACCUM_CLEAR_VALUE, GLfloat, 4)		\
      F(GL_ACCUM_RED_BITS, GLint, 1)			\
      F(GL_ACTIVE_TEXTURE, GLint, 1)			\
      F(GL_ALIASED_POINT_SIZE_RANGE, GLint, 2)		\
      F(GL_ALIASED_LINE_WIDTH_RANGE, GLint, 2)		\
      F(GL_ALPHA_BIAS, GLfloat, 1)			\
      F(GL_ALPHA_BITS, GLint, 1)			\
      F(GL_ALPHA_SCALE, GLfloat, 1)			\
      F(GL_ALPHA_TEST, GLboolean, 1)			\
      F(GL_ALPHA_TEST_FUNC, GLenum, 1)			\
      F(GL_ALPHA_TEST_REF, GLfloat, 1)			\
      F(GL_ARRAY_BUFFER_BINDING, GLint, 1)		\
      F(GL_ATTRIB_STACK_DEPTH, GLint, 1)		\
      F(GL_AUTO_NORMAL, GLboolean, 1)			\
      F(GL_AUX_BUFFERS, GLint, 1)			\
      F(GL_BLEND, GLboolean, 1)				\
      F(GL_BLEND_COLOR, GLfloat, 4)			\
      F(GL_BLEND_DST_ALPHA, GLint, 1)			\
      F(GL_BLEND_DST_RGB, GLint, 1)			\
      							\
      							\
      F(GL_MAX_VERTEX_ATTRIBS, GLint, 1)

      template<size_t width, class T> struct return_type { typedef std::tr1::array<T, width> Type; };
      template<class T> struct return_type<1, T> { typedef T Type; };

      
      template<GLenum val> struct glget_enum_to_type {};
#define GL_GET_TO_C_TYPE(GL_ENUM_VAL, C_TYPE, WIDTH)		\
      template<> struct glget_enum_to_type<GL_ENUM_VAL>		\
      {								\
	typedef C_TYPE Type;					\
	typedef return_type<WIDTH, C_TYPE>::Type ReturnType;	\
      };
      
      GL_GET_ENUM_TYPE_TRAIT_FACTORY(GL_GET_TO_C_TYPE)
#undef GL_ENUM_TO_C_TYPE
#undef GL_GET_ENUM_TYPE_TRAIT_FACTORY

      inline void glGetWorker(GLenum val, GLboolean* ptr)
      {
	glGetBooleanv(val, ptr);
      }

      inline void glGetWorker(GLenum val, GLdouble* ptr)
      {
	glGetDoublev(val, ptr);
      }

      inline void glGetWorker(GLenum val, GLfloat* ptr)
      {
	glGetFloatv(val, ptr);
      }

      inline void glGetWorker(GLenum val, GLint* ptr)
      {
	glGetIntegerv(val, ptr);
      }

      /*! \brief A type safe glGet command, to fetch parameters of the
       * current OpenGL state.
       */
      template<GLenum ENUM>
      typename glget_enum_to_type<ENUM>::ReturnType glGet()
      {	
	typename glget_enum_to_type<ENUM>::ReturnType retval;
	glGetWorker(ENUM, reinterpret_cast<typename glget_enum_to_type<ENUM>::Type*>(&retval));
	return retval;
      }
    }
  }
}
