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

#define BASE_GL_TYPE_FACTORY(F)		\
  F(GL_BYTE, GLbyte)			\
  F(GL_UNSIGNED_BYTE, GLubyte)		\
  F(GL_SHORT, GLshort)			\
  F(GL_UNSIGNED_SHORT, GLushort)	\
  F(GL_INT, GLint)			\
  F(GL_UNSIGNED_INT, GLuint)		\
  F(GL_FLOAT, GLfloat)			\
  F(GL_DOUBLE, GLdouble)


namespace magnet {
  namespace GL {
    namespace detail {
      /*! \brief Type trait to convert GL type enums into actual GL C types.
       */
      template <GLenum T> struct gl_enum_to_c_type {};

#ifndef DOXYGEN_SHOULD_IGNORE_THIS

#define GL_ENUM_TO_C_TYPE(GL_ENUM_VAL, C_TYPE)				\
      template <> struct gl_enum_to_c_type<GL_ENUM_VAL> { typedef C_TYPE Type; };

      BASE_GL_TYPE_FACTORY(GL_ENUM_TO_C_TYPE)
#undef GL_ENUM_TO_C_TYPE


      /*! \brief Type trait to convert GL C types into GL enum values.
       */
      template <class T> struct c_type_to_gl_enum {};

#define C_TYPE_TO_GL_ENUM(GL_ENUM_VAL, C_TYPE)				\
      template <> struct c_type_to_gl_enum<C_TYPE> { static const GLenum val = GL_ENUM_VAL; };

      BASE_GL_TYPE_FACTORY(C_TYPE_TO_GL_ENUM)
#undef C_TYPE_TO_GL_ENUM

#endif
    }
  }
}
#undef BASE_GL_TYPE_FACTORY
