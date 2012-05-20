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

#include <magnet/exception.hpp>

namespace magnet {
  namespace GL {
    namespace detail {
      inline void errorCheck()
      {
#ifdef MAGNET_DEBUG
	GLenum errcode = glGetError();
	switch (errcode)
	  {
	  case GL_NO_ERROR: return;
	  case GL_INVALID_ENUM: M_throw() << "glGetError() returned GL_INVALID_ENUM";
	  case GL_INVALID_VALUE: M_throw() << "glGetError() returned GL_INVALID_VALUE";
	  case GL_INVALID_OPERATION: M_throw() << "glGetError() returned GL_INVALID_OPERATION";
	  case GL_OUT_OF_MEMORY: M_throw() << "glGetError() returned GL_OUT_OF_MEMORY";
	  case GL_INVALID_FRAMEBUFFER_OPERATION: M_throw() << "glGetError() returned GL_INVALID_FRAMEBUFFER_OPERATION";
	  default: M_throw() << "glGetError() returned " << errcode;
	  }
#endif       
      }
      
    }
  }
}
