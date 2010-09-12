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

#include <string>
#include <magnet/detail/exception.hpp>

namespace magnet {
  namespace GL {
    namespace detail {
      
      /* This is a CRTP base class that builds shaders.
       *
       * It requires that the type that inherits it, specifies its own
       * type in the template parameter (T) and defines static member
       * functions called  T::vertexShaderSource() and T::fragmentShaderSource().
       */
      template<class T>
      class shader {

      protected:
	inline void build()
	{
	  std::string vertexShaderSrc = format_code(T::vertexShaderSource());
	  std::string fragmentShaderSrc = format_code(T::fragmentShaderSource());
	  
	  if (!(_vertexShaderHandle = glCreateShaderObjectARB(GL_VERTEX_SHADER)))
	    M_throw() << "Failed to create vertex shader handle";
	  
	  if (!(_fragmentShaderHandle = glCreateShaderObjectARB(GL_FRAGMENT_SHADER)))
	    M_throw() << "Failed to create fragment shader handle";
	  
	  GLint result;

	  const GLcharARB* src = vertexShaderSrc.c_str();
	  glShaderSourceARB(_vertexShaderHandle, 1, &src, NULL);	 
	  glCompileShaderARB(_vertexShaderHandle);	  
	  glGetObjectParameterivARB(_vertexShaderHandle, GL_OBJECT_COMPILE_STATUS_ARB, &result);
	  if (!result)
	    M_throw() << "Vertex shader compilation failed, build log follows\n"
		      << getShaderBuildlog(_vertexShaderHandle);

	  src = fragmentShaderSrc.c_str();
	  glShaderSourceARB(_fragmentShaderHandle, 1, &src, NULL);
	  glCompileShaderARB(_fragmentShaderHandle);	  
	  glGetObjectParameterivARB(_fragmentShaderHandle, GL_OBJECT_COMPILE_STATUS_ARB, &result);

	  if (!result)
	    M_throw() << "Fragment shader compilation failed, build log follows\n"
		      << getShaderBuildlog(_fragmentShaderHandle);


	  //Now we've build both shaders, combine them into a program
	  _shaderID = glCreateProgramObjectARB();
	  glAttachObjectARB(_shaderID,_vertexShaderHandle);
	  glAttachObjectARB(_shaderID,_fragmentShaderHandle);
	  glLinkProgramARB(_shaderID);

	  //Done, now the inheriting shader should grab the locations of its uniforms
	}
	
	GLhandleARB _vertexShaderHandle;
	GLhandleARB _fragmentShaderHandle;
	GLhandleARB _shaderID;

	inline std::string getShaderBuildlog(GLhandleARB shaderHandle)
	{
	  std::string retVal;
	  GLint errorLoglength;
	  glGetObjectParameterivARB(shaderHandle, GL_OBJECT_INFO_LOG_LENGTH_ARB, &errorLoglength);
	  retVal.resize(errorLoglength);
	  
	  GLsizei actualErrorLogLength;
	  glGetInfoLogARB(shaderHandle, errorLoglength, &actualErrorLogLength, &retVal[0]);
	
	  return retVal;
	}
	
	/*! \brief Can search and replace elements in a std::string. */
	inline std::string 
	format_code(std::string in)
	{
	  return search_replace(in,";",";\n");
	}
	
	inline std::string 
	search_replace(std::string in, const std::string& from, const std::string& to)
	{
	  if (!in.empty())
	    {
	      std::string::size_type toLen = to.length();
	      std::string::size_type frLen = from.length();
	      std::string::size_type loc = 0;
	      
	      while (std::string::npos != (loc = in.find(from, loc)))
		{
		  in.replace(loc, frLen, to);
		  loc += toLen;
		  
		  if (loc >= in.length())
		    break;
		}
	    }
	  return in;
	}
	
      };
    }
  }
}
