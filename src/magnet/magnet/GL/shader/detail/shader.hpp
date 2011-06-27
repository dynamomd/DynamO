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

#define GL_GLEXT_PROTOTYPES
#include <GL/glew.h>
#include <magnet/exception.hpp>
#include <magnet/string/searchreplace.hpp>
#include <tr1/array>
#include <string>

#define STRINGIFY(A) #A

namespace magnet {
  namespace GL {
    namespace shader {
      namespace detail {
	/*! \brief A OpenGL shader object.
	 *!
	 *! 
	 */
	class Shader 
	{
	public:

	protected:
	  /*! \brief This class is returned from \ref Shader::operator[]
	   *! calls to handle type based assignments of the shader
	   */
	  class ShaderUniform
	  {
	  public:
	    inline ShaderUniform(std::string uniformName, GLhandleARB programHandle):
	      _programHandle(programHandle)
	    {
	      _uniformHandle = glGetUniformLocationARB(_programHandle, uniformName.c_str());
	      if (_uniformHandle == -1) M_throw() << "Uniform " << uniformName << " not found in this shader";
	    }

	    inline void operator=(GLfloat val)
	    { 
	      glUseProgramObjectARB(_programHandle);
	      glUniform1f(_uniformHandle, val);
	    }

	    inline void operator=(const GLint& val)
	    { 
	      glUseProgramObjectARB(_programHandle);
	      glUniform1i(_uniformHandle, val);
	    }

	    template<class T>
	    inline void operator=(const std::tr1::array<T, 1>& val)
	    { operator=(val[0]); }

	    inline void operator=(const std::tr1::array<GLfloat, 2>& val)
	    { 
	      glUseProgramObjectARB(_programHandle);
	      glUniform2fv(_uniformHandle, 1, &(val[0]));
	    }

	    inline void operator=(const std::tr1::array<GLfloat, 3>& val)
	    { 
	      glUseProgramObjectARB(_programHandle);
	      glUniform3fv(_uniformHandle, 1, &(val[0]));
	    }

	    inline void operator=(const std::tr1::array<GLfloat, 4>& val)
	    { 
	      glUseProgramObjectARB(_programHandle);
	      glUniform4fv(_uniformHandle, 1, &(val[0]));
	    }

	    inline void operator=(const std::tr1::array<GLint, 2>& val)
	    { 
	      glUseProgramObjectARB(_programHandle);
	      glUniform2iv(_uniformHandle, 1, &(val[0]));
	    }

	    inline void operator=(const std::tr1::array<GLint, 3>& val)
	    { 
	      glUseProgramObjectARB(_programHandle);
	      glUniform3iv(_uniformHandle, 1, &(val[0]));
	    }

	    inline void operator=(const std::tr1::array<GLint, 4>& val)
	    { 
	      glUseProgramObjectARB(_programHandle);
	      glUniform4iv(_uniformHandle, 1, &(val[0]));
	    }
	  
	  private:
	    GLint _uniformHandle;
	    GLhandleARB _programHandle;
	  };

	public:
	  inline Shader():_built(false) {}

	  inline ~Shader() { deinit(); }

	  inline void deinit()
	  {
	    if (_built)
	      {
		glDeleteProgram(_shaderID);
		glDeleteShader(_vertexShaderHandle);
		glDeleteShader(_fragmentShaderHandle);
	      }

	    _built = false;
	  }

	  inline void attach() { glUseProgramObjectARB(_shaderID); }

	  inline ShaderUniform operator[](std::string uniformName)
	  {
	    return ShaderUniform(uniformName, _shaderID);
	  }

	  inline void build()
	  {
	    //Check for the ability to use fragment and vertex shaders
	    if (!GLEW_ARB_fragment_program || !GLEW_ARB_vertex_program
		|| !GLEW_ARB_fragment_shader || !GLEW_ARB_vertex_shader)
	      M_throw() << "This OpenGL context/driver does not support GLSL (programmable shaders)";

	    std::string vertexShaderSrc = format_code(vertexShaderSource());
	    std::string fragmentShaderSrc = format_code(fragmentShaderSource());
	  
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
	    _built = true;
	  }

	  void drawScreenQuad()
	  {
	    //Save the matrix state
	    glMatrixMode(GL_PROJECTION);
	    glPushMatrix();
	    glLoadIdentity();
	  
	    glMatrixMode(GL_MODELVIEW);
	    glPushMatrix();
	    glLoadIdentity();
	  
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
	  
	    //Restore the matrix state
	    glMatrixMode(GL_PROJECTION);
	    glPopMatrix();
	  
	    glMatrixMode(GL_MODELVIEW);
	    glPopMatrix();
	  }

	protected:
	  GLhandleARB _vertexShaderHandle;
	  GLhandleARB _fragmentShaderHandle;
	  GLhandleARB _shaderID;

	  bool _built;

	  virtual std::string vertexShaderSource()
	  {
	    return STRINGIFY( 
void main()
{
  gl_Position = ftransform();
  gl_TexCoord[0] = gl_MultiTexCoord0;
});
	  }

	  virtual std::string fragmentShaderSource() = 0;
	
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
	  { return magnet::string::search_replace(in,";",";\n"); }
	};
      }
    }
  }
}

#undef STRINGIFY
