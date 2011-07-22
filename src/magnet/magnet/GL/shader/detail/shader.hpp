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
#include <magnet/GL/context.hpp>
#include <magnet/exception.hpp>
#include <magnet/string/formatcode.hpp>
#include <tr1/array>
#include <tr1/unordered_map>
#include <string>

#define STRINGIFY(A) #A

namespace magnet {
  namespace GL {
    namespace shader {
      class Shader;
      namespace detail {

	/*! \brief Type handling for \ref Shader uniform (AKA
	 *  argument) assignment.
	 * 
	 * This class is returned from \ref Shader::operator[]() calls
	 * to handle type based assignments of the shader.
	 *
	 * Please do not copy, reference or store this class in any
	 * way, it does not track the currently bound program object
	 * and so it should only be returned as a temporary from the
	 * \ref Shader::operator[]() calls.
	 */
	class ShaderUniform
	{
	public:
	  inline ShaderUniform(GLint uniformHandle):
	    _uniformHandle(uniformHandle)
	  {}

	  inline void operator=(GLfloat val)
	  { glUniform1f(_uniformHandle, val); }

	  inline void operator=(const GLint& val)
	  { glUniform1i(_uniformHandle, val); }

	  template<class T>
	  inline void operator=(const std::tr1::array<T, 1>& val)
	  { operator=(val[0]); }

	  inline void operator=(const std::tr1::array<GLfloat, 2>& val)
	  { glUniform2fv(_uniformHandle, 1, &(val[0])); }

	  inline void operator=(const std::tr1::array<GLfloat, 3>& val)
	  { glUniform3fv(_uniformHandle, 1, &(val[0])); }

	  inline void operator=(const std::tr1::array<GLfloat, 4>& val)
	  { glUniform4fv(_uniformHandle, 1, &(val[0])); }

	  inline void operator=(const std::tr1::array<GLint, 2>& val)
	  { glUniform2iv(_uniformHandle, 1, &(val[0])); }

	  inline void operator=(const std::tr1::array<GLint, 3>& val)
	  { glUniform3iv(_uniformHandle, 1, &(val[0])); }

	  inline void operator=(const std::tr1::array<GLint, 4>& val)
	  { glUniform4iv(_uniformHandle, 1, &(val[0])); }

	  inline void operator=(const std::tr1::array<GLfloat, 9>& val)
	  { glUniformMatrix3fv(_uniformHandle, 1, GL_FALSE, &(val[0])); }

	  inline void operator=(const std::tr1::array<GLfloat, 16>& val)
	  { glUniformMatrix4fv(_uniformHandle, 1, GL_FALSE, &(val[0])); }

	  inline void operator=(const Vector& vec)
	  { 
	    std::tr1::array<GLfloat, 3> val = {{vec[0], vec[1], vec[2]}};
	    operator=(val);
	  }

	  inline void operator=(const Matrix& mat)
	  { 
	    std::tr1::array<GLfloat, 9> val;
	    for (size_t i(0); i < 3 * 3; ++i)
	      val[i] = mat(i);
	    operator=(val);
	  }
	  
	private:
	  GLint _uniformHandle;
	};


	/*! \brief A OpenGL shader object.
	 *
	 * This class maintains the GL objects associated to a
	 * complete shader program, including the vertex, fragment and
	 * geometry shaders. After the shaders have been \ref built(),
	 * the shader can be \ref attach() ed, or \ref deinit() ed to
	 * release the associated GL resources.
	 *
	 * The shader source can be changed at any point, and if the
	 * shader is already built, it will be recompiled.
	 *
	 * The shader performs some cacheing of the uniform locations
	 * and the uniform values, so you may redundantly assign
	 * values to the shader uniforms without an additional OpenGL
	 * library call cost.
	 *
	 *
	 * There are several default bindings for attributes in the
	 * shader. These default bindings (indices from 0 to 6) may be
	 * used by your shader, but be warned that they are used by
	 * the GL Context as aliases for some common state variables.
	 *
	 * The table of indices is as follows:
	 * \li "vPosition" = \ref Context::vertexPositionAttrIndex
	 * \li "vColor" = \ref Context::vertexColorAttrIndex
	 * \li "vNormal" = \ref Context::vertexNormalAttrIndex
	 * \li "iOrigin" = \ref Context::instanceOriginAttrIndex
	 * \li "iOrientation" = \ref Context::instanceOrientationAttrIndex
	 * \li "iScale" = \ref Context::instanceScaleAttrIndex
	 * \li "vTexCoord" = \ref Context::vertexTexCoordAttrIndex
	 *
	 */
	class Shader 
	{
	public:
	  /*! \brief Constructor for Shader objects.
	   */
	  Shader(): _built(false) {}

	  //! \brief Destructor
	  inline ~Shader() { deinit(); }

	  //! \brief Cause the shader to release its OpenGL resources.
	  inline void deinit()
	  {
	    if (_built)
	      {
		glDeleteProgram(_programHandle);
		if (!(_vertexShaderCode.empty()))
		  glDeleteShader(_vertexShaderHandle);
		if (!(_fragmentShaderCode.empty()))
		  glDeleteShader(_fragmentShaderHandle);

		if (!(_geometryShaderCode.empty()))
		  glDeleteShader(_geometryShaderHandle);
	      }

	    _vertexShaderCode.clear();
	    _fragmentShaderCode.clear();
	    _geometryShaderCode.clear();
	    _built = false;
	  }

	  /*! \brief Attach the shader, so it is used for the next
	   *rendering of OpenGL objects.
	   *
	   * This function also registers the matrix callbacks if
	   * either the view or projection matricies are used inside
	   * the shader.
	   */
	  inline void attach() 
	  {
	    if (!_built) M_throw() << "Cannot attach a Shader which has not been built()";
	    _context->setShader(_programHandle);
	  }

	  /*! \brief Used to set values of \ref Shader uniforms (AKA
           *   arguments).
	   *
	   * This function lets you assign values to uniforms quite
	   * easily:
	   * \code Shader A;
	   * A.build();
	   * A["ShaderVariable"] = 1.0; \endcode
	   *
	   * \param uniformName The name of the uniform to assign a
	   * value to.
	   *
	   * \return A ShaderUniform object representing a uniform.
	   */
	  inline ShaderUniform operator[](std::string uniformName)
	  {
	    if (!_built) M_throw() << "Cannot set the uniforms of a shader which has not been built()";

	    _context->setShader(_programHandle);

	    std::tr1::unordered_map<std::string, GLuint>::const_iterator it = _uniformAddressCache.find(uniformName);	    
	    if (it != _uniformAddressCache.end()) return ShaderUniform(it->second);

	    GLint uniformHandle = glGetUniformLocationARB(_programHandle, uniformName.c_str());
	    if (uniformHandle == -1)
	      M_throw() << "Uniform " << uniformName << " not found in this shader";	    
	    _uniformAddressCache[uniformName] = uniformHandle;
	    return ShaderUniform(uniformHandle);
	  }

	  /*! \brief Builds the shader and allocates the associated
	    OpenGL objects.
	    *
	    * This function will throw an exception if compilation
	    * fails.
	    */
	  inline void build()
	  {
	    _context = &(Context::getContext());

	    if (_vertexShaderCode.empty()) 
	      _vertexShaderCode = magnet::string::format_code(initVertexShaderSource());
	    if (_fragmentShaderCode.empty()) 
	      _fragmentShaderCode = magnet::string::format_code(initFragmentShaderSource());
	    if (_geometryShaderCode.empty()) 
	      _geometryShaderCode = magnet::string::format_code(initGeometryShaderSource());
	  
	    GLint result;

	    //Vertex shader
	    if (!_vertexShaderCode.empty())
	      {
		if (!GLEW_ARB_vertex_shader)
		  M_throw() << "Vertex shaders are not supported by your OpenGL driver.";

		if (!(_vertexShaderHandle = glCreateShaderObjectARB(GL_VERTEX_SHADER)))
		  M_throw() << "Failed to create vertex shader handle";
		const GLcharARB* src = _vertexShaderCode.c_str();
		glShaderSourceARB(_vertexShaderHandle, 1, &src, NULL);	 
		glCompileShaderARB(_vertexShaderHandle);	  
		glGetObjectParameterivARB(_vertexShaderHandle, GL_OBJECT_COMPILE_STATUS_ARB, &result);
		if (!result)
		  M_throw() << "Vertex shader compilation failed, build log follows\n"
			    << getShaderBuildlog(_vertexShaderHandle);
	      }

	    //Fragment shader
	    if (!_fragmentShaderCode.empty())
	      {
		if (!GLEW_ARB_fragment_shader)
		  M_throw() << "Fragment shaders are not supported by your OpenGL driver.";

		if (!(_fragmentShaderHandle = glCreateShaderObjectARB(GL_FRAGMENT_SHADER)))
		  M_throw() << "Failed to create fragment shader handle";
		const GLcharARB* src = _fragmentShaderCode.c_str();
		glShaderSourceARB(_fragmentShaderHandle, 1, &src, NULL);
		glCompileShaderARB(_fragmentShaderHandle);	  
		glGetObjectParameterivARB(_fragmentShaderHandle, GL_OBJECT_COMPILE_STATUS_ARB, &result);
		if (!result)
		  M_throw() << "Fragment shader compilation failed, build log follows\n"
			    << getShaderBuildlog(_fragmentShaderHandle);
	      }

	    //Geometry shader
	    if (!(_geometryShaderCode.empty()))
	      {
		if (!GL_EXT_geometry_shader4)
		  M_throw() << "Geometry shaders are not supported by your OpenGL driver.";

		if (!(_geometryShaderHandle = glCreateShaderObjectARB(GL_GEOMETRY_SHADER_EXT)))
		  M_throw() << "Failed to create geometry shader handle";
		const GLcharARB* src = _fragmentShaderCode.c_str();
		glShaderSourceARB(_geometryShaderHandle, 1, &src, NULL);
		glCompileShaderARB(_geometryShaderHandle);
		glGetObjectParameterivARB(_fragmentShaderHandle, GL_OBJECT_COMPILE_STATUS_ARB, &result);
		if (!result)
		  M_throw() << "Geometry shader compilation failed, build log follows\n"
			    << getShaderBuildlog(_fragmentShaderHandle);
	      }

	    //Now we've build both shaders, combine them into a program
	    _programHandle = glCreateProgramObjectARB();

	    if (!(_vertexShaderCode.empty()))
	      glAttachObjectARB(_programHandle,_vertexShaderHandle);

	    if (!(_fragmentShaderCode.empty()))
	      glAttachObjectARB(_programHandle,_fragmentShaderHandle);

	    if (!(_geometryShaderCode.empty()))
	      glAttachObjectARB(_programHandle,_geometryShaderHandle);
	    
	    //Bind the default shader variables to the indices
	    //specified in the \ref Context class.
	    glBindAttribLocation(_programHandle, Context::vertexPositionAttrIndex, "vPosition");
	    glBindAttribLocation(_programHandle, Context::vertexColorAttrIndex, "vColor");
	    glBindAttribLocation(_programHandle, Context::vertexNormalAttrIndex, "vNormal");
	    glBindAttribLocation(_programHandle, Context::vertexTexCoordAttrIndex, "vTexCoord");
	    glBindAttribLocation(_programHandle, Context::instanceOriginAttrIndex, "iOrigin");
	    glBindAttribLocation(_programHandle, Context::instanceOrientationAttrIndex, "iOrientation");
	    glBindAttribLocation(_programHandle, Context::instanceScaleAttrIndex, "iScale");

	    glLinkProgramARB(_programHandle);

	    //Done, now the inheriting shader should grab the locations of its uniforms
	    _built = true;
	  }
	  	  
	  //! \brief Fetch the vertex shader source code.
	  const std::string getVertexShaderSource() const 
	  { return _vertexShaderCode; }

	  /*! \brief Set vertex shader source code.
	   *
	   * If the shader has already been built, this will force a
	   * recompilation of all the shaders source
	   */
	  void setVertexShaderSource(std::string src)
	  { _vertexShaderCode = src; if (_built) { deinit(); build(); } }

	  //! \brief Fetch the vertex shader source code.
	  const std::string getFragmentShaderSource() const 
	  { return _fragmentShaderCode; }

	  /*! \brief Set fragment shader source code.
	   *
	   * If the shader has already been built, this will force a
	   * recompilation of all the shaders source
	   */
	  void setFragmentShaderSource(std::string src)
	  { _fragmentShaderCode = src; if (_built) { deinit(); build(); } }

	  //! \brief Fetch the vertex shader source code.
	  const std::string getGeometryShaderSource() const 
	  { return _geometryShaderCode; }

	  /*! \brief Set fragment shader source code.
	   *
	   * If the shader has already been built, this will force a
	   * recompilation of all the shaders source
	   */
	  void setGeometryShaderSource(std::string src)
	  { _geometryShaderCode = src; if (_built) { deinit(); build(); } }

	protected:
	  GLhandleARB _vertexShaderHandle;
	  GLhandleARB _fragmentShaderHandle;
	  GLhandleARB _geometryShaderHandle;
	  GLhandleARB _programHandle;
	  bool _built;
	  Context* _context;

	  std::string _vertexShaderCode;
	  std::string _fragmentShaderCode;
	  std::string _geometryShaderCode;

	  std::tr1::unordered_map<std::string, GLuint> _uniformAddressCache;

	  /*! \brief Specifies the initial source of the geometry
	   * shader.
	   *
	   * Derived \ref Shader classes only need to override this if
	   * they want to specify a geometry shader.
	   */
	  virtual std::string initGeometryShaderSource() { return ""; }
	  
	  /*! \brief Specifies the initial source of the vertex
	   * shader.
	   *
	   * Derived \ref Shader classes only need to override this if
	   * they want a non-trivial vertex shader.
	   */
	  virtual std::string initVertexShaderSource() { return ""; }

	  /*! \brief Specifies the initial source of the fragment
	   * shader.
	   *
	   * Every derived \ref Shader class needs to override this
	   * and specify the fragment shader.
	   */
	  virtual std::string initFragmentShaderSource() { return ""; }
	
	  /*! \brief Fetches the build log for the passed shader
	   * handle.
	   */
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
       	};
      }
    }
  }
}

#undef STRINGIFY
