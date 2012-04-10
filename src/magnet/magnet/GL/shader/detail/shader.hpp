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
#include <magnet/GL/matrix.hpp>
#include <magnet/exception.hpp>
#include <magnet/string/formatcode.hpp>
#include <magnet/string/line_number.hpp>
#include <boost/any.hpp>
#include <tr1/array>
#include <tr1/unordered_map>
#include <string>
#include <cstring>

#define STRINGIFY(A) #A

namespace magnet {
  namespace GL {
    namespace shader {
      namespace detail {
	class Shader;
	/*! \brief This class is used to store the assigned value of a
	  shader uniform and facilitate updating or retrieving shader
	  uniform values.
	 
	  The stored value is used to optimise redundant assignments
	  of shader uniforms and to allow fast, type-safe access to
	  the currently assigned uniform value.
	 
	  As shader uniforms may have several types, we must store the
	  type information ourselves. All standard uniform types
	  passed to a shader can be reduced into either an array of
	  GLfloat's or GLint's. This class reduces the data to a
	  std::tr1::array of floats or ints and places the data in a
	  boost::any containter.
	 
	  This class is returned from \ref Shader::operator[]() calls
	  to handle type based assignments of the shader.
	 
	  Please do not copy, reference or store this class in any
	  way, it does not track the currently bound program object
	  and so it should only be returned as a temporary from the
	  \ref Shader::operator[]() calls.
	 */
	class ShaderUniformValue
	{
	  friend class Shader;
	  /*! \brief Set the uniform handle corresponding to this class.
	   
	    Only the Shader is allowed to update this value.
	   */
	  inline void setHandle(GLint uniformHandle) { _uniformHandle = uniformHandle; }
	public:	 
	  inline ShaderUniformValue(): _uniformHandle(-1) {}

	  /*! \brief Test the current value of the uniform.*/
	  template<class T>
	  bool operator==(const T& val) const
	  {
	    return ((!_data.empty())
		    && (typeid(T) == _data.type())
		    && (boost::any_cast<T>(_data) == val));
	  }

	  /*! \brief Retrieve the current value of the uniform.
	   
	    This function can only return std::tr1::array types! All
	    values passed to the shader are converted to either
	    std::tr1::array<GLfloat,X> or std::tr1::array<GLint,X>
	    (where X is the number of elements) before being passed to
	    this class. So you must fetch them back in exactly this
	    form.
	   */
	  template<class T> const T as()
	  {
	    if (_data.empty()) M_throw() << "Uniform hasn't been assigned yet! Cannot retrieve its value";
	    if (typeid(T) != _data.type()) M_throw() << "Invalid as() cast for uniform value";
	    return boost::any_cast<T>(_data);
	  }

	  inline void operator=(const GLint& val)
	  { if (test_assign(val)) uniform(1, 1, &val); }

	  inline void operator=(const GLfloat& val)
	  { if (test_assign(val)) uniform(1, 1, &val); }

	  inline void operator=(const std::tr1::array<GLfloat, 1>& val)
	  { if (test_assign(val)) uniform(1, 1, &(val[0])); }

	  inline void operator=(const std::tr1::array<GLfloat, 2>& val)
	  { if (test_assign(val)) uniform(2, 1,&(val[0])); }

	  inline void operator=(const std::tr1::array<GLfloat, 3>& val)
	  { if (test_assign(val)) uniform(3, 1,&(val[0])); }

	  inline void operator=(const std::tr1::array<GLfloat, 4>& val)
	  { if (test_assign(val)) uniform(4, 1,&(val[0])); }

	  inline void operator=(const std::tr1::array<GLint, 1>& val)
	  { if (test_assign(val)) uniform(1, 1,&(val[0])); }

	  inline void operator=(const std::tr1::array<GLint, 2>& val)
	  { if (test_assign(val)) uniform(2, 1,&(val[0])); }

	  inline void operator=(const std::tr1::array<GLint, 3>& val)
	  { if (test_assign(val)) uniform(3, 1,&(val[0])); }

	  inline void operator=(const std::tr1::array<GLint, 4>& val)
	  { if (test_assign(val)) uniform(4, 1,&(val[0])); }

	  inline void operator=(const GLMatrix& val)
	  { 
	    if (test_assign(val)) 
	      { 
		glUniformMatrix4fv(_uniformHandle, 1, GL_FALSE, &(val[0])); 
		GL::detail::errorCheck(); 
	      }
	  }
	
	  inline void operator=(const math::Matrix& mat)
	  { 
	    std::tr1::array<GLfloat, 9> val;
	    for (size_t i(0); i < 3 * 3; ++i)
	      val[i] = mat(i);
	  
	    if (test_assign(mat))
	      { 
		glUniformMatrix3fv(_uniformHandle, 1, GL_FALSE, &(val[0])); 
		GL::detail::errorCheck(); 
	      }
	  }

	  inline void operator=(const math::Vector& vec)
	  { 
	    if (test_assign(vec))
	      {
		std::tr1::array<GLfloat, 3> val = {{vec[0], vec[1], vec[2]}};
		uniform(3, 1, &(val[0]));
	      }
	  }

	  inline void operator=(const std::vector<math::Vector>& val)
	  {
	    if (val.empty())
	      M_throw() << "Cannot assign a uniform from an empty vector";
	    
	    std::vector<GLfloat> data;
	    data.resize(val.size()*3);
	    
	    for (size_t i(0); i < val.size(); ++i)
	      for (size_t j(0); j < 3; ++j)
		data[i * 3 + j] = val[i][j];
	    
	    if (test_assign(val)) uniform(3, val.size(), &data[0]);
	  }

	  /**@}*/

	private:
	  inline void uniform(size_t width, size_t count, const GLfloat* val) 
	  { 
	    switch (width)
	      {
	      case 1: glUniform1fv(_uniformHandle, count, val); break;
	      case 2: glUniform2fv(_uniformHandle, count, val); break;
	      case 3: glUniform3fv(_uniformHandle, count, val); break;
	      case 4: glUniform4fv(_uniformHandle, count, val); break;
	      default: M_throw() << "Invalid uniform width";
	      }
	    GL::detail::errorCheck(); 
	  }

	  inline void uniform(size_t width, size_t count, const GLint* val)
	  { 
	    switch (width)
	      {
	      case 1: glUniform1iv(_uniformHandle, count, val); break;
	      case 2: glUniform2iv(_uniformHandle, count, val); break;
	      case 3: glUniform3iv(_uniformHandle, count, val); break;
	      case 4: glUniform4iv(_uniformHandle, count, val); break;
	      default: M_throw() << "Invalid uniform width";
	      }
	    GL::detail::errorCheck(); 
	  }

	  /*! \brief Returns true if (val != current value), and
	    updates the cached value of the uniform.
	   
	    This function is used to test if an update of the uniform
	    is actually required, and if it is it updates the cached
	    value before returning true.
	   */
	  template<class T> 
	  bool test_assign(const T& val)
	  {
	    //If this uniform does not exist in the code, don't ever
	    //try to assign it
	    if (_uniformHandle == -1) return false;
	    
	    //If we're in debug mode we always reset the value of the uniform
#ifndef MAGNET_DEBUG
	    //Check if the value is already set to the new value
	    if (*this == val) return false;
#endif
	    
	    _data = boost::any(val);
	    return true;
	  }

	  GLint _uniformHandle;
	  boost::any _data;
	};

	/*! \brief An object to store the values of defines.
	 */
	class ShaderDefineValue
	{
	  std::string _value;
	  bool _needsRecompilation;

	  friend class Shader;

	  bool checkForRecompilation()
	  { 
	    bool retval = _needsRecompilation;
	    _needsRecompilation = false;
	    return retval;
	  }

	public:
	  ShaderDefineValue(): _needsRecompilation(false) {}

	  /*! \brief Test the current value of the define.*/
	  template<class T>
	  bool operator==(const T& val) const
	  { return !_value.compare(boost::lexical_cast<std::string>(val)); }

	  template<class T>
	  inline void operator=(const T& val)
	  { 
	    if (*this == val) return;
	    _value = boost::lexical_cast<std::string>(val);
	    _needsRecompilation = true;
	  }

	  operator std::string() const { return _value; }
	};

	/*! \brief A OpenGL shader object.
	 
	  This class maintains the GL objects associated to a
	  complete shader program, including the vertex, fragment and
	  geometry shaders. After the shaders have been \ref built(),
	  the shader can be \ref attach() ed, or \ref deinit() ed to
	  release the associated GL resources.
	 
	  The shader source can be changed at any point, and if the
	  shader is already built, it will be recompiled.
	 
	  The shader performs some cacheing of the uniform locations
	  and the uniform values, so you may redundantly assign
	  values to the shader uniforms without an additional OpenGL
	  library call cost.
	 
	 
	  There are several default bindings for attributes in the
	  shader. These default bindings (indices from 0 to 6) may be
	  used by your shader, but be warned that they are used by
	  the GL Context as aliases for some common state variables.
	 
	  The table of indices is as follows:
	  \li "vPosition" = \ref Context::vertexPositionAttrIndex
	  \li "vColor" = \ref Context::vertexColorAttrIndex
	  \li "vNormal" = \ref Context::vertexNormalAttrIndex
	  \li "iOrigin" = \ref Context::instanceOriginAttrIndex
	  \li "iOrientation" = \ref Context::instanceOrientationAttrIndex
	  \li "iScale" = \ref Context::instanceScaleAttrIndex
	  \li "vTexCoord" = \ref Context::vertexTexCoordAttrIndex
	 
	 */
	class Shader 
	{
	private:
	  Shader(const Shader& );
	  void operator=(const Shader& );
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
		GL::detail::errorCheck();
	      }

	    _uniformCache.clear();
	    _built = false;
	    _context.reset();
	  }

	  /*! \brief Attach the shader, so it is used for the next
	    rendering of OpenGL objects.
	   
	    This function optimises away redundant attach() calls, and
	    updates the GL Context to mark the shader as attached.
	   */
	  inline void attach() 
	  {
	    if (!_built) M_throw() << "Cannot attach a Shader which has not been built()";

	    typedef std::tr1::unordered_map<std::string, ShaderDefineValue>::iterator it;
	    bool rebuild = false;
	    for (it iPtr = _defineCache.begin(); iPtr != _defineCache.end(); ++iPtr)
	      if (iPtr->second.checkForRecompilation())
		rebuild = true;

	    if (rebuild)
	      {
		deinit();
		build();
	      }
	    
	    _context->_shaderStack.push_back(this);

	    glUseProgramObjectARB(_programHandle);
	    GL::detail::errorCheck();
	  }

	  inline void detach()
	  {
	    //This is to help catch assymetric attach/detach calls or
	    //other objects messing with the shader stack
	    if (_context->_shaderStack.empty())
	      M_throw() << "Detaching a shader from an empty shader stack!";

	    if (_context->_shaderStack.back() != this)
	      M_throw() << "Error, detaching a shader which is not the current shader!";
	    _context->_shaderStack.pop_back();
	    if (_context->_shaderStack.empty())
	      glUseProgramObjectARB(0);
	    else
	      glUseProgramObjectARB(_context->_shaderStack.back()->_programHandle);
	    GL::detail::errorCheck();
	  }

	  /*! \brief Used to set and retrieve values of \ref Shader
              uniforms (AKA shader arguments).
	   
	    This function lets you assign values to uniforms easily:
	    \code Shader A;
	    A.build();
	    //Assign a single integer uniform value
	    A["intShaderVariable"] = 1;
	    //Assign a vec3 uniform
	    std::tr1::array<GLfloat,3> vec3val
	    A["vec3ShaderVariable"] = vec3val; 
	    \endcode
	   
	    You may also retrieve the value of uniforms:
	    \code
	    std::tr1::array<GLint, 1> value = A["ShaderVariable"].as<std::tr1::array<GLint, 1> >();
	    std::tr1::array<GLfloat, 3> value = A["vec3ShaderVariable"].as<std::tr1::array<GLfloat, 3> >(); 
	    \endcode
	   
	    Due to the way the cached value is stored, ALL variable types
	    must be returned using the above format (wrapped in a std::tr1::array).
	   
	    \param uniformName The name of the uniform to assign a
	    value to.
	   
	    \return A ShaderUniformValue object representing a
	    uniform.
	   */
	  inline ShaderUniformValue& operator[](const std::string& uniformName)
	  {
	    if (!_built) M_throw() << "Cannot set the uniforms of a shader which has not been built()";

	    if (_context->_shaderStack.back() != this)
	      M_throw() << "You must attach() a shader before you can change/read its uniform's values";

	    //In non-debug mode, we cache the uniform address,
	    //otherwise we always redetermine it in case a program
	    //like gDebugger has changed it.
#ifndef MAGNET_DEBUG
	    std::tr1::unordered_map<std::string, ShaderUniformValue>::iterator it 
	      = _uniformCache.find(uniformName);

	    if (it != _uniformCache.end()) return it->second;
#endif

	    GLint uniformHandle = glGetUniformLocationARB(_programHandle, uniformName.c_str());
	    GL::detail::errorCheck();
#ifdef MAGNET_DEBUG
	    if (uniformHandle == -1)
	      std::cerr << "\nMAGNET WARNING: Uniform " << uniformName 
			<< " not found in this shader, returning dummy uniform\n";
#endif	      
	    _uniformCache[uniformName].setHandle(uniformHandle);

	    return _uniformCache[uniformName];
	  }

	  /*! \brief Used to set and retrieve values of \ref Shader
	    defines.
	    
	    \param defineName The name of the define to assign a value
	    to.
	    
	    \return A ShaderUniformValue object representing a
	    uniform.
	  */
	  inline ShaderDefineValue& defines(const std::string& defineName)
	  { return _defineCache[defineName]; }

	  /*! \brief Builds the shader and allocates the associated
	    OpenGL objects.
	    
	    This function will throw an exception if compilation
	    fails.
	  */
	  inline void build()
	  {
	    _context = Context::getContext();

	    if (_vertexShaderCode.empty())
	      _vertexShaderCode = magnet::string::format_code(initVertexShaderSource());
	    if (_fragmentShaderCode.empty())
	      _fragmentShaderCode = magnet::string::format_code(initFragmentShaderSource());
	    if (_geometryShaderCode.empty())
	      _geometryShaderCode = magnet::string::format_code(initGeometryShaderSource());
	  
	    GLint result;

	    _programHandle = glCreateProgramObjectARB();
	    GL::detail::errorCheck();

	    std::string defines = genDefines();
	    
	    //Vertex shader
	    if (!_vertexShaderCode.empty())
	      {
		if (!_context->testExtension("GL_ARB_vertex_program"))
		  M_throw() << "GL-Context " << _context 
			    << ": Critical OpenGL dependency: Vertex programs are not supported";
		
		if (!_context->testExtension("GL_ARB_vertex_shader"))
		  M_throw() << "GL-Context " << _context 
			    << ": Critical OpenGL dependency: Vertex shaders are not supported";

		GLhandleARB _vertexShaderHandle = glCreateShaderObjectARB(GL_VERTEX_SHADER);
		GL::detail::errorCheck();
		if (!_vertexShaderHandle)
		  M_throw() << "Failed to create vertex shader handle";
		const GLcharARB* src[2] = {defines.c_str(), _vertexShaderCode.c_str()};
		glShaderSourceARB(_vertexShaderHandle, 2, src, NULL);	 
		GL::detail::errorCheck();
		glCompileShaderARB(_vertexShaderHandle);	  
		GL::detail::errorCheck();
		glGetObjectParameterivARB(_vertexShaderHandle, GL_OBJECT_COMPILE_STATUS_ARB, &result);
		GL::detail::errorCheck();
		if (!result)
		  M_throw() << "Vertex shader compilation failed, build log follows\n"
			    << getShaderBuildlog(_vertexShaderHandle)
			    << "\n Source code:\n"
			    << magnet::string::add_line_numbers(defines + _vertexShaderCode)
			    << "\n";

		glAttachObjectARB(_programHandle,_vertexShaderHandle);
		GL::detail::errorCheck();
		glDeleteShader(_vertexShaderHandle);
		GL::detail::errorCheck();
	      }

	    //Fragment shader
	    if (!_fragmentShaderCode.empty())
	      {	
		if (!_context->testExtension("GL_ARB_fragment_program"))
		  M_throw() << "GL-Context " << _context 
			    << ": Critical OpenGL dependency: Fragment programs are not supported";
		
		if (!_context->testExtension("GL_ARB_fragment_shader"))
		  M_throw() << "GL-Context " << _context 
			    << ": Critical OpenGL dependency: Fragment shaders are not supported";

		GLhandleARB _fragmentShaderHandle = glCreateShaderObjectARB(GL_FRAGMENT_SHADER);
		GL::detail::errorCheck();

		if (!_fragmentShaderHandle)
		  M_throw() << "Failed to create fragment shader handle";

		const GLcharARB* src[2] = {defines.c_str(), _fragmentShaderCode.c_str()};
		glShaderSourceARB(_fragmentShaderHandle, 2, src, NULL);
		GL::detail::errorCheck();
		glCompileShaderARB(_fragmentShaderHandle);
		GL::detail::errorCheck();
		glGetObjectParameterivARB(_fragmentShaderHandle, GL_OBJECT_COMPILE_STATUS_ARB, &result);
		GL::detail::errorCheck();
		if (!result)
		  M_throw() << "Fragment shader compilation failed, build log follows\n"
			    << getShaderBuildlog(_fragmentShaderHandle)
			    << "\n Source code:\n"
			    << magnet::string::add_line_numbers(defines + _fragmentShaderCode)
			    << "\n";

		glAttachObjectARB(_programHandle,_fragmentShaderHandle);
		GL::detail::errorCheck();
		glDeleteShader(_fragmentShaderHandle);
		GL::detail::errorCheck();
	      }

	    //Geometry shader
	    if (!(_geometryShaderCode.empty()))
	      {
		if (!_context->testExtension("GL_EXT_geometry_shader4"))
		  M_throw() << "Geometry shaders are not supported by your OpenGL driver.";

		GLhandleARB _geometryShaderHandle = glCreateShaderObjectARB(GL_GEOMETRY_SHADER_EXT);
		GL::detail::errorCheck();

		if (!_geometryShaderHandle)
		  M_throw() << "Failed to create geometry shader handle";
		const GLcharARB* src[2] = {defines.c_str(), _geometryShaderCode.c_str()};
		glShaderSourceARB(_geometryShaderHandle, 2, src, NULL);
		GL::detail::errorCheck();
		glCompileShaderARB(_geometryShaderHandle);
		GL::detail::errorCheck();
		glGetObjectParameterivARB(_geometryShaderHandle, GL_OBJECT_COMPILE_STATUS_ARB, &result);
		GL::detail::errorCheck();
		if (!result)
		  M_throw() << "Geometry shader compilation failed, build log follows\n"
			    << getShaderBuildlog(_geometryShaderHandle)
			    << "\n Source code:\n"
			    << magnet::string::add_line_numbers(defines + _geometryShaderCode)
			    << "\n";

		glAttachObjectARB(_programHandle, _geometryShaderHandle);
		GL::detail::errorCheck();
		glDeleteShader(_geometryShaderHandle);
		GL::detail::errorCheck();
	      }
	    
	    //Bind the default shader variables to the indices
	    //specified in the \ref Context class.
	    glLinkProgramARB(_programHandle);
	    GL::detail::errorCheck();

	    //Done, now the inheriting shader should grab the
	    //locations of its uniforms
	    _built = true;
	  }
	  	  
	  //! \brief Fetch the vertex shader source code.
	  const std::string getVertexShaderSource() const 
	  { return _vertexShaderCode; }

	  /*! \brief Set vertex shader source code.
	    
	    If the shader has already been built, this will force a
	    recompilation of all the shaders source
	   */
	  void setVertexShaderSource(std::string src)
	  { _vertexShaderCode = src; if (_built) { deinit(); build(); } }

	  //! \brief Fetch the vertex shader source code.
	  const std::string getFragmentShaderSource() const 
	  { return _fragmentShaderCode; }

	  /*! \brief Set fragment shader source code.
	   
	    If the shader has already been built, this will force a
	    recompilation of all the shaders source
	   */
	  void setFragmentShaderSource(std::string src)
	  { _fragmentShaderCode = src; if (_built) { deinit(); build(); } }

	  //! \brief Fetch the vertex shader source code.
	  const std::string getGeometryShaderSource() const 
	  { return _geometryShaderCode; }

	  /*! \brief Set fragment shader source code.
	   
	    If the shader has already been built, this will force a
	    recompilation of all the shaders source
	   */
	  void setGeometryShaderSource(std::string src)
	  { _geometryShaderCode = src; if (_built) { deinit(); build(); } }

	protected:
	  GLhandleARB _programHandle;
	  bool _built;
	  Context::ContextPtr _context;

	  std::string _vertexShaderCode;
	  std::string _fragmentShaderCode;
	  std::string _geometryShaderCode;

	  std::tr1::unordered_map<std::string, ShaderUniformValue> _uniformCache;
	  std::tr1::unordered_map<std::string, ShaderDefineValue> _defineCache;

	  std::string genDefines()
	  {
	    std::ostringstream os;

	    typedef std::tr1::unordered_map<std::string, ShaderDefineValue>::iterator it;
	    for (it iPtr = _defineCache.begin(); iPtr != _defineCache.end(); ++iPtr)
	      os << "#define " << iPtr->first << " " << std::string(iPtr->second) << "\n";

	    return os.str();
	  }

	  /*! \brief Specifies the initial source of the geometry
	    shader.
	   
	    Derived \ref Shader classes only need to override this if
	    they want to specify a geometry shader.
	   */
	  virtual std::string initGeometryShaderSource() { return ""; }
	  
	  /*! \brief Specifies the initial source of the vertex
	    shader.
	   
	    Derived \ref Shader classes only need to override this if
	    they want a non-trivial vertex shader.
	   */
	  virtual std::string initVertexShaderSource() { return ""; }

	  /*! \brief Specifies the initial source of the fragment
	    shader.
	   
	    Every derived \ref Shader class needs to override this and
	    specify the fragment shader.
	   */
	  virtual std::string initFragmentShaderSource() { return ""; }
	
	  /*! \brief Fetches the build log for the passed shader
	    handle.
	   */
	  inline std::string getShaderBuildlog(GLhandleARB shaderHandle)
	  {
	    GLint errorLoglength;
	    glGetObjectParameterivARB(shaderHandle, GL_OBJECT_INFO_LOG_LENGTH_ARB, &errorLoglength);
	    GL::detail::errorCheck();
	    char* buffer = new char[errorLoglength];
	  
	    GLsizei actualErrorLogLength;
	    glGetInfoLogARB(shaderHandle, errorLoglength, &actualErrorLogLength, buffer);
	    GL::detail::errorCheck();

	    std::string retval(buffer, buffer + actualErrorLogLength);
	    delete[] buffer;
	    return retval;
	  }
       	};
      }
    }
  }
}

#undef STRINGIFY
