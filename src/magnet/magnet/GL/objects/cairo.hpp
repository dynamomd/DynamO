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
#include <magnet/GL/texture.hpp>
#include <magnet/GL/shader/detail/shader.hpp>
#include <cairommconfig.h>
#include <cairomm/context.h>
#include <cairomm/surface.h>

namespace magnet {
  namespace GL {
    namespace objects {
      /*! \brief A textured quad generated from cairo drawing commands.
       *
       * This class is used to form a base class for rendering cairo
       * surfaces into an OpenGL scene.
       */
      class CairoSurface
      {
	class CairoShader: public GL::shader::detail::Shader
	{
	  bool _alpha_testing;
	public:
	  CairoShader(): Shader(true, true), _alpha_testing(false) {}

	  void build(bool alpha_testing)
	  {
	    _alpha_testing = alpha_testing;
	    Shader::build();
	  }

#define STRINGIFY(A) #A
	  virtual std::string initVertexShaderSource()
	  {
	    std::ostringstream os;
	    os << "const bool ALPHA_TESTING = " << _alpha_testing << ";"
	       << STRINGIFY(
uniform mat4 ProjectionMatrix;
uniform mat4 ViewMatrix;

attribute vec4 vPosition;
attribute vec4 vColor;
attribute vec4 iOrigin;
attribute vec4 iOrientation;
attribute vec4 iScale;

varying vec2 texCoord;
varying vec4 color;

vec3 qrot(vec4 q, vec3 v)
{ return v + 2.0 * cross(cross(v,q.xyz) + q.w * v, q.xyz); }

void main()
{
  //Rotate the vertex according to the instance transformation, and
  //then move it to the instance origin.
  vec4 vVertex = ViewMatrix * vec4(qrot(iOrientation, vPosition.xyz * iScale.xyz)
				   + iOrigin.xyz, 1.0);
  gl_Position = ProjectionMatrix * vVertex;
  texCoord = 0.5 + 0.5 * vPosition.xy;
  if (ALPHA_TESTING)
    color = vColor;
});
	    return os.str();
	  }
	
	  virtual std::string initFragmentShaderSource()
	  {
	    std::ostringstream os;
	    os << "const bool ALPHA_TESTING = " << _alpha_testing << ";"
	       << STRINGIFY(
uniform sampler2D cairoTexture;
varying vec2 texCoord;
varying vec4 color;
void main() 
{ 
  if (ALPHA_TESTING)
    {
      if (texture2D(cairoTexture, texCoord).r > 0.5) discard;
      gl_FragColor = color; 
    }
  else
    gl_FragColor = texture2D(cairoTexture, texCoord); 
}); 
	    return os.str();
	  }
	};

#undef STRINGIFY
	
      public:
	//! \brief Destructor
	inline ~CairoSurface() { deinit(); }

	//! \brief Release any associated OpenGL resources.
	inline void deinit() 
	{ 
	  _cairoSurface.clear();
	  _cairoContext.clear();
	  _vertexData.deinit();
	  _shader.deinit();
	  _width = _height = 0;
	}

	/*! \brief Sets up the vertex buffer objects for the quad.
	 */
	inline void init(size_t width, size_t height, bool alpha_testing = true)
	{
	  _alpha_testing = alpha_testing;

	  _width = width;
	  _height = height;

	  { ///////////////////Vertex Data
	    // Single quad, in pre-transformed screen coordinates
	    // (rotate, translate and scale using the instancing
	    // mechanism)
	    std::vector<GLfloat> vertexdata(4 * 2);
	    vertexdata[2 * 0 + 0] = -1; vertexdata[2 * 0 + 1] = -1;
	    vertexdata[2 * 1 + 0] =  1; vertexdata[2 * 1 + 1] = -1;
	    vertexdata[2 * 2 + 0] =  1; vertexdata[2 * 2 + 1] =  1;
	    vertexdata[2 * 3 + 0] = -1; vertexdata[2 * 3 + 1] =  1;
	    _vertexData.init(vertexdata);
	    magnet::GL::Buffer<GLfloat> _vertexData;
	  }

	  _shader.build(_alpha_testing);
	  _surface.init(_width, _height);
	  
	  _cairoSurface = Cairo::ImageSurface::create(_alpha_testing ? Cairo::FORMAT_A8 : Cairo::FORMAT_ARGB32, 
						      _width, _height);
	  _cairoContext = Cairo::Context::create(_cairoSurface);
	}
	
	void redraw()
	{
	  //Clear the surface
	  _cairoContext->save();
	  _cairoContext->set_operator(Cairo::OPERATOR_SOURCE);
	  _cairoContext->set_source_rgba(0, 0, 0, 0);
	  _cairoContext->paint();
	  _cairoContext->set_operator(Cairo::OPERATOR_OVER);
	  _cairoContext->set_source_rgba(1, 0, 0, 1);
	  
	  drawCommands();
	  _cairoContext->restore();
	  
	  if (_alpha_testing)
	    signedDistanceTransform(_cairoSurface->get_data());

	  //Send the cairo surface to the GL texture
	  _surface.subImage(_cairoSurface->get_data(), _alpha_testing ? GL_RED : GL_BGRA, _width, _height);
	  _surface.parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	  _surface.parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	  _surface.parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	  _surface.parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	}
	
	/*! \brief Attaches the vertex buffer and renders the quad.
	 */
	inline void glRender()
	{
	  _vertexData.getContext().cleanupAttributeArrays();
	  _surface.bind(6);
	  _shader["cairoTexture"] = 6;
	  _shader.attach();
	  
	  _vertexData.getContext().color(0,0,0,1);
	  _vertexData.drawArray(magnet::GL::element_type::QUADS, 2); 
	}

      protected:
	virtual void drawCommands() 
	{
	  _cairoContext->scale(_width,_height);
	  _cairoContext->move_to(0.5,0.5);
	  _cairoContext->set_font_size(0.1);
	  _cairoContext->show_text("Hello!");
	}

	struct P: public std::tr1::array<int,2> 
	{ P() { (*this)[0] = (*this)[1] = -1; } };

	inline size_t iPos(size_t x, size_t y) { return y * _width + x; }

	inline void check(int x, int y, double delta, 
			  std::vector<P>& p, std::vector<double>& d, size_t i2)
	{
	  int i1 = iPos(x,y);
	  if (d[i1] + delta < d[i2]) 
	    {						      
	      p[i2] = p[i1];			              
	      const double deltaX = (p[i1][0] - x);
	      const double deltaY = (p[i1][1] - y);
	      d[i2] = std::sqrt(deltaX*deltaX + deltaY*deltaY);
	    }
	}

	void signedDistanceTransform(unsigned char* I)
	{
	  //The border of the image must be 0
	  for (size_t x(0); x < _width; ++x)
	    { I[iPos(x, 0)] = 0; I[iPos(x, _height - 1)] = 0; }
	  for (size_t y(0); y < _height; ++y)
	    { I[iPos(0, y)] = 0; I[iPos(_width - 1, y)] = 0; }
	  
	  double floatInf = std::numeric_limits<float>::max() / 2;

	  std::vector<P> p; p.resize(_width * _height);
	  std::vector<double> d; d.resize(_width * _height, floatInf);

	  //Iterate over the image's non-border pixels, finding color
	  //edge pixels.
	  for (size_t y(1); y < _height - 1; ++y)
	    for (size_t x(1); x < _width - 1; ++x)
	      if ((   I[iPos(x-1, y  )] != I[iPos(x, y)])
		  || (I[iPos(x+1, y  )] != I[iPos(x, y)])
		  || (I[iPos(x  , y-1)] != I[iPos(x, y)])
		  || (I[iPos(x  , y+1)] != I[iPos(x, y)]))
		{
		  d[iPos(x,y)] = 0;
		  p[iPos(x,y)][0] = x;
		  p[iPos(x,y)][1] = y;
		}

#define _check(X,Y,Delta)				      \
	  { int i1 = iPos((X),(Y));			      \
	  if (d[i1]+(Delta) < d[i2]) {			      \
	    p[i2] = p[i1];			              \
	    const double deltaX = (p[i1][0] - x);	      \
	    const double deltaY = (p[i1][1] - y);	      \
	    d[i2] = std::sqrt(deltaX*deltaX + deltaY*deltaY); \
	  }}

	  //First pass
	  for (int y(1); y < int(_height) - 1; ++y)
	    for (int x(1); x < int(_width) - 1; ++x)
	      {
		int i2 = iPos(x,y);
		_check( x-1, y-1, std::sqrt(2));
		_check( x,   y-1, 1);
		_check( x+1, y-1, std::sqrt(2));
		_check( x-1, y,   1);
	      }

	  //Second pass
	  for (int y(_height - 2); y > 0; --y)
	    for (int x(_width - 2); x > 0; --x)
	      {
		size_t i2 = iPos(x,y);
                _check(x + 1, y    ,            1);
                _check(x - 1, y + 1, std::sqrt(2));
                _check(x    , y + 1,            1);
                _check(x + 1, y + 1, std::sqrt(2));
	      }

#undef _check

	  //Transform the output to a scaled range
	  for (size_t x(0); x < _width; ++x)
	    for (size_t y(0); y < _height; ++y)
	      {
		size_t i = iPos(x,y);
		double locd = I[i] ? -d[i] : d[i];
		I[i] = std::min(255.0, std::max(0.0, 128 + locd));
	      }
	}

	Texture2D _surface;
	size_t _width;
	size_t _height;
	bool _alpha_testing;
	Cairo::RefPtr<Cairo::ImageSurface> _cairoSurface;
	Cairo::RefPtr<Cairo::Context> _cairoContext;
	magnet::GL::Buffer<GLfloat> _vertexData;
	CairoShader _shader;
      };
    }
  }
}
