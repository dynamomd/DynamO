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
#include <magnet/image/signed_distance.hpp>
#include <cairommconfig.h>
#include <cairomm/context.h>
#include <cairomm/surface.h>

namespace magnet {
  namespace GL {
    namespace objects {
      /*! \brief A quad textured with a 2D image generated from cairo
       * drawing commands.
       *
       * This class is used to form a base class for rendering cairo
       * surfaces into an OpenGL scene.
       *
       * It also provides an alpha-tested magnification routine and
       * corresponding shader to help fake "vectorised" bitmap
       * graphics. The technique is briefly described in the paper
       * "Improved Alpha-Tested Magnification for math::Vector Textures and
       * Special Effects," by Chris Green from Valve.
       */
      class CairoSurface
      {
	/*! \brief An alpha-testing shader for painting Cario
	 * generated textures.
	 */
	class CairoShader: public GL::shader::detail::Shader
	{
	  size_t _alpha_testing;
	public:
	  CairoShader(): _alpha_testing(0) {}

	  /*! \brief Builds the shader and sets the draw mode.
	   *
	   * \param alpha_testing Controls the mode of the shader,
	   * current supported modes are:
	   *
	   * \li 0 : Standard texturing of the quad with the passed texture.
	   *
	   * \li 1 : Use the red channel of the texture to perform
	   * alpha testing for a value of r 0.5. The color of the
	   * object is taken from the GL state.
	   */
	  void build(size_t alpha_testing)
	  {
	    _alpha_testing = alpha_testing;
	    Shader::build();
	  }

#define STRINGIFY(A) #A
	  virtual std::string initVertexShaderSource()
	  {
	    std::ostringstream os;
	    os << "const int ALPHA_TESTING = " << _alpha_testing << ";"
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
  vec4 vVertex = ViewMatrix * vec4(qrot(iOrientation, vPosition.xyz * iScale.xyz)
				   + iOrigin.xyz, 1.0);
  gl_Position = ProjectionMatrix * vVertex;
  texCoord = 0.5 + 0.5 * vPosition.xy * vec2(1.0, -1.0);
  if (ALPHA_TESTING > 0) color = vColor;
});
	    return os.str();
	  }
	
	  virtual std::string initFragmentShaderSource()
	  {
	    std::ostringstream os;
	    os << "const int ALPHA_TESTING = " << _alpha_testing << ";"
	       << STRINGIFY(
uniform sampler2D cairoTexture;
varying vec2 texCoord;
varying vec4 color;
void main() 
{ 
  if (ALPHA_TESTING > 0)
    {
      if (texture2D(cairoTexture, texCoord).r <= 0.5) discard;
      gl_FragColor = color;
    }
  else
    {
      vec4 sample = texture2D(cairoTexture, texCoord);
      if (sample.a == 0) discard;
      gl_FragColor = sample;
    }
}); 
	    return os.str();
	  }
	};

#undef STRINGIFY
	
      public:
	//! \brief Destructor
	inline ~CairoSurface() { deinit(); }

	//! \brief Release any associated OpenGL resources.
	virtual void deinit() 
	{ 
	  _cairoSurface.clear();
	  _cairoContext.clear();
	  _surface.deinit();
	  _vertexData.deinit();
	  _shader.deinit();
	  _pango.clear();
	  _width = _height = 0;
	}

	/*! \brief Resizes the cairo texture if required.
	 */
	virtual void resize(size_t width, size_t height)
	{
	  if ((width == _width) && (height == _height))
	    return;

	  init(width, height, _alpha_testing);
	}

	/*! \brief Sets up the vertex buffer objects for the quad and
	 * the Cairo backend for rendering the texture.
	 *
	 * \param width The width of the final texture.
	 *
	 * \param height The height of the final texture.
	 *
	 * \param alpha_testing If alpha_testing > 0, this enables the alpha-tested
	 * texture generation and sets the relative pixel size of the
	 * Cario scene. For a value of 0, this class simply renders a
	 * cairo scene and pastes it into an OpenGL texture. See the class
	 * description for more general information.
	 *
	 */
	virtual void init(size_t width, size_t height, size_t alpha_testing = 0)
	{
	  deinit();
	  _alpha_testing = alpha_testing;

	  _width = width * (alpha_testing + !alpha_testing);
	  _height = height * (alpha_testing + !alpha_testing);

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
	  _surface.init(_width / (_alpha_testing + !_alpha_testing), 
			_height / (_alpha_testing + !_alpha_testing));
	  _surface.parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	  _surface.parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	  _surface.parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	  _surface.parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	  
	  _cairoSurface = Cairo::ImageSurface::create(_alpha_testing ? Cairo::FORMAT_A8 : Cairo::FORMAT_ARGB32, 
						      _width, _height);
	  _cairoContext = Cairo::Context::create(_cairoSurface);
	  _pango = Pango::Layout::create(_cairoContext);
	  Pango::FontDescription font("sans 12");
	  _pango->set_font_description(font);	 
	}
	
	/*! \brief Forces the underlying Cairo scene to be rerendered
	 * and the texture to be updated.
	 */
	void redraw()
	{
	  //Clear the surface
	  _cairoContext->save();
	  _cairoContext->set_operator(Cairo::OPERATOR_SOURCE);
	  //The clear alpha must be 0 for the alpha masking effect
	  _cairoContext->set_source_rgba(0, 0, 0, 0);
	  _cairoContext->paint();
	  //The draw alpha must be >0 for the alpha masking effect
	  _cairoContext->set_operator(Cairo::OPERATOR_OVER);
	  _cairoContext->set_source_rgba(1, 1, 1, 1);
	  
	  drawCommands();
	  _cairoContext->restore();
	  
	  //Send the cairo surface to the GL texture
	  if (_alpha_testing)
	    {
	      //Calculate the distance texture
	      image::SignedDistanceTransform(_cairoSurface->get_data(), _width, _height);

	      //Downsample to the actual size	      
	      size_t texXSize = _width / _alpha_testing;
	      size_t texYSize = _height / _alpha_testing;
	      std::vector<unsigned char> downsampled;
	      downsampled.resize(texXSize * texYSize);
	      for (size_t y(0); y < texYSize; ++y)
		for (size_t x(0); x < texXSize; ++x)
		  downsampled[y * texXSize + x] 
		    = _cairoSurface->get_data()[y * _alpha_testing * _width + x * _alpha_testing];
	      
	      //Send the data to the texture
	      _surface.subImage(downsampled, GL_RED);
	    }
	  else
	    _surface.subImage(_cairoSurface->get_data(), GL_BGRA, _width, _height);
	  
	}
	
	/*! \brief Renders the Cairo scene.
	 *
	 * The position, orientation and size of the scene can be
	 * controlled through the \ref Shader instance attributes.  Or
	 * alternately through the passed modelview and projection matrix.
	 */
	inline void glRender(const GLMatrix& projection = GLMatrix::identity(),
			     const GLMatrix& modelview = GLMatrix::identity())
	{
	  _shader.attach();
	  _surface.bind(6);
	  _shader["cairoTexture"] = 6;
	  _shader["ProjectionMatrix"] = projection;
	  _shader["ViewMatrix"] = modelview;
	  _vertexData.drawArray(magnet::GL::element_type::QUADS, 2);
	  _shader.detach();
	}

      protected:
	virtual void drawCommands() = 0;

	void drawCursor(double x, double y, double size)
	{
	  _cairoContext->move_to(x - 2 * size, y + 0.5 * size);
	  _cairoContext->line_to(x - 0.5 * size, y + 0.5 * size);
	  _cairoContext->line_to(x - 0.5 * size, y + 2 * size);
	  
	  _cairoContext->move_to(x + 0.5 * size, y - 2   * size);
	  _cairoContext->line_to(x + 0.5 * size, y - 0.5 * size);
	  _cairoContext->line_to(x + 2   * size, y - 0.5 * size);
	}

	void textBox(double x, double y, std::string text, double padding = 5)
	{
	  _pango->set_text(text.c_str());
	  //Fetch the text dimensions, use pango's built in extents
	  //calculator as it is very fast
	  double topleft[2] = {x, y};
	  int pango_width, pango_height;
	  _pango->get_size(pango_width, pango_height);
	  double bottomright[2];
	  bottomright[0] = topleft[0] + double(pango_width) / Pango::SCALE + 2 * padding;
	  bottomright[1] = topleft[1] + double(pango_height) / Pango::SCALE + 2 * padding;

	  //Make sure the box doesn't overlap the sides. The left hand
	  //side takes priority over the right
	  double dimensions[2] = {_width, _height};
	  for (size_t i(0); i < 2; ++i)
	    {
	      //right/bottom edge
	      double shift = std::min(0.0, dimensions[i] - bottomright[i]);
	      topleft[i] += shift;
	      bottomright[i] += shift;
	      //left/top edge
	      shift = std::max(-topleft[i], 0.0);
	      topleft[i] += shift;
	      bottomright[i] += shift;
	    }
	  
	  //Background box
	  _cairoContext->begin_new_path();
	  _cairoContext->rectangle(topleft[0], topleft[1], 
	  			   bottomright[0] - topleft[0],
	  			   bottomright[1] - topleft[1]);
	  
	  _cairoContext->set_source_rgba(0.5, 0.70588, 0.94118, 0.7);
	  _cairoContext->fill();
	  
	  //Main text
	  _cairoContext->set_source_rgba(0, 0, 0, 1);
	  _cairoContext->move_to(topleft[0] + padding, topleft[1] + padding);
	  _pango->show_in_cairo_context(_cairoContext);
	}

	Texture2D _surface;
	size_t _width;
	size_t _height;
	size_t _alpha_testing;
	Cairo::RefPtr<Cairo::ImageSurface> _cairoSurface;
	Cairo::RefPtr<Cairo::Context> _cairoContext;
	Glib::RefPtr<Pango::Layout> _pango;
	magnet::GL::Buffer<GLfloat> _vertexData;
	CairoShader _shader;
      };
    }
  }
}
