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
#include <ctime>

namespace magnet {
  namespace GL {
    namespace objects {
      /*! \brief A quad textured with a 2D image generated from cairo
        drawing commands.
       
        This class is used to form a base class for rendering cairo
        surfaces into an OpenGL scene.
       
        It also provides an alpha-tested magnification routine and
        corresponding shader to help fake "vectorised" bitmap
        graphics. The technique is briefly described in the paper
        "Improved Alpha-Tested Magnification for math::Vector Textures
        and Special Effects," by Chris Green from Valve.
       */
      class CairoSurface
      {
	/*! \brief An alpha-testing shader for painting Cario
	  generated textures.
	 */
	class CairoShader: public GL::shader::detail::Shader
	{
	  size_t _alpha_testing;
	public:
	  CairoShader(): _alpha_testing(0) {}

	  /*! \brief Builds the shader and sets the draw mode.
	   
	    \param alpha_testing Controls the mode of the shader,
	    current supported modes are:
	   
	    \li 0 : Standard texturing of the quad with the passed
	    texture.
	   
	    \li 1 : Use the red channel of the texture to perform
	    alpha testing for a value of r 0.5. The color of the
	    object is taken from the GL state.
	   */
	  void build(size_t alpha_testing)
	  {
	    _alpha_testing = alpha_testing;
	    Shader::build();
	  }

#define STRINGIFY(A) #A
	  virtual std::string initVertexShaderSource()
	  { 
	    return STRINGIFY(
layout (location = 1) in vec4 vColor;

out vec4 color;

void main()
{ color = vColor; }
);
	  }

	  virtual std::string initGeometryShaderSource()
	  {
	    std::ostringstream os;
	    os << "const int ALPHA_TESTING = " << _alpha_testing << ";"
	       << STRINGIFY(
layout(points) in;
layout(triangle_strip) out;
layout(max_vertices = 4) out;

uniform mat4 ProjectionMatrix;
uniform mat4 ViewMatrix;

in vec4 color[];

smooth out vec2 texCoord;
flat out vec4 frag_color;

void main()
{
  frag_color = color[0];

  //The y coordinate of the texture is flipped here as cairo provides
  //upsided down images

  texCoord = vec2(0.0, 1.0);
  gl_Position = ProjectionMatrix * (ViewMatrix * vec4(-1.0, -1.0, 0.0, 1.0));
  EmitVertex();

  texCoord = vec2(0.0, 0.0);
  gl_Position = ProjectionMatrix * (ViewMatrix * vec4(-1.0, +1.0, 0.0, 1.0));
  EmitVertex();

  texCoord = vec2(1.0, 1.0);
  gl_Position = ProjectionMatrix * (ViewMatrix * vec4(+1.0, -1.0, 0.0, 1.0));
  EmitVertex();

  texCoord = vec2(1.0, 0.0);
  gl_Position = ProjectionMatrix * (ViewMatrix * vec4(+1.0, +1.0, 0.0, 1.0));
  EmitVertex();

  EndPrimitive();
});
	    return os.str();
	  }
	
	  virtual std::string initFragmentShaderSource()
	  {
	    std::ostringstream os;
	    os << "const int ALPHA_TESTING = " << _alpha_testing << ";"
	       << STRINGIFY(
uniform sampler2D cairoTexture;
smooth in vec2 texCoord;
flat in vec4 frag_color;

layout (location = 0) out vec4 color_out;

void main() 
{ 
  if (ALPHA_TESTING > 0)
    {
      if (texture(cairoTexture, texCoord).r <= 0.5) discard;
      color_out = frag_color;
    }
  else
    {
      vec4 sample = texture(cairoTexture, texCoord);
      if (sample.a == 0.0) discard;
      color_out = sample;
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
	  the Cairo backend for rendering the texture.
	 
	  \param width The width of the final texture.
	 
	  \param height The height of the final texture.
	 
	  \param alpha_testing If alpha_testing > 0, this enables the
	  alpha-tested texture generation and sets the relative pixel
	  size of the Cario scene. For a value of 0, this class simply
	  renders a cairo scene and pastes it into an OpenGL
	  texture. See the class description for more general
	  information.
	 
	 */
	virtual void init(size_t width, size_t height, size_t alpha_testing = 0)
	{
	  deinit();
	  _alpha_testing = alpha_testing;

	  _width = width * (alpha_testing + !alpha_testing);
	  _height = height * (alpha_testing + !alpha_testing);

	  _shader.build(_alpha_testing);
	  _surface.init(_width / (_alpha_testing + !_alpha_testing), 
			_height / (_alpha_testing + !_alpha_testing));
	  _surface.parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	  _surface.parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	  _surface.parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	  _surface.parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	  
	  _cairoSurface = Cairo::ImageSurface::create(_alpha_testing ? Cairo::FORMAT_A8 : Cairo::FORMAT_ARGB32, 
						      _width, _height);
	  _cairoContext = Cairo::Context::create(_cairoSurface);
	  _pango = Pango::Layout::create(_cairoContext);
	  Pango::FontDescription font("sans 12");
	  _pango->set_font_description(font);
	  _cairoContext->save();
	}
	
	void clear()
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
	  _cairoContext->restore();
	}
	
	/* \brief Copies the cairo image to the OpenGL texture. */
	void syncCairoGL()
	{
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
	  
	  The position, orientation and size of the scene can be
	  controlled through the \ref Shader instance attributes.  Or
	  alternately through the passed modelview and projection
	  matrix.
	 */
	inline void glRender(const GLMatrix& projection = GLMatrix::identity(),
			     const GLMatrix& modelview = GLMatrix::identity())
	{
	  _shader.attach();
	  _surface.bind(6);
	  _shader["cairoTexture"] = 6;
	  _shader["ProjectionMatrix"] = projection;
	  _shader["ViewMatrix"] = modelview;
	  
	  glDrawArrays(element_type::POINTS, 0, 1);
	  
	  _shader.detach();
	}

      protected:

	Texture2D _surface;
	size_t _width;
	size_t _height;
	size_t _alpha_testing;
	Cairo::RefPtr<Cairo::ImageSurface> _cairoSurface;
	Cairo::RefPtr<Cairo::Context> _cairoContext;
	Glib::RefPtr<Pango::Layout> _pango;
	CairoShader _shader;
      };
    }
  }
}
