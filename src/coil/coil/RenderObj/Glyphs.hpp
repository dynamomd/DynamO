/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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

#include <coil/RenderObj/DataSet.hpp>
#include <coil/RenderObj/AttributeColorSelector.hpp>
#include <coil/RenderObj/AttributeOrientationSelector.hpp>
#include <magnet/GL/shader/sphere.hpp>
#include <magnet/GL/shader/cylinder.hpp>
#include <magnet/GL/shader/simple_render.hpp>
#include <magnet/GL/shader/render.hpp>
#include <magnet/GL/buffer.hpp>

#define STRINGIFY(A) #A

namespace coil {
  /*! \brief A shader for transform feedback.  
   */
  class DumbbellShader: public magnet::GL::shader::detail::Shader
  {
  public:
    DumbbellShader() {
      _tfVaryings = {"gl_Position", "g_color", "g_orientation", "g_scale"};
    }

    virtual std::string initVertexShaderSource()
    { 
      return STRINGIFY(
layout (location = 0) in vec4 vPosition;
layout (location = 1) in vec4 vColor;
layout (location = 4) in vec4 iOrientation;
layout (location = 5) in vec4 iScale;

out vec4 f_color;
out vec3 f_director;
out vec4 f_orientation;
out vec4 f_scale;

vec3 qrot(vec4 q, vec3 v)
{ return v + 2.0 * cross(q.xyz, cross(q.xyz, v) + q.w * v); }

//In this case, iScale we have (radius1, radius2, distance1, distance2)
void main()
{
  gl_Position = vPosition;
  f_color = vColor;
  f_director = qrot(iOrientation, vec3(0.0, 0.0, 1.0));
  f_orientation = iOrientation;
  f_scale = iScale + vec4(equal(iScale, vec4(0.0))) * iScale.x;
});
    }

    virtual std::string initGeometryShaderSource()
    {
      return STRINGIFY(
layout(points) in;
layout(points, max_vertices = 2) out;

in vec4 f_color[];
in vec3 f_director[];
in vec4 f_orientation[];
in vec4 f_scale[];

//vec4 gl_position is also collected
flat out vec4 g_color;
flat out vec4 g_orientation;
flat out float g_scale;

void main()
{  
  g_color = f_color[0];
  g_orientation = f_orientation[0];
  g_scale = f_scale[0].x;
  gl_Position = vec4(gl_in[0].gl_Position.xyz + f_scale[0].z * f_director[0], 0.0);
  EmitVertex();
  EndPrimitive();

  g_color = f_color[0];
  g_orientation = f_orientation[0];
  g_scale = f_scale[0].y;
  gl_Position = vec4(gl_in[0].gl_Position.xyz - f_scale[0].w * f_director[0], 0.0);
  EmitVertex();
  EndPrimitive();
});
	}

  };

  class Glyphs : public DataSetChild
  {
    enum GLYPH_TYPE
      {
	SPHERE_GLYPH=0,
	ARROW_GLYPH=1,
	CYLINDER_GLYPH=2,
	ROD_GLYPH=3,
	LINE_GLYPH=4,
	CUBE_GLYPH=5,
	DUMBBELL_GLYPH=6
      };

  public:
    Glyphs(magnet::GL::Context::ContextPtr context, std::string pointsName, DataSet& ds);

    inline ~Glyphs() { deinit(); }

    virtual void glRender(const magnet::GL::Camera&, RenderMode, const uint32_t);
    
    virtual void init(const std::shared_ptr<magnet::thread::TaskQueue>&);
    
    virtual void deinit();

    virtual void showControls(Gtk::ScrolledWindow* win);

    virtual Glib::RefPtr<Gdk::Pixbuf> getIcon();

    virtual uint32_t pickableObjectCount()
    { 
      if (visible())
	return _N * (2 * _xperiodicimages->get_value_as_int() + 1)
	  * (2 * _yperiodicimages->get_value_as_int() + 1)
	  * (2 * _zperiodicimages->get_value_as_int() + 1);
	
      return 0; 
    }

    virtual magnet::math::NVector<GLfloat, 4> getCursorPosition(uint32_t objID);

    virtual std::string getCursorText(uint32_t objID);

    virtual magnet::math::Vector getMaxCoord() const;
    virtual magnet::math::Vector getMinCoord() const;

  protected:
    void glyph_type_changed();
    void guiUpdate();

    virtual magnet::GL::element_type::Enum  getElementType();
    
    std::vector<GLfloat> getPrimitiveVertices();   
    std::vector<GLfloat> getPrimitiveNormals();
    std::vector<GLuint>  getPrimitiveIndicies();

    magnet::GL::Buffer<GLfloat> _primitiveVertices;
    magnet::GL::Buffer<GLfloat> _primitiveNormals;
    magnet::GL::Buffer<GLuint>  _primitiveIndices;

    std::unique_ptr<Gtk::VBox> _gtkOptList;
    std::unique_ptr<AttributeSelector> _scaleSel; 
    std::unique_ptr<AttributeColorSelector> _colorSel;
    std::unique_ptr<AttributeOrientationSelector> _orientSel;
    std::unique_ptr<Gtk::ComboBoxText> _glyphType;
    std::unique_ptr<Gtk::SpinButton> _glyphLOD;
    std::unique_ptr<Gtk::HBox> _glyphBox;
    std::unique_ptr<Gtk::CheckButton> _glyphRaytrace;
    std::unique_ptr<Gtk::SpinButton> _xperiodicimages;
    std::unique_ptr<Gtk::SpinButton> _yperiodicimages;
    std::unique_ptr<Gtk::SpinButton> _zperiodicimages;

    std::unique_ptr<Gtk::CheckButton> _drawbillboards;

    std::unique_ptr<Gtk::HBox>  _scaleFactorBox;
    std::unique_ptr<Gtk::Label> _scaleLabel;
    std::unique_ptr<Gtk::Entry> _scaleFactor;
    
    size_t _N;
    float _scale;
    int _initGlyphType;
    std::string _pointsName;
    magnet::GL::Context::ContextPtr _context;
    magnet::GL::shader::RenderShader _renderShader;
    magnet::GL::shader::RenderVSMShader _renderVSMShader;
    magnet::GL::shader::SphereShader _sphereShader;
    magnet::GL::shader::SphereVSMShader _sphereVSMShader;
    magnet::GL::shader::CylinderShader _cylinderShader;
    magnet::GL::shader::CylinderVSMShader _cylinderVSMShader;
    magnet::GL::shader::SimpleRenderShader _simpleRenderShader;
    DumbbellShader _dumbbellShader;
 };
}

#undef STRINGIFY
