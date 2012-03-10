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

#include <coil/RenderObj/DataSet.hpp>
#include <coil/RenderObj/AttributeColorSelector.hpp>
#include <coil/RenderObj/AttributeOrientationSelector.hpp>
#include <magnet/GL/objects/instanced.hpp>
#include <magnet/GL/shader/sphere.hpp>
#include <magnet/GL/objects/primitives/sphere.hpp>
#include <magnet/GL/objects/primitives/cylinder.hpp>
#include <magnet/GL/objects/primitives/arrow.hpp>


namespace coil {  
  class Glyphs : public DataSetChild, public magnet::GL::objects::Instanced
  {
  public:
    inline Glyphs(std::string name, DataSet& ds): DataSetChild(name, ds), _scale(1) {}

    inline virtual void clTick(const magnet::GL::Camera& cam) {}

    virtual void glRender(const magnet::GL::Camera&, RenderMode);
    
    virtual void init(const std::tr1::shared_ptr<magnet::thread::TaskQueue>&);
    
    virtual void deinit();

    virtual void showControls(Gtk::ScrolledWindow* win);

    virtual Glib::RefPtr<Gdk::Pixbuf> getIcon();

    virtual uint32_t pickableObjectCount()
    { return visible() ? _N : 0; }

    virtual void pickingRender(magnet::GL::FBO& fbo,
			       const magnet::GL::Camera& cam, 
			       const uint32_t offset);

  protected:
    void glyph_type_changed();
    void guiUpdate();

    inline virtual magnet::GL::element_type::Enum  getElementType()
    { return magnet::GL::element_type::TRIANGLES; }
    
    virtual std::vector<GLfloat> getPrimitiveVertices();   
    virtual std::vector<GLfloat> getPrimitiveNormals();
    virtual std::vector<GLuint>  getPrimitiveIndicies();

    std::auto_ptr<Gtk::VBox> _gtkOptList;
    std::auto_ptr<AttributeSelector> _scaleSel; 
    std::auto_ptr<AttributeColorSelector> _colorSel;
    std::auto_ptr<AttributeOrientationSelector> _orientSel;
    std::auto_ptr<Gtk::ComboBoxText> _glyphType;
    std::auto_ptr<Gtk::SpinButton> _glyphLOD;
    std::auto_ptr<Gtk::HBox> _glyphBox;
    std::auto_ptr<Gtk::CheckButton> _glyphRaytrace;

    std::auto_ptr<Gtk::HBox>  _scaleFactorBox;
    std::auto_ptr<Gtk::Label> _scaleLabel;
    std::auto_ptr<Gtk::Entry> _scaleFactor;
    
    bool _raytraceable;
    bool _GL_ARB_sample_shading;
    float _scale;
    magnet::GL::Context::ContextPtr _context;
    magnet::GL::shader::SphereShader<false> _sphereShader;
    magnet::GL::shader::SphereVSMShader<false> _sphereVSMShader;
 };
}
