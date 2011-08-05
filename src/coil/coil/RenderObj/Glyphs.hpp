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
#include <magnet/GL/objects/instanced.hpp>
#include <magnet/GL/primatives/Sphere.hpp>

namespace coil {  
  class Glyphs : public DataSetChild, public magnet::GL::objects::Instanced
  {
  public:
    inline Glyphs(std::string name, DataSet& ds): DataSetChild(name, ds) {}

    inline virtual void clTick(const magnet::GL::Camera& cam) {}

    inline virtual void glRender(magnet::GL::FBO& fbo, const magnet::GL::Camera& cam, RenderMode mode) 
    {
      //Do not allow a glRender if uninitialised
      if (!_primitiveVertices.size()) return;

      _primitiveVertices.getContext().resetInstanceTransform();
      _positionSel->bindAttribute();
      _scaleSel->bindAttribute();
      _colorSel->bindAttribute();
      _orientSel->bindAttribute();

      Instanced::glRender();
    }

    inline virtual void init(const std::tr1::shared_ptr<magnet::thread::TaskQueue>& systemQueue)
    {
      RenderObj::init(systemQueue);            
      //Initialise the Gtk controls
      _gtkOptList.reset(new Gtk::VBox);
      _gtkOptList->show();

      Gtk::HSeparator* separator;

      //Glyph selection and level of detail
      _glyphBox.reset(new Gtk::HBox); _glyphBox->show();
      
      {
	Gtk::Label* label = Gtk::manage(new Gtk::Label("Glyph Type")); label->show();
	_glyphBox->pack_start(*label, false, false, 5);
	
	_glyphType.reset(new Gtk::ComboBoxText); _glyphType->show();

	_glyphType->append_text("Sphere");
	_glyphType->append_text("Arrows");
	_glyphType->append_text("Cylinder");
	_glyphType->append_text("Rod");
	_glyphType->set_active(0);

	_glyphBox->pack_start(*_glyphType, false, false, 5);
	_glyphType->signal_changed()
	  .connect(sigc::mem_fun(*this, &Glyphs::glyph_type_changed));
      }
      
      {
	_glyphLOD.reset(new Gtk::SpinButton(1.0, 0)); _glyphLOD->show();
	_glyphLOD->get_adjustment()->configure(1, 1, 1, 1.0, 1.0, 0.0); //A temporary low LOD setting
	_glyphLOD->set_numeric(true);
	_glyphBox->pack_end(*_glyphLOD, false, false, 5);
	_glyphLOD->signal_value_changed()
	  .connect(sigc::mem_fun(*this, &Glyphs::glyph_LOD_changed));

	Gtk::Label* label = Gtk::manage(new Gtk::Label("Level of Detail")); label->show();
	_glyphBox->pack_end(*label, false, false, 5);
      }

      _gtkOptList->pack_start(*_glyphBox, false, false, 5);

      separator = Gtk::manage(new Gtk::HSeparator); separator->show(); _gtkOptList->pack_start(*separator, false, false, 0);

      //The attribute selectors
      _positionSel.reset(new AttributeSelector(magnet::GL::Context::instanceOriginAttrIndex,
					       false));
      
      _positionSel->buildEntries("Position Data Field:", _ds, 3, 3, Attribute::COORDINATE, 0);
      _gtkOptList->pack_start(*_positionSel, false, false);

      separator = Gtk::manage(new Gtk::HSeparator); separator->show(); _gtkOptList->pack_start(*separator, false, false, 0);
      _scaleSel.reset(new AttributeSelector(magnet::GL::Context::instanceScaleAttrIndex));

      _scaleSel->buildEntries("Scale Data Field:", _ds, 1, 4,
			      Attribute::INTENSIVE | Attribute::EXTENSIVE, 3);
      _gtkOptList->pack_start(*_scaleSel, false, false);

      separator = Gtk::manage(new Gtk::HSeparator); separator->show(); _gtkOptList->pack_start(*separator, false, false, 0);
      _colorSel.reset(new AttributeColorSelector);
      _colorSel->buildEntries("Color Data Field:", _ds, 1, 4, 
			      Attribute::INTENSIVE | Attribute::EXTENSIVE, 4);
      _gtkOptList->pack_start(*_colorSel, false, false);

      separator = Gtk::manage(new Gtk::HSeparator); separator->show(); _gtkOptList->pack_start(*separator, false, false, 0);

      _orientSel.reset(new AttributeOrientationSelector);
      _orientSel->buildEntries("Orientation Data Field:", _ds, 3, 4, 
			       Attribute::INTENSIVE | Attribute::EXTENSIVE, 4);
      _gtkOptList->pack_start(*_orientSel, false, false);

      glyph_type_changed();
    }
    
    inline virtual void deinit()
    {
      Instanced::deinit();
      RenderObj::deinit();
      _gtkOptList.reset();
      _positionSel.reset();
      _scaleSel.reset(); 
      _colorSel.reset();
      _orientSel.reset();
      _glyphType.reset();
      _glyphLOD.reset();
    }

    inline virtual void showControls(Gtk::ScrolledWindow* win)
    {
      win->remove();
      _gtkOptList->unparent();
      win->add(*_gtkOptList);
      win->show();
    }


    inline virtual Gtk::TreeModel::iterator addViewRows(RenderObjectsGtkTreeView& view, 
							Gtk::TreeModel::iterator& parent_iter)
    {
      Gtk::TreeModel::iterator iter = view._store->append(parent_iter->children());
      (*iter)[view._columns->m_name] = getName();
      (*iter)[view._columns->m_visible] = visible();
      (*iter)[view._columns->m_shadowcasting] = shadowCasting();
      (*iter)[view._columns->m_obj] = this;
      return iter;
    }

  protected:
    inline void glyph_type_changed()
    {
      int type = _glyphType->get_active_row_number();
      switch (type)
	{
	case 0: //Spheres
	  _glyphLOD->get_adjustment()->configure(3, 0.0, 8.0, 1.0, 1.0, 0.0);
	  break;
	case 1: //Arrows
	case 2: //Cylinder
	case 3: //Rod
	default:
	  _glyphLOD->get_adjustment()->configure(6, 3.0, 32.0, 1.0, 5.0, 0.0);
	  break;
	}
    }
    
    inline void glyph_LOD_changed() { Instanced::init(_ds.size()); }

    inline virtual magnet::GL::element_type::Enum  getElementType()
    { return magnet::GL::element_type::TRIANGLES; }
    
    inline virtual std::vector<GLfloat> getPrimitiveVertices()
    {
      int LOD = _glyphLOD->get_value_as_int();
      int type = _glyphType->get_active_row_number();

      switch (type)
	{
	case 0: //Spheres
	  {
	    magnet::GL::primatives::Sphere sph(magnet::GL::primatives::Sphere::tetrahedron, LOD);
	    return std::vector<GLfloat>(sph.getVertices(), sph.getVertices() + sph.getVertexCount() * 3);
	  }
	case 1: //Arrows
	case 2: //Cylinder
	case 3: //Rod
	default:
	  {
	    std::vector<GLfloat> vertices(2 * LOD * 3);

	    for (int vert = 0; vert < 2 * LOD; ++vert)
	      {
		vertices[3 * vert + 0] = 0.5f * std::sin((vert / 2) * 2.0f * M_PI / LOD);
		vertices[3 * vert + 1] = 0.5f * std::cos((vert / 2) * 2.0f * M_PI / LOD);
		vertices[3 * vert + 2] = vert % 2;
	      }
	    return vertices;
	  }
	}
    }
    
    inline virtual std::vector<GLfloat> getPrimitiveNormals()
    {
      int LOD = _glyphLOD->get_value_as_int();
      int type = _glyphType->get_active_row_number();

      switch (type)
	{
	case 0: //Spheres
	  {
	    magnet::GL::primatives::Sphere sph(magnet::GL::primatives::Sphere::tetrahedron, LOD);
	    return std::vector<GLfloat>(sph.getVertices(), sph.getVertices() + sph.getVertexCount() * 3);
	  }
	case 1: //Arrows
	case 2: //Cylinder
	case 3: //Rod
	default:
	  {
	    std::vector<GLfloat> normals(2 * LOD * 3);

	    for (int vert = 0; vert < 2 * LOD; ++vert)
	      {
		GLfloat x = 0.5f * std::sin((vert / 2) * 2.0f * M_PI / LOD);
		GLfloat y = 0.5f * std::cos((vert / 2) * 2.0f * M_PI / LOD);
		GLfloat scale = 1.0f / std::sqrt(x * x + y * y);
		normals[3 * vert + 0] = x * scale;
		normals[3 * vert + 1] = y * scale;
		normals[3 * vert + 2] = 0;
	      }
	    return normals;
	  }
	}
    }
    
    inline virtual std::vector<GLuint>  getPrimitiveIndicies()
    {
      int LOD = _glyphLOD->get_value_as_int();
      int type = _glyphType->get_active_row_number();

      switch (type)
	{
	case 0: //Spheres
	  {
	    magnet::GL::primatives::Sphere sph(magnet::GL::primatives::Sphere::tetrahedron, LOD);
	    return std::vector<GLuint>(sph.getFaces(), sph.getFaces() + sph.getFaceCount() * 3);
	  }
	case 1: //Arrows
	case 2: //Cylinder
	case 3: //Rod
	default:
	  {
	    //2 triangles per face (3 indices per triangle)
	    std::vector<GLuint> indices(6 * LOD);      
	    for (int vert = 0; vert < LOD; ++vert)
	      {
		indices[6 * vert + 0] = (2 * vert + 0) % (2 * LOD);
		indices[6 * vert + 1] = (2 * vert + 1) % (2 * LOD);
		indices[6 * vert + 2] = (2 * vert + 2) % (2 * LOD);
		indices[6 * vert + 3] = (2 * vert + 1) % (2 * LOD);
		indices[6 * vert + 4] = (2 * vert + 3) % (2 * LOD);
		indices[6 * vert + 5] = (2 * vert + 2) % (2 * LOD);
	      }
	    return indices;
	  }
	}
    }
    
    std::auto_ptr<Gtk::VBox> _gtkOptList;
    std::auto_ptr<AttributeSelector> _positionSel;
    std::auto_ptr<AttributeSelector> _scaleSel; 
    std::auto_ptr<AttributeColorSelector> _colorSel;
    std::auto_ptr<AttributeOrientationSelector> _orientSel;
    std::auto_ptr<Gtk::ComboBoxText> _glyphType;
    std::auto_ptr<Gtk::SpinButton> _glyphLOD;
    std::auto_ptr<Gtk::HBox> _glyphBox;
 };
}
