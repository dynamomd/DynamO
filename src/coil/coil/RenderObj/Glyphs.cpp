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

#include <coil/RenderObj/Glyphs.hpp>
#include <magnet/GL/objects/primitives/sphere.hpp>
#include <magnet/GL/objects/primitives/cylinder.hpp>
#include <magnet/GL/objects/primitives/arrow.hpp>
#include <magnet/GL/objects/primitives/cube.hpp>
#include <coil/images/images.hpp>

namespace coil {  
  Glib::RefPtr<Gdk::Pixbuf> 
  Glyphs::getIcon()
  { return coil::images::Glyphs_Icon(); }

  magnet::GL::element_type::Enum  
  Glyphs::getElementType()
  { 
    switch (_glyphType->get_active_row_number())
      {
      case LINE_GLYPH:
	return magnet::GL::element_type::LINES;
      default:
	return magnet::GL::element_type::TRIANGLES;
      }
  }


  void 
  Glyphs::glRender(const magnet::GL::Camera& cam, RenderMode mode) 
  {    
    _primitiveVertices.getContext()->resetInstanceTransform();

    switch (_glyphType->get_active_row_number())
      {
      case SPHERE_GLYPH:
	{
	  if (_raytraceable && _glyphRaytrace->get_active())
	    {
	      if (_context->testExtension("GL_ARB_sample_shading"))
		{
		  _primitiveVertices.getContext()->setSampleShading(true);
		  glMinSampleShadingARB(1.0);
		}

	      if (mode & RenderObj::SHADOW)
		{
		  _sphereVSMShader.attach();
		  _sphereVSMShader["ProjectionMatrix"] = cam.getProjectionMatrix();
		  _sphereVSMShader["ViewMatrix"] = cam.getViewMatrix();
		  _sphereVSMShader["global_scale"] = _scale;
		  _scaleSel->bindAttribute(magnet::GL::Context::instanceScaleAttrIndex, 0);
		  _colorSel->bindAttribute(magnet::GL::Context::vertexColorAttrIndex, 0);
		  _ds.getPositionBuffer().drawArray(magnet::GL::element_type::POINTS);
		  _sphereVSMShader.detach();

		}
	      else
		{
		  _sphereShader.attach();
		  _sphereShader["ProjectionMatrix"] = cam.getProjectionMatrix();
		  _sphereShader["ViewMatrix"] = cam.getViewMatrix();
		  _sphereShader["global_scale"] = _scale;
		  _scaleSel->bindAttribute(magnet::GL::Context::instanceScaleAttrIndex, 0);
		  _colorSel->bindAttribute(magnet::GL::Context::vertexColorAttrIndex, 0);
		  _ds.getPositionBuffer().drawArray(magnet::GL::element_type::POINTS);
		  _sphereShader.detach();
		}

	      if (_context->testExtension("GL_ARB_sample_shading"))
		_primitiveVertices.getContext()->setSampleShading(false);
	      return;
	    }
	}
	break;
      case ARROW_GLYPH:
      case CYLINDER_GLYPH:
      case LINE_GLYPH:
      case CUBE_GLYPH:
      default:
	break;
      }
    
    if (!_primitiveVertices.size()) return;

    _renderShader.attach();
    _renderShader["ProjectionMatrix"] = cam.getProjectionMatrix();
    _renderShader["ViewMatrix"] = cam.getViewMatrix();

    _ds.getPositionBuffer().attachToAttribute(magnet::GL::Context::instanceOriginAttrIndex, 1);
    _scaleSel->bindAttribute(magnet::GL::Context::instanceScaleAttrIndex, 1);
    _orientSel->bindAttribute(magnet::GL::Context::instanceOrientationAttrIndex, 1);
    _colorSel->bindAttribute(magnet::GL::Context::vertexColorAttrIndex, 1);
    _primitiveVertices.attachToVertex();
    _primitiveNormals.attachToNormal();
    _primitiveIndices.drawInstancedElements(getElementType(), _N);
    _renderShader.detach();
  }

  void 
  Glyphs::deinit()
  {
    _N = 0;
    _primitiveVertices.deinit();
    _primitiveNormals.deinit();
    _primitiveIndices.deinit();
    RenderObj::deinit();
    _sphereShader.deinit();
    _sphereVSMShader.deinit();
    _renderShader.deinit();
    _simpleRenderShader.deinit();
    _gtkOptList.reset();
    _scaleSel.reset(); 
    _colorSel.reset();
    _orientSel.reset();
    _glyphType.reset();
    _glyphLOD.reset();
    _scaleFactorBox.reset();
    _scaleLabel.reset();
    _scaleFactor.reset();
    _glyphRaytrace.reset();
  }

  void 
  Glyphs::showControls(Gtk::ScrolledWindow* win)
  {
    win->remove();
    _gtkOptList->unparent();
    win->add(*_gtkOptList);
    win->show();
  }

  void 
  Glyphs::init(const std::tr1::shared_ptr<magnet::thread::TaskQueue>& systemQueue)
  {
    RenderObj::init(systemQueue);            
    //Initialise the Gtk controls
    _gtkOptList.reset(new Gtk::VBox);
    _gtkOptList->show();

    Gtk::HSeparator* separator;

    //Glyph selection and level of detail
    _glyphBox.reset(new Gtk::HBox); _glyphBox->show();
    
    _context = magnet::GL::Context::getContext();
    _raytraceable = _context->testExtension("GL_EXT_geometry_shader4");

    _renderShader.build();
    _simpleRenderShader.build();

    if (_raytraceable) 
      {
	_sphereShader.build();
	_sphereVSMShader.build();
      }	

    {
      Gtk::Label* label = Gtk::manage(new Gtk::Label("Glyph Type")); label->show();
      _glyphBox->pack_start(*label, false, false, 5);
	
      _glyphType.reset(new Gtk::ComboBoxText);
      _glyphType->show();

      _glyphType->append_text("Sphere");      
      _glyphType->append_text("Arrow");
      _glyphType->append_text("Cylinder");
      _glyphType->append_text("Line");
      _glyphType->append_text("Cube");
      _glyphType->set_active(_initGlyphType);

      _glyphBox->pack_start(*_glyphType, false, false, 5);
    }
    
    {
      _glyphRaytrace.reset(new Gtk::CheckButton("RayTrace")); ;
      if (_raytraceable) 
	_glyphRaytrace->show();
      _glyphRaytrace->set_active(_raytraceable);
      _glyphRaytrace->set_sensitive(_raytraceable);
      _glyphBox->pack_start(*_glyphRaytrace, false, false, 5);
    }

    {
      _glyphLOD.reset(new Gtk::SpinButton(1.0, 0)); _glyphLOD->show();
      _glyphLOD->set_numeric(true);
      _glyphBox->pack_end(*_glyphLOD, false, false, 5);
      Gtk::Label* label = Gtk::manage(new Gtk::Label("Level of Detail")); label->show();
      _glyphBox->pack_end(*label, false, false, 5);
    }

    _gtkOptList->pack_start(*_glyphBox, false, false, 5);

    separator = Gtk::manage(new Gtk::HSeparator);
    separator->show();
    _gtkOptList->pack_start(*separator, false, false, 0);

    //The attribute selectors
    _scaleSel.reset(new AttributeSelector(magnet::GL::Context::instanceScaleAttrIndex));

    _scaleSel->buildEntries("Scale Data Field:", _ds, 1, 4,
			    Attribute::INTENSIVE | Attribute::EXTENSIVE, 3,
			    Attribute::DEFAULT_GLYPH_SCALING);
    _gtkOptList->pack_start(*_scaleSel, false, false);

    _scaleFactorBox.reset(new Gtk::HBox);
    _scaleFactorBox->show();
    _gtkOptList->pack_start(*_scaleFactorBox, false, false, 5);
    _scaleLabel.reset(new Gtk::Label("Scale factor", 1.0, 0.5f));
    _scaleLabel->show();
    _scaleFactorBox->pack_start(*_scaleLabel, true, true, 5);
    _scaleFactor.reset(new Gtk::Entry);
    _scaleFactor->show();
    _scaleFactorBox->pack_start(*_scaleFactor, false, false, 5);
    _scaleFactor->set_text("1.0");
    _scaleFactor->set_width_chars(5);


    separator = Gtk::manage(new Gtk::HSeparator); 
    separator->show(); 
    _gtkOptList->pack_start(*separator, false, false, 0);

    _colorSel.reset(new AttributeColorSelector);
    _colorSel->buildEntries("Color Data Field:", _ds, 1, 4, 
			    Attribute::INTENSIVE | Attribute::EXTENSIVE, 4);
    _gtkOptList->pack_start(*_colorSel, false, false);

    separator = Gtk::manage(new Gtk::HSeparator); 
    separator->show(); 
    _gtkOptList->pack_start(*separator, false, false, 0);

    _orientSel.reset(new AttributeOrientationSelector);
    _orientSel->buildEntries("Orientation Data Field:", _ds, 3, 4, 
			     Attribute::INTENSIVE | Attribute::EXTENSIVE, 4,
			     Attribute::DEFAULT_GLYPH_ORIENTATION);
    _gtkOptList->pack_start(*_orientSel, false, false);

    glyph_type_changed();

    _glyphRaytrace->signal_toggled()
      .connect(sigc::mem_fun(this, &Glyphs::guiUpdate));
    _glyphLOD->signal_value_changed()
      .connect(sigc::mem_fun(*this, &Glyphs::guiUpdate));
    _glyphType->signal_changed()
      .connect(sigc::mem_fun(*this, &Glyphs::glyph_type_changed));
    _scaleFactor->signal_changed()
      .connect(sigc::mem_fun(*this, &Glyphs::guiUpdate));
  }

  void 
  Glyphs::glyph_type_changed()
  {
    _glyphRaytrace->set_sensitive(false);
    _glyphRaytrace->set_active(false);
    _glyphLOD->set_sensitive(true);

    switch (_glyphType->get_active_row_number())
      {
      case SPHERE_GLYPH:
	{	  
	  if (_raytraceable)
	    {
	      _glyphRaytrace->set_sensitive(true);
	      _glyphRaytrace->set_active(true);
	    }
	    
	  _glyphLOD->get_adjustment()->configure(1, 0.0, 4.0, 1.0, 1.0, 0.0);
	}
	break;
      case ARROW_GLYPH:
      case CYLINDER_GLYPH:
	_glyphLOD->get_adjustment()->configure(6.0, 6.0, 32.0, 1.0, 5.0, 0.0);
	break;
      case LINE_GLYPH:
      case CUBE_GLYPH:
	_glyphLOD->set_sensitive(false);
	break;
      default:
	break;
      }

    guiUpdate();
  }


  void 
  Glyphs::guiUpdate()
  {
    magnet::gtk::forceNumericEntry(_scaleFactor.get());
    try {
      _scale = boost::lexical_cast<double>(_scaleFactor->get_text());
    } catch (...) {}

    switch (_glyphType->get_active_row_number())
      {
      case SPHERE_GLYPH:
	_glyphLOD->set_sensitive(true);
	if (_raytraceable)
	  if (_glyphRaytrace->get_active())
	    _glyphLOD->set_sensitive(false);
	break;
      case ARROW_GLYPH:
      case CYLINDER_GLYPH:
      case LINE_GLYPH:
      case CUBE_GLYPH:
      default:
	break;
      }
    
    //No need to deinitialise, we'll just initialise over the top of the old data
    //deinit(); 
    _N = _ds.size();
    
    //Load the primitive data into the VBO's
    _primitiveVertices.init(getPrimitiveVertices(), 3, magnet::GL::buffer_usage::STATIC_DRAW);
    _primitiveNormals.init(getPrimitiveNormals(), 3, magnet::GL::buffer_usage::STATIC_DRAW);
    _primitiveIndices.init(getPrimitiveIndicies(), 3, magnet::GL::buffer_usage::STATIC_DRAW);
  }

  std::vector<GLfloat> 
  Glyphs::getPrimitiveVertices()
  {
    int LOD = _glyphLOD->get_value_as_int();
    int type = _glyphType->get_active_row_number();

    std::vector<GLfloat> vertices;

    switch (type)
      {
      case SPHERE_GLYPH:
	{
	  magnet::GL::objects::primitives::Sphere sph(magnet::GL::objects::primitives::Sphere::icosahedron, LOD);
	  vertices = std::vector<GLfloat>(sph.getVertices(), sph.getVertices() + sph.getVertexCount() * 3);
	  break;
	}
      case ARROW_GLYPH:
	vertices = magnet::GL::objects::primitives::Arrow::getVertices(LOD);
	break;
      case CYLINDER_GLYPH:
	vertices = magnet::GL::objects::primitives::Cylinder::getVertices(LOD);
	break;
      case LINE_GLYPH:
	vertices.clear();
	vertices.push_back(0); vertices.push_back(0); vertices.push_back(0.5);
	vertices.push_back(0); vertices.push_back(0); vertices.push_back(-0.5);
	break;
      case CUBE_GLYPH:
	vertices = magnet::GL::objects::primitives::Cube::getVertices();
	break;
      default:
	M_throw() << "Unrecognised glyph type";
      }

    for (std::vector<GLfloat>::iterator iPtr = vertices.begin();
	 iPtr != vertices.end(); ++iPtr)
      *iPtr *= _scale;
      
    return vertices;
  }
    
  std::vector<GLfloat> 
  Glyphs::getPrimitiveNormals()
  {
    int LOD = _glyphLOD->get_value_as_int();
    int type = _glyphType->get_active_row_number();

    switch (type)
      {
      case SPHERE_GLYPH:
	{
	  magnet::GL::objects::primitives::Sphere sph(magnet::GL::objects::primitives::Sphere::icosahedron, LOD);
	  return std::vector<GLfloat>(sph.getVertices(), sph.getVertices() + sph.getVertexCount() * 3);
	}
      case ARROW_GLYPH:
	return magnet::GL::objects::primitives::Arrow::getNormals(LOD);
      case CYLINDER_GLYPH:
	return magnet::GL::objects::primitives::Cylinder::getNormals(LOD);
      case LINE_GLYPH: 
	{
	  //(Normals are 0's to stop them being shaded)
	  std::vector<GLfloat> normals;
	  normals.push_back(0); normals.push_back(0); normals.push_back(0);
	  normals.push_back(0); normals.push_back(0); normals.push_back(0);
	  return normals;
	}
      case CUBE_GLYPH: //Cube
	return magnet::GL::objects::primitives::Cube::getNormals();
      default:
	M_throw() << "Unrecognised glyph type";
      }
  }

  std::vector<GLuint> 
  Glyphs::getPrimitiveIndicies()
  {
    int LOD = _glyphLOD->get_value_as_int();
    int type = _glyphType->get_active_row_number();

    switch (type)
      {
      case SPHERE_GLYPH:
	{
	  magnet::GL::objects::primitives::Sphere 
	    sph(magnet::GL::objects::primitives::Sphere::icosahedron, LOD);
	  return std::vector<GLuint>(sph.getFaces(), sph.getFaces() + sph.getFaceCount() * 3);
	}
      case ARROW_GLYPH:
	return magnet::GL::objects::primitives::Arrow::getIndices(LOD);
      case CYLINDER_GLYPH:
	return magnet::GL::objects::primitives::Cylinder::getIndices(LOD);
      case LINE_GLYPH:
	{
	  std::vector<GLuint> indices;
	  indices.push_back(0); indices.push_back(1);
	  return indices;
	}
      case CUBE_GLYPH:
	return magnet::GL::objects::primitives::Cube::getIndices();
      default:
	M_throw() << "Unrecognised glyph type";
      }
  }

  void 
  Glyphs::pickingRender(const magnet::GL::Camera& cam, 
			const uint32_t offset)
  {
    _primitiveVertices.getContext()->resetInstanceTransform();

    magnet::GL::Buffer<GLubyte> colorbuf;
    {//Send unique color id's to colorbuf
      std::vector<GLubyte> colors;
      colors.resize(4 * _N);
      for (uint32_t i(0); i < _N; ++i)
	*reinterpret_cast<uint32_t*>(&(colors[4 * i])) = offset + i;
      colorbuf.init(colors, 4);
    }

    switch (_glyphType->get_active_row_number())
      {
      case SPHERE_GLYPH:
	{
	  if (_raytraceable && _glyphRaytrace->get_active())
	    {
	      _sphereShader.attach();
	      _sphereShader["ProjectionMatrix"] = cam.getProjectionMatrix();
	      _sphereShader["ViewMatrix"] = cam.getViewMatrix();
	      _sphereShader["global_scale"] = _scale;
	      _scaleSel->bindAttribute(magnet::GL::Context::instanceScaleAttrIndex, 0);
	      colorbuf.attachToColor();
	      _ds.getPositionBuffer().drawArray(magnet::GL::element_type::POINTS);
	      _sphereShader.detach();
	      return;
	    }
	}
	break;
      case ARROW_GLYPH:
      case CYLINDER_GLYPH:
      case LINE_GLYPH:
      case CUBE_GLYPH:
      default:
	break;
      }
    
    if (!_primitiveVertices.size()) return;

    _simpleRenderShader.attach();
    _simpleRenderShader["ProjectionMatrix"] = cam.getProjectionMatrix();
    _simpleRenderShader["ViewMatrix"] = cam.getViewMatrix();
    _ds.getPositionBuffer().attachToAttribute(magnet::GL::Context::instanceOriginAttrIndex, 1);
    _scaleSel->bindAttribute(magnet::GL::Context::instanceScaleAttrIndex, 1);
    _orientSel->bindAttribute(magnet::GL::Context::instanceOrientationAttrIndex, 1);
    colorbuf.attachToAttribute(magnet::GL::Context::vertexColorAttrIndex, 1, true);
    _primitiveVertices.attachToVertex();
    _primitiveNormals.attachToNormal();
    _primitiveIndices.drawInstancedElements(getElementType(), _N);
    _simpleRenderShader.detach();
  }
}
