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

namespace coil {  
  class Glyphs : public DataSetChild, public magnet::GL::objects::Instanced
  {
    static const size_t _LOD = 6;

  public:
    inline Glyphs(std::string name, DataSet& ds): DataSetChild(name, ds) {}

    inline virtual void clTick(const magnet::GL::Camera& cam) {}

    inline virtual void glRender(magnet::GL::FBO& fbo, const magnet::GL::Camera& cam, RenderMode mode) 
    {
      //Do not allow a glRender if uninitialised
      if (!_primitiveVertices.size()) return;

      _primitiveVertices.getContext().resetInstanceTransform();
      _positionSel->bindAttribute(magnet::GL::Context::instanceOriginAttrIndex);
      _scaleSel->bindAttribute(magnet::GL::Context::instanceScaleAttrIndex);
      Instanced::glRender();
    }

    inline virtual void init(const std::tr1::shared_ptr<magnet::thread::TaskQueue>& systemQueue)
    {
      RenderObj::init(systemQueue);
      Instanced::init(_ds.size());
            
      //Initialise the Gtk controls
      _gtkOptList.reset(new Gtk::VBox);
      _gtkOptList->show();
      
      _positionSel.reset(new AttributeSelector);
      _positionSel->buildEntries("Position Data Field:", _ds, Attribute::COORDINATE, 3);
      _gtkOptList->pack_start(*_positionSel, false, false);

      _scaleSel.reset(new AttributeSelector);
      _scaleSel->buildEntries("Scale Data Field:", _ds, Attribute::COORDINATE, 1);
      _gtkOptList->pack_start(*_scaleSel, false, false);
    }
    
    inline virtual void deinit()
    {
      Instanced::deinit();
      RenderObj::deinit();
      _gtkOptList.reset();
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


    inline virtual magnet::GL::element_type::Enum  getElementType()
    { return magnet::GL::element_type::TRIANGLE_STRIP; }
    
    inline virtual std::vector<GLfloat> getPrimitiveVertices()
    {
      std::vector<GLfloat> vertices(2 * _LOD * 3);
      
      for (size_t vert = 0; vert < 2 * _LOD; ++vert)
	{
	  vertices[3 * vert + 0] = 0.5f * std::sin((vert / 2) * 2.0f * M_PI / _LOD);
	  vertices[3 * vert + 1] = 0.5f * std::cos((vert / 2) * 2.0f * M_PI / _LOD);
	  vertices[3 * vert + 2] = vert % 2;
	}
      return vertices;
    }
    
    inline virtual std::vector<GLfloat> getPrimitiveNormals()
    {
      std::vector<GLfloat> normals(2 * _LOD * 3);
      
      for (size_t vert = 0; vert < 2 * _LOD; ++vert)
	{
	  GLfloat x = 0.5f * std::sin((vert / 2) * 2.0f * M_PI / _LOD);
	  GLfloat y = 0.5f * std::cos((vert / 2) * 2.0f * M_PI / _LOD);
	  GLfloat scale = 1.0f / std::sqrt(x * x + y * y);
	  normals[3 * vert + 0] = x * scale;
	  normals[3 * vert + 1] = y * scale;
	  normals[3 * vert + 2] = 0;
	}
      return normals;
    }
    
    inline virtual std::vector<GLuint>  getPrimitiveIndicies()
    {
      std::vector<GLuint> indices(2 * _LOD + 2);
      
      //Main vertices
      for (size_t vert = 0; vert < 2 * _LOD; ++vert)
	indices[vert] = vert;
      
      //Rejoin the end vertex
      indices[2 * _LOD + 0] =  0;
      indices[2 * _LOD + 1] =  1;
      
      return indices;
    }
    
  protected:
    std::auto_ptr<Gtk::VBox> _gtkOptList;
    std::auto_ptr<AttributeSelector> _positionSel;
    std::auto_ptr<AttributeSelector> _scaleSel;
  };
}
