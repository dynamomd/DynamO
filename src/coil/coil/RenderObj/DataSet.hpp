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

#pragma once
#include <coil/RenderObj/RenderObj.hpp>
#include <coil/RenderObj/Attribute.hpp>
#include <magnet/GL/buffer.hpp>
#include <magnet/gtk/numericEntry.hpp>
#include <magnet/gtk/colorMapSelector.hpp>
#include <boost/lexical_cast.hpp>
#include <vector>
#include <memory>

namespace coil {
  class DataSet;
  class AttributeSelector;

  class DataSetChild: public RenderObj
  {
  public:
    inline DataSetChild(std::string name, DataSet& ds): RenderObj(name), _ds(ds) {}
    
    virtual bool deletable() { return true; }

    /*! \brief Called when the object should be deleted. */
    virtual void request_delete();

    virtual std::array<GLfloat, 4> getCursorPosition(uint32_t objID);

    virtual std::string getCursorText(uint32_t objID);

  protected:

    DataSet& _ds;
  };
    
  /*! \brief A container class for a collection of \ref Attribute
   * instances forming a dataset, and any active filters/glyphs or any
   * other type derived from \ref DataSetChild.
   */
  class DataSet: public RenderObj
  {
    std::map<std::string, std::shared_ptr<Attribute>> _attributes;
    std::map<std::string, magnet::GL::Buffer<GLuint>> _pointSets;
    std::map<std::string, magnet::GL::Buffer<GLuint>> _linkSets;
    
  public:
    DataSet(std::string name, size_t N, int defaultGlyphType = 0):
      RenderObj(name), 
      _N(N),
      _defaultGlyphType(defaultGlyphType)
    {}
    
    virtual void init(const std::shared_ptr<magnet::thread::TaskQueue>& systemQueue);

    virtual void deinit();

    std::map<std::string, std::shared_ptr<Attribute> >& getAttributes() { return _attributes; }
    const std::map<std::string, std::shared_ptr<Attribute> >& getAttributes() const { return _attributes; }

    void addPoints(std::string name, const std::vector<GLuint>&);

    /** @name The host code interface. */
    /**@{*/

    /*! \brief Add an Attribute to the DataSet.
     
      \param name The name of the Attribute.
      \param type The type of the attribute.
      \param components The number of components per value of the attribute.
     */
    void addAttribute(std::string name, int type, size_t components);

    /*! \brief Looks up an attribute by its name.
     */
    inline Attribute& operator[](const std::string& name)
    {
      auto iPtr = _attributes.find(name);
      if (iPtr == _attributes.end())
	M_throw() << "No attribute named " << name << " in Data set";
      
      return *(iPtr->second);
    }

    inline size_t size() const { return _N; }

    void setPeriodicVectors(Vector x, Vector y, Vector z)
    {
      _periodicImageX = x; 
      _periodicImageY = y; 
      _periodicImageZ = z;
    }

    /**@}*/
        
    virtual void glRender(const magnet::GL::Camera& cam, RenderMode mode, const uint32_t offset)
    {
      if (mode == PICKING)
	{
	  uint32_t obj_offset = offset;
	  
	  for (auto& child : _children)
	    {
	      uint32_t objs = child->pickableObjectCount();
	      if (objs) child->glRender(cam, RenderObj::PICKING, obj_offset);
	      obj_offset += objs;
	    }
	}
      else
	{
	  for (auto& child : _children)
	    if (child->visible() && (!(mode & SHADOW) || child->shadowCasting()))
	      child->glRender(cam, mode);
	  
	  for (auto& data : _attributes)
	    data.second->renderComplete();
	  
	  rebuildGui();
	}
    }
    
    virtual Gtk::TreeModel::iterator addViewRows(RenderObjectsGtkTreeView& view, Gtk::TreeModel::iterator& iter)
    {
      _view = &view;
      
      _iter = iter;
      RenderObj::addViewRows(view, _iter);

      for (auto& child : _children)
	{
	  Gtk::TreeModel::iterator child_iter = view._store->append(_iter->children());
	  child->addViewRows(view, child_iter);
	}
      return _iter;
    }

    virtual void showControls(Gtk::ScrolledWindow* win);

    virtual Glib::RefPtr<Gdk::Pixbuf> getIcon();

    void deleteChild(DataSetChild* child)
    {
      _context->queueTask(std::bind(&DataSet::deleteChildWorker, this, child));
    }

    magnet::GL::Context::ContextPtr getContext() { return _context; }

    virtual uint32_t pickableObjectCount()
    {
      if (!visible()) return 0;

      uint32_t sum = 0;
      for (auto& child : _children) sum += child->pickableObjectCount();

      return sum;
    }

    virtual std::shared_ptr<RenderObj>
    getPickedObject(uint32_t& objID, const std::shared_ptr<RenderObj>& my_ptr)
    { 
      uint32_t obj_offset = 0;

      for (auto& child : _children)
	{
	  const uint32_t n_objects = child->pickableObjectCount();
	  
	  if ((objID >= obj_offset) && (objID - obj_offset) < n_objects)
	    {
	      objID -= obj_offset;
	      return child->getPickedObject(objID, child);
	    }

	  obj_offset += n_objects;
	}

      M_throw() << "The selected object was not drawn by this RenderObj";
    }

    virtual magnet::math::Vector getMinCoord() const;
    virtual magnet::math::Vector getMaxCoord() const;

    magnet::GL::Buffer<GLfloat>& getPositionBuffer();

    std::shared_ptr<AttributeSelector>& getPositionSelector() { return _positionSel; }

    void addGlyphs();
        
    virtual std::string getCursorText(uint32_t objID);

    virtual std::array<GLfloat, 4> getCursorPosition(uint32_t objID);

    Vector getPeriodicVectorX() const { return _periodicImageX; }
    Vector getPeriodicVectorY() const { return _periodicImageY; }
    Vector getPeriodicVectorZ() const { return _periodicImageZ; }

  protected:
    void deleteChildWorker(DataSetChild* child);

    void addPointsWorker(std::string name, std::vector<GLuint>);

    /*! \brief An iterator to this DataSet's row in the Render object
      treeview.
     */
    Gtk::TreeModel::iterator _iter;
    RenderObjectsGtkTreeView* _view;
    magnet::GL::Context::ContextPtr _context;
    std::unique_ptr<Gtk::VBox> _gtkOptList;
    size_t _N;
    std::vector<std::shared_ptr<DataSetChild> > _children;
    std::shared_ptr<AttributeSelector> _positionSel;

    int _defaultGlyphType;

    void initGtk();
    void rebuildGui();

    struct ModelColumns : Gtk::TreeModelColumnRecord
    {
      ModelColumns()
      { add(name); add(components); add(type); add(min); add(max); }
      
      Gtk::TreeModelColumn<Glib::ustring> name;
      Gtk::TreeModelColumn<size_t> components;
      Gtk::TreeModelColumn<std::string> min;
      Gtk::TreeModelColumn<std::string> max;
      Gtk::TreeModelColumn<Attribute::AttributeType> type;
    };
    
    std::unique_ptr<ModelColumns> _attrcolumns;
    Glib::RefPtr<Gtk::TreeStore> _attrtreestore;
    std::unique_ptr<Gtk::TreeView> _attrview;
    std::unique_ptr<Gtk::Label> _infolabel;
    std::unique_ptr<Gtk::ComboBoxText> _comboPointSet;

    Vector _periodicImageX;
    Vector _periodicImageY;
    Vector _periodicImageZ;
  };
}
