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
#include <coil/RenderObj/RenderObj.hpp>
#include <coil/RenderObj/Attribute.hpp>
#include <magnet/GL/buffer.hpp>
#include <magnet/gtk/numericEntry.hpp>
#include <magnet/gtk/colorMapSelector.hpp>
#include <boost/lexical_cast.hpp>
#include <vector>
#include <tr1/memory>

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

  protected:

    DataSet& _ds;
  };
    
  /*! \brief A container class for a collection of \ref Attribute
   * instances forming a dataset, and any active filters/glyphs or any
   * other type derived from \ref DataSetChild.
   */
  class DataSet: public RenderObj, public std::map<std::string, std::tr1::shared_ptr<Attribute> >
  {
  public:
    DataSet(std::string name, size_t N): 
      RenderObj(name), 
      _N(N)
    {}
    
    virtual void init(const std::tr1::shared_ptr<magnet::thread::TaskQueue>& systemQueue);

    virtual void deinit();

    /** @name The host code interface. */
    /**@{*/

    /*! \brief Add an Attribute to the DataSet.
     *
     * \param name The name of the Attribute.
     * \param type The type of the attribute.
     * \param components The number of components per value of the attribute.
     */
    void addAttribute(std::string name, int type, size_t components);

    /*! \brief Looks up an attribute by its name.
     */
    inline Attribute& operator[](const std::string& name)
    {
      iterator iPtr = find(name);
      if (iPtr == end())
	M_throw() << "No attribute named " << name << " in Data set";
      
      return *(iPtr->second);
    }

    inline size_t size() const { return _N; }

    /**@}*/
        
    virtual void clTick(const magnet::GL::Camera& cam)
    {
      for (std::vector<std::tr1::shared_ptr<DataSetChild> >::iterator iPtr = _children.begin();
	   iPtr != _children.end(); ++iPtr)
	(*iPtr)->clTick(cam);
    }

    virtual void glRender(magnet::GL::FBO& fbo, const magnet::GL::Camera& cam, RenderMode mode = DEFAULT)
    {
      for (std::vector<std::tr1::shared_ptr<DataSetChild> >::iterator iPtr = _children.begin();
	   iPtr != _children.end(); ++iPtr)
	if ((*iPtr)->visible() && (!(mode & SHADOW) || (*iPtr)->shadowCasting()))
	  (*iPtr)->glRender(fbo, cam, mode);
      
      for (iterator iPtr = begin(); iPtr != end(); ++iPtr)
	iPtr->second->renderComplete();
    }
    
    virtual Gtk::TreeModel::iterator addViewRows(RenderObjectsGtkTreeView& view, Gtk::TreeModel::iterator& iter)
    {
      _view = &view;
      
      _iter = iter;
      RenderObj::addViewRows(view, _iter);

      for (std::vector<std::tr1::shared_ptr<DataSetChild> >::iterator iPtr = _children.begin();
	   iPtr != _children.end(); ++iPtr)
	{
	  Gtk::TreeModel::iterator child_iter = view._store->append(_iter->children());
	  (*iPtr)->addViewRows(view, child_iter);
	}
      return _iter;
    }

    virtual void showControls(Gtk::ScrolledWindow* win);

    virtual Glib::RefPtr<Gdk::Pixbuf> getIcon();

    void deleteChild(DataSetChild* child)
    {
      _context->queueTask(magnet::function::Task::makeTask(&DataSet::deleteChildWorker, this, child));
    }

    magnet::GL::Context::ContextPtr getContext()
    { return _context; }
    

    virtual uint32_t pickableObjectCount()
    {
      if (!visible()) return 0;

      uint32_t sum = 0;
      for (std::vector<std::tr1::shared_ptr<DataSetChild> >::iterator iPtr = _children.begin();
	   iPtr != _children.end(); ++iPtr)
	sum += (*iPtr)->pickableObjectCount();

      return sum;
    }

    virtual void pickingRender(magnet::GL::FBO& fbo, 
			       const magnet::GL::Camera& cam, 
			       const uint32_t offset)
    {
      uint32_t obj_offset = offset;
      for (std::vector<std::tr1::shared_ptr<DataSetChild> >::iterator 
	     iPtr = _children.begin(); iPtr != _children.end(); ++iPtr)
	{
	  uint32_t objs = (*iPtr)->pickableObjectCount();
	  if (objs)
	    (*iPtr)->pickingRender(fbo, cam, obj_offset);

	  obj_offset += objs;
	}
    }

    RenderObj* getPickedObject(uint32_t objID)
    { 
      uint32_t obj_offset = 0;

      for (std::vector<std::tr1::shared_ptr<DataSetChild> >::iterator 
	     iPtr = _children.begin(); iPtr != _children.end(); ++iPtr)
	{
	  const uint32_t n_objects = (*iPtr)->pickableObjectCount();
	  
	  if ((objID >= obj_offset) && (objID - obj_offset) < n_objects)
	    return (*iPtr)->getPickedObject(objID - obj_offset);

	  obj_offset += n_objects;
	}

      return NULL;
    }

    magnet::GL::Buffer<GLfloat>& getPositionBuffer();

    void addGlyphs();
        
    virtual std::string getCursorText(uint32_t objID);

    virtual std::tr1::array<GLfloat, 4> getCursorPosition(uint32_t objID);

  protected:
    void deleteChildWorker(DataSetChild* child);

    /*! \brief An iterator to this DataSet's row in the Render object
      treeview.
     */
    Gtk::TreeModel::iterator _iter;
    RenderObjectsGtkTreeView* _view;
    magnet::GL::Context::ContextPtr _context;
    std::auto_ptr<Gtk::VBox> _gtkOptList;
    size_t _N;
    std::vector<std::tr1::shared_ptr<DataSetChild> > _children;
    std::tr1::shared_ptr<AttributeSelector> _positionSel;

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
    
    std::auto_ptr<ModelColumns> _attrcolumns;
    Glib::RefPtr<Gtk::TreeStore> _attrtreestore;
    std::auto_ptr<Gtk::TreeView> _attrview;
  };
}
