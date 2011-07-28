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
#include <gtkmm.h>
#include <magnet/GL/context.hpp>
#include <magnet/thread/taskQueue.hpp>
#include <magnet/GL/camera.hpp>
#include <magnet/GL/FBO.hpp>
#include <tr1/memory>

namespace Gtk { class ScrolledWindow; }

namespace coil {
  class RenderObjectsGtkTreeView;

  class RenderObj
  {
  public:
    RenderObj(std::string name): _name(name), _visible(true) {}
  
    virtual void init(const std::tr1::shared_ptr<magnet::thread::TaskQueue>& systemQueue) 
    { _systemQueue = systemQueue; }
    virtual void deinit() {}
    virtual void clTick(const magnet::GL::Camera& cam) = 0;
    virtual void glRender(magnet::GL::FBO& fbo, const magnet::GL::Camera& cam) = 0;
    virtual void interfaceRender(const magnet::GL::Camera& camera) {}
    virtual void initPicking(cl_uint& offset) {}
    virtual void pickingRender(magnet::GL::FBO& fbo, const magnet::GL::Camera& cam) {}
    virtual void finishPicking(cl_uint& offset, const cl_uint val) {}
    virtual void showControls(Gtk::ScrolledWindow* win) {}
    virtual inline void addViewRows(RenderObjectsGtkTreeView& view);
    inline void setVisible(bool val) { _visible = val; }
    inline bool isVisible() const { return _visible; }
    inline const std::string& getName() const { return _name; }
    std::tr1::shared_ptr<magnet::thread::TaskQueue> getQueue() { return _systemQueue; }

  protected:
    std::string _name;
    bool _visible;
    std::tr1::shared_ptr<magnet::thread::TaskQueue> _systemQueue;
  };

  class RenderObjectsGtkTreeView
  {
  public:
    void init(Gtk::TreeView* tree)
    { 
      _columns.reset(new ModelColumns); 
      _view = tree;
      _store = Gtk::TreeStore::create(*_columns);
      _view->set_model(_store);

      { //The cell renderer 
	Gtk::CellRendererToggle* renderer =
	  Gtk::manage(new Gtk::CellRendererToggle());
	int visible_col = _view->append_column("Visible", *renderer);
	Gtk::TreeViewColumn* column = _view->get_column(visible_col -1);
	if (column)
	  column->add_attribute(renderer->property_active(), _columns->m_visible);
	renderer->signal_toggled().connect(sigc::mem_fun(*this, &RenderObjectsGtkTreeView::visibleToggled));
      }

      _view->append_column("Object Name", _columns->m_name);

    }
    
    struct ModelColumns : Gtk::TreeModelColumnRecord
    {
      ModelColumns()
      { add(m_name); add(m_visible); add(m_obj);}
      
      Gtk::TreeModelColumn<Glib::ustring> m_name;
      Gtk::TreeModelColumn<bool> m_visible;
      Gtk::TreeModelColumn<RenderObj*> m_obj;
    };
    
    std::auto_ptr<ModelColumns> _columns;
    Glib::RefPtr<Gtk::TreeStore> _store;
    Gtk::TreeView* _view;
    
  protected:
    
    void visibleToggled(const Glib::ustring& path_string)
    {
      Gtk::TreeModel::iterator iter 
	= _store->get_iter(path_string);

      bool visible = !(*iter)[_columns->m_visible];
      (*iter)[_columns->m_visible] = visible;

      RenderObj* obj = (*iter)[_columns->m_obj];
      obj->setVisible(visible);
    }
  };

  void RenderObj::addViewRows(RenderObjectsGtkTreeView& view)
  {
    Gtk::TreeModel::iterator iter = view._store->append();
    
    (*iter)[view._columns->m_name] = getName();
    (*iter)[view._columns->m_visible] = isVisible();
    (*iter)[view._columns->m_obj] = this;
  }
}
