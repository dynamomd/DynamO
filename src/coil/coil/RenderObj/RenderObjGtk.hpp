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
#include <memory>

namespace coil {
  class RenderObj;

  struct RenderObjectsGtkTreeView
  {
    void init(Gtk::TreeView* tree)
    { 
      _columns.reset(new ModelColumns); 
      _view = tree;
      _store = Gtk::TreeStore::create(*_columns);
      _view->set_model(_store);
      _view->append_column("Visible", _columns->m_visible);
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
  };
}
