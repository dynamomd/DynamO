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
#include <coil/RenderObj/RenderObj.hpp>

extern const guint8 default_RObj_Icon[];
extern const size_t default_RObj_Icon_size;

namespace coil {
  Glib::RefPtr<Gdk::Pixbuf> 
  RenderObj::getIcon()
  {
    return Gdk::Pixbuf::create_from_inline(default_RObj_Icon_size, default_RObj_Icon);
  }

  Gtk::TreeModel::iterator 
  RenderObj::addViewRows(RenderObjectsGtkTreeView& view, Gtk::TreeModel::iterator& iter)
  {
    (*iter)[view._columns->m_name] = getName();
    (*iter)[view._columns->m_visible] = visible();
    (*iter)[view._columns->m_shadowcasting] = shadowCasting();
    (*iter)[view._columns->m_obj] = this;
    (*iter)[view._columns->m_icon] = getIcon();
    return iter;
  }

  void
  RenderObjectsGtkTreeView::init(Gtk::TreeView* tree)
  { 
    _columns.reset(new ModelColumns); 
    _view = tree;
    _store = Gtk::TreeStore::create(*_columns);
    _view->set_model(_store);

    {
      Gtk::TreeViewColumn* column
	= Gtk::manage(new Gtk::TreeViewColumn("Object"));
      _view->append_column(*column);

      if (column)
	{
	  column->pack_start(_columns->m_icon, false);
	  column->pack_start(_columns->m_name, true);
	  column->set_sizing(Gtk::TREE_VIEW_COLUMN_AUTOSIZE);
	  column->set_expand(true);
	}
    }

    { //The cell renderer 
      Gtk::CellRendererToggle* renderer =
	Gtk::manage(new Gtk::CellRendererToggle());
      int visible_col = _view->append_column("Visible", *renderer);
      Gtk::TreeViewColumn* column = _view->get_column(visible_col -1);
      if (column)
	{
	  column->add_attribute(renderer->property_active(), _columns->m_visible);
	  column->set_sizing(Gtk::TREE_VIEW_COLUMN_AUTOSIZE);
	  column->set_expand(false);
	}
      renderer->signal_toggled().connect(sigc::mem_fun(*this, &RenderObjectsGtkTreeView::visibleToggled));
    }

    { //The cell renderer 
      Gtk::CellRendererToggle* renderer =
	Gtk::manage(new Gtk::CellRendererToggle());
      int shadow_col = _view->append_column("Shadow Casting", *renderer);
      Gtk::TreeViewColumn* column = _view->get_column(shadow_col -1);
      if (column)
	{
	  column->add_attribute(renderer->property_active(), _columns->m_shadowcasting);
	  column->set_sizing(Gtk::TREE_VIEW_COLUMN_AUTOSIZE);
	  column->set_expand(false);
	}
      renderer->signal_toggled().connect(sigc::mem_fun(*this, &RenderObjectsGtkTreeView::shadowingToggled));
    }

    _view->set_enable_tree_lines(true);
    //_view->set_level_indentation(5);
  }

  void 
  RenderObjectsGtkTreeView::buildRenderView()
  {
    _store->clear();
    
    for (std::vector<std::tr1::shared_ptr<RenderObj> >::iterator iPtr = _renderObjects.begin();
	 iPtr != _renderObjects.end(); ++iPtr)
      {
	Gtk::TreeModel::iterator iter = _store->append();
	(*iPtr)->addViewRows(*this, iter);
      }
  }

  void 
  RenderObjectsGtkTreeView::visibleToggled(const Glib::ustring& path_string)
  {
    Gtk::TreeModel::iterator iter 
      = _store->get_iter(path_string);
    
    bool visible = !(*iter)[_columns->m_visible];
    (*iter)[_columns->m_visible] = visible;
    
    RenderObj* obj = (*iter)[_columns->m_obj];
    obj->setVisible(visible);
  }

  void 
  RenderObjectsGtkTreeView::shadowingToggled(const Glib::ustring& path_string)
  {
    Gtk::TreeModel::iterator iter 
      = _store->get_iter(path_string);

    bool casting = !(*iter)[_columns->m_shadowcasting];
    (*iter)[_columns->m_shadowcasting] = casting;

    RenderObj* obj = (*iter)[_columns->m_obj];
    obj->setShadowCasting(casting);
  }
}
