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

extern const guint8 visible_on_Icon[];
extern const size_t visible_on_Icon_size;
extern const guint8 visible_off_Icon[];
extern const size_t visible_off_Icon_size;

extern const guint8 shadow_on_Icon[];
extern const size_t shadow_on_Icon_size;
extern const guint8 shadow_off_Icon[];
extern const size_t shadow_off_Icon_size;

extern const guint8 delete_Icon[];
extern const size_t delete_Icon_size;

extern const guint8 delete_off_Icon[];
extern const size_t delete_off_Icon_size;

namespace {
  Glib::RefPtr<Gdk::Pixbuf> 
  getVisibleIcon(bool enabled)
  {
    if (enabled) 
      return Gdk::Pixbuf::create_from_inline(visible_on_Icon_size, visible_on_Icon);
    else
      return Gdk::Pixbuf::create_from_inline(visible_off_Icon_size, visible_off_Icon);
  }

  Glib::RefPtr<Gdk::Pixbuf> 
  getShadowIcon(bool enabled)
  {
    if (enabled) 
      return Gdk::Pixbuf::create_from_inline(shadow_on_Icon_size, shadow_on_Icon);
    else
      return Gdk::Pixbuf::create_from_inline(shadow_off_Icon_size, shadow_off_Icon);
  }

  Glib::RefPtr<Gdk::Pixbuf> 
  getDeleteIcon(bool enabled)
  {
    if (enabled) 
      return Gdk::Pixbuf::create_from_inline(delete_Icon_size, delete_Icon);
    else
      return Gdk::Pixbuf::create_from_inline(delete_off_Icon_size, delete_off_Icon);
  }
}

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
    (*iter)[view._columns->m_obj] = this;
    (*iter)[view._columns->m_icon] = getIcon();
    (*iter)[view._columns->m_visible] = getVisibleIcon(visible());
    (*iter)[view._columns->m_shadowcasting] = getShadowIcon(shadowCasting());
    (*iter)[view._columns->m_delete] = getDeleteIcon(deletable());
    return iter;
  }

  void
  RenderObjectsGtkTreeView::init(Gtk::TreeView* tree)
  { 
    _columns.reset(new ModelColumns); 
    _view = tree;
    _store = Gtk::TreeStore::create(*_columns);
    _view->set_model(_store);

    _view->append_column("Icon", _columns->m_icon);
      
    {
      Gtk::CellRendererText* renderer =
	Gtk::manage(new Gtk::CellRendererText());
      int cols_count = _view->append_column("Name", *renderer);
      Gtk::TreeViewColumn* column = _view->get_column(cols_count - 1);
      
      if (column)
	{
	  column->add_attribute(renderer->property_text(), _columns->m_name);
	  renderer->property_editable() = true;
	  renderer->signal_edited()
	    .connect(sigc::mem_fun(*this, &RenderObjectsGtkTreeView::name_edited));
	  column->set_sizing(Gtk::TREE_VIEW_COLUMN_AUTOSIZE);
	  column->set_expand(true);
	}
    }

    _view->append_column("Visible", _columns->m_visible);
    _view->append_column("Shadow", _columns->m_shadowcasting);
    _view->append_column("Delete", _columns->m_delete);
    
    _view->set_headers_visible(false);
    _view->set_enable_tree_lines(true);
    _view->add_events(::Gdk::BUTTON_PRESS_MASK);
    _view->signal_button_press_event()
      .connect(sigc::mem_fun(*this, &RenderObjectsGtkTreeView::button_press), false);
  }

  void
  RenderObjectsGtkTreeView::name_edited(const Glib::ustring& path, const Glib::ustring& string)
  {
    Gtk::TreeModel::iterator iter = _store->get_iter(path);
    if (iter)
      {
	RenderObj* obj = (*iter)[_columns->m_obj];
	obj->setName(string);
	(*iter)[_columns->m_name] = string;
      }
  }
  

  bool 
  RenderObjectsGtkTreeView::button_press(GdkEventButton* btn)
  {
    if (btn)
      {
	if ((btn->button == 1) && (btn->type == GDK_BUTTON_PRESS))
	  {
	    Gtk::TreeModel::Path path;
	    Gtk::TreeViewColumn* clm;
	    int cell_x, celly;
	    _view->get_path_at_pos(btn->x, btn->y, path, clm, cell_x, celly);
	    Gtk::TreeModel::iterator iter = _store->get_iter(path);
	    if (clm && iter)
	      {
		if (clm->get_title() == "Visible")
		  {		    
		    RenderObj* obj = (*iter)[_columns->m_obj];
		    bool visible = !obj->visible();
		    obj->setVisible(visible);
		    (*iter)[_columns->m_visible] = getVisibleIcon(visible);
		    
		    obj->setVisible(visible);
		  }
		else if (clm->get_title() == "Shadow")
		  {
		    RenderObj* obj = (*iter)[_columns->m_obj];
		    bool casting = !obj->shadowCasting();
		    obj->setShadowCasting(casting);
		    (*iter)[_columns->m_shadowcasting] = getShadowIcon(casting);		    
		  }
		else if (clm->get_title() == "Delete")
		  {
		    RenderObj* obj = (*iter)[_columns->m_obj];
		    obj->request_delete();
		  }
	      }
	  }
      }

    return false;
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
}
