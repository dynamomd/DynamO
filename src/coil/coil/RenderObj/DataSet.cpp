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
#include <coil/RenderObj/Glyphs.hpp>

namespace coil {  
  void 
  DataSet::init(const std::tr1::shared_ptr<magnet::thread::TaskQueue>& systemQueue)
  {
    RenderObj::init(systemQueue); 
    initGtk(); 
    
    for (std::vector<std::tr1::shared_ptr<DataSetChild> >::iterator iPtr = _children.begin();
	 iPtr != _children.end(); ++iPtr)
      (*iPtr)->init(systemQueue);
    
    //We don't initialise the attributes, as they're initialised on access
    _context = &(magnet::GL::Context::getContext());
  }


  void 
  DataSet::showControls(Gtk::ScrolledWindow* win)
  {
    win->remove();
    _gtkOptList->unparent();
    win->add(*_gtkOptList);
    win->show();
  }

  void 
  DataSet::initGtk()
  {
    _gtkOptList.reset(new Gtk::VBox);

    { _attrcolumns.reset(new ModelColumns);
      _attrtreestore = Gtk::TreeStore::create(*_attrcolumns);
      _attrtreestore->set_sort_column(_attrcolumns->components, Gtk::SORT_DESCENDING);

      _attrview.reset(new Gtk::TreeView);
      _attrview->set_model(_attrtreestore);
      _attrview->append_column("Name", _attrcolumns->name);
      _attrview->append_column("Components", _attrcolumns->components);
      _attrview->show();
      Gtk::ScrolledWindow* win = Gtk::manage(new Gtk::ScrolledWindow);
      win->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);
      win->add(*_attrview);
      win->set_size_request(-1, 150);

      Gtk::Frame* frame = Gtk::manage(new Gtk::Frame("Available Attributes")); frame->show();
      frame->add(*win);
      _gtkOptList->pack_start(*frame, false, false, 5);
      win->show();
    }

    {
      Gtk::Button* btn = Gtk::manage(new Gtk::Button("Add Glyph"));
      btn->signal_clicked().connect(sigc::mem_fun(*this, &DataSet::addGlyphs));
      btn->show();
      _gtkOptList->pack_start(*btn, false, false, 5);
    }

    _gtkOptList->show();
    rebuildGui();
  }

  void 
  DataSet::addGlyphs()
  {
    std::tr1::shared_ptr<Glyphs> glyph(new Glyphs("Glyphs", *this));
    _children.push_back(glyph);
        
    if (_context)
      _children.back()->init(_systemQueue);

    if (_iter)
      {
	Gtk::TreeModel::iterator child_iter = _children.back()->addViewRows(*_view, _iter);
	_view->_view->expand_to_path(_view->_store->get_path(child_iter));
      }
  }

  void 
  DataSet::addAttribute(std::string name, Attribute::AttributeType type, size_t components)
  {
    if (find(name) != end())
      M_throw() << "Trying to add an Attribute with a existing name, " << name;

    //Spinlock to force that the Data set is initialised before the attribute is created
    for (;;) if (_context) break;
    
    mapped_type ptr(new Attribute(_N, type, components, _context));
    insert(value_type(name, ptr));

    //If we're initialised, we should rebuild the view of attributes
    if (_context)
      _context->queueTask(magnet::function::Task::makeTask(&DataSet::rebuildGui, this));
  }
  
  void
  DataSet::rebuildGui()
  {
    _attrtreestore->clear();
    for (iterator iPtr = begin(); iPtr != end(); ++iPtr)
      {
	Gtk::TreeModel::iterator iter = _attrtreestore->append();
	(*iter)[_attrcolumns->name] = iPtr->first;
	(*iter)[_attrcolumns->components] = iPtr->second->components();
      }
  }

  void
  DataSet::deinit()
  {
    _gtkOptList.reset();
    _attrcolumns.reset();
    _attrview.reset();
    _attrtreestore.reset();
    for (std::vector<std::tr1::shared_ptr<DataSetChild> >::iterator iPtr = _children.begin();
	 iPtr != _children.end(); ++iPtr)
      (*iPtr)->deinit();

    for (iterator iPtr = begin(); iPtr != end(); ++iPtr)
      iPtr->second->deinit();

    _context = NULL;
    RenderObj::deinit();
  }
}
