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
#include <iostream>

namespace coil {  
  void 
  DataSet::init(const std::tr1::shared_ptr<magnet::thread::TaskQueue>& systemQueue)
  {
    _context = &(magnet::GL::Context::getContext());
    RenderObj::init(systemQueue); 
    initGtk(); 
    
    for (std::vector<DataSetChild>::iterator iPtr = _children.begin();
	 iPtr != _children.end(); ++iPtr)
      iPtr->init(systemQueue);
    
    //We don't initialise the attributes, as they're initialised on access
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

      _attrview.reset(new Gtk::TreeView);
      _attrview->set_model(_attrtreestore);
      _attrview->append_column("Name", _attrcolumns->name);
      _attrview->append_column("Components", _attrcolumns->components);
      _attrview->show();
      Gtk::ScrolledWindow* win = Gtk::manage(new Gtk::ScrolledWindow);
      win->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);
      win->add(*_attrview);
      _gtkOptList->pack_start(*win, true, true);
      win->show();
    }

    {
      Gtk::Button* btn = Gtk::manage(new Gtk::Button("Add Glyph"));
      btn->signal_clicked().connect(sigc::mem_fun(*this, &DataSet::addGlyphs));
      btn->show();
      _gtkOptList->pack_start(*btn, false, false);
    }

    _gtkOptList->show();
    rebuildGui();
  }


  void 
  DataSet::addGlyphs()
  {
    std::cerr << "Adding glyphs!";
  }

  void 
  DataSet::addAttribute(std::string name, Attribute::AttributeType type, size_t components)
  {
    if (find(name) != end())
      M_throw() << "Trying to add an Attribute with a existing name, " << name;
      
    insert(value_type(name, Attribute(_N, type, components)));

    //If we're initialised, we should rebuild the view of attributes
    if (_context)
      _context->queueTask(magnet::function::Task::makeTask(&DataSet::rebuildGui, this));
  }
  
  void
  DataSet::rebuildGui()
  {
    for (iterator iPtr = begin(); iPtr != end(); ++iPtr)
      {
	Gtk::TreeModel::iterator iter = _attrtreestore->append();
	(*iter)[_attrcolumns->name] = iPtr->first;
	(*iter)[_attrcolumns->components] = iPtr->second.getNComponents();
      }
  }

  void
  DataSet::deinit()
  {
    _gtkOptList.reset();
    _attrcolumns.reset();
    _attrview.reset();
    _attrtreestore.reset();
    for (std::vector<DataSetChild>::iterator iPtr = _children.begin();
	 iPtr != _children.end(); ++iPtr)
      iPtr->deinit();

    for (iterator iPtr = begin(); iPtr != end(); ++iPtr)
      iPtr->second.deinit();

    _context = NULL;
    RenderObj::deinit();
  }
}
