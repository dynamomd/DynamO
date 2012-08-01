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

#include <coil/RenderObj/DataSet.hpp>
#include <coil/RenderObj/Glyphs.hpp>

extern const guint8 DataSet_Icon[];
extern const size_t DataSet_Icon_size;

namespace coil {  
  Glib::RefPtr<Gdk::Pixbuf> 
  DataSet::getIcon()
  {
    return Gdk::Pixbuf::create_from_inline(DataSet_Icon_size, DataSet_Icon);
  }

  void 
  DataSetChild::request_delete()
  { _ds.deleteChild(this); }

  void 
  DataSet::deleteChildWorker(DataSetChild* child)
  {
    for (std::vector<std::tr1::shared_ptr<DataSetChild> >::iterator iPtr = _children.begin();
	 iPtr != _children.end(); ++iPtr)
      if (iPtr->get() == child)
	{
	  //Found the child to delete
	  (*iPtr)->deinit();
	  _children.erase(iPtr);

	  //Rebuild this objects gui
	  rebuildGui();	  
	  //And the render object tree view
	  _view->buildRenderView();
	  return;
	}
  }

  std::tr1::array<GLfloat, 4> 
  DataSetChild::getCursorPosition(uint32_t objID)
  {
    return _ds.getCursorPosition(objID); 
  }

  std::string 
  DataSetChild::getCursorText(uint32_t objID)
  { return _ds.getCursorText(objID); }

  void 
  DataSet::init(const std::tr1::shared_ptr<magnet::thread::TaskQueue>& systemQueue)
  {
    RenderObj::init(systemQueue); 
    initGtk(); 
    
    for (std::vector<std::tr1::shared_ptr<DataSetChild> >::iterator iPtr = _children.begin();
	 iPtr != _children.end(); ++iPtr)
      (*iPtr)->init(systemQueue);
    
    //We don't initialise the attributes, as they're initialised on access
    _context = magnet::GL::Context::getContext();
  }


  void 
  DataSet::showControls(Gtk::ScrolledWindow* win)
  {
    win->remove();
    _gtkOptList->unparent();
    win->add(*_gtkOptList);
    win->show();
    //Force a rebuild of the Gui when this object is selected, to
    //allow stuff like the data set statistics to be updated
    rebuildGui();
  }

  void 
  DataSet::initGtk()
  {
    _gtkOptList.reset(new Gtk::VBox);

    {//The heading of the data set window
      Gtk::HBox* box = Gtk::manage(new Gtk::HBox); box->show();
      Gtk::Frame* frame = Gtk::manage(new Gtk::Frame("Data Set Statistics")); frame->show();
      box->pack_start(*frame, true, true, 5);
      Gtk::VBox* vbox = Gtk::manage(new Gtk::VBox); vbox->show();
      frame->add(*vbox);

      {
	Gtk::Label* label 
	  = Gtk::manage(new Gtk::Label("Points: " 
				       + boost::lexical_cast<std::string>(_N))); 
	label->show();
	vbox->pack_start(*label, false, false, 5);	
      }

      {
	Gtk::Button* btn = Gtk::manage(new Gtk::Button("Add Glyph"));
	btn->signal_clicked().connect(sigc::mem_fun(*this, &DataSet::addGlyphs));
	btn->show();
	box->pack_start(*btn, false, false, 5);
      }
      
      _gtkOptList->pack_start(*box, false, false, 5);
    }
    
    { 
      _attrcolumns.reset(new ModelColumns);
      _attrtreestore = Gtk::TreeStore::create(*_attrcolumns);
      _attrtreestore->set_sort_column(_attrcolumns->components, Gtk::SORT_DESCENDING);

      _attrview.reset(new Gtk::TreeView);
      _attrview->set_model(_attrtreestore);
      _attrview->append_column("Name", _attrcolumns->name);
      _attrview->append_column("Components", _attrcolumns->components);
      _attrview->append_column("Min Values", _attrcolumns->min);
      _attrview->append_column("Max Values", _attrcolumns->max);
      _attrview->show();
      Gtk::ScrolledWindow* win = Gtk::manage(new Gtk::ScrolledWindow);
      win->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_AUTOMATIC);
      win->add(*_attrview);

      Gtk::Frame* frame = Gtk::manage(new Gtk::Frame("Available Attributes")); frame->show();
      frame->add(*win);
      _gtkOptList->pack_start(*frame, true, true, 5);
      win->show();
    }

    {
      _positionSel.reset(new AttributeSelector(false));
      _gtkOptList->pack_start(*_positionSel, false, false);
    }

    _gtkOptList->show();
    rebuildGui();
  }

  void 
  DataSet::addGlyphs()
  {    
    if (!_context) M_throw() << "Cannot add glyphs before the Dataset is initialised";
   
    std::tr1::shared_ptr<Glyphs> glyph(new Glyphs("Glyphs", *this, _defaultGlyphType));
    _children.push_back(glyph); 
    _children.back()->init(_systemQueue);

    if (_iter)
      {
	Gtk::TreeModel::iterator child_iter = _view->_store->append(_iter->children());
	_children.back()->addViewRows(*_view, child_iter);
	_view->_view->expand_to_path(_view->_store->get_path(child_iter));
      }
  }

  void 
  DataSet::addAttribute(std::string name, int type, size_t components)
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
	const std::vector<GLfloat>& mins = iPtr->second->minVals();
	const std::vector<GLfloat>& maxs = iPtr->second->maxVals();
	if (!mins.empty() && !maxs.empty())
	  {
	    std::ostringstream os;
	    os << "[" << mins[0];
	    for (size_t i(1); i < mins.size(); ++i)
	      os << ", " << mins[i];
	    os << "]";
	    (*iter)[_attrcolumns->min] = os.str();
	    os.str("");

	    os << "[" << maxs[0];
	    for (size_t i(1); i < maxs.size(); ++i)
	      os << ", " << maxs[i];
	    os << "]";

	    (*iter)[_attrcolumns->max] = os.str();
	  }
	else
	  {
	    (*iter)[_attrcolumns->min] = "N/A";
	    (*iter)[_attrcolumns->max] = "N/A";
	  }
      }

    _positionSel->buildEntries("Position Attribute:", *this, 3, 3, 
			       Attribute::COORDINATE, 0,
			       Attribute::DEFAULT_GLYPH_POSITION);
  }

  void
  DataSet::deinit()
  {
    _positionSel.reset();
    _gtkOptList.reset();
    _attrcolumns.reset();
    _attrview.reset();
    _attrtreestore.reset();
    for (std::vector<std::tr1::shared_ptr<DataSetChild> >::iterator iPtr 
	   = _children.begin(); iPtr != _children.end(); ++iPtr)
      (*iPtr)->deinit();

    for (iterator iPtr = begin(); iPtr != end(); ++iPtr)
      iPtr->second->deinit();

    _context.reset();
    RenderObj::deinit();
  }
  
  magnet::GL::Buffer<GLfloat>&
  DataSet::getPositionBuffer()
  { 
    return _positionSel->getBuffer(); 
  }

  std::string
  DataSet::getCursorText(uint32_t objID)
  {
    std::ostringstream os;
    for (const_iterator iPtr = begin(); iPtr != end();)
      {
	size_t comps = iPtr->second->components();
	os << iPtr->first << " <"; 
	
	for (size_t i(0); i < comps-1; ++i)
	  os << (*(iPtr->second))[objID * comps + i] << ", ";
	os << (*(iPtr->second))[(objID + 1) * comps - 1] << ">";
	if (++iPtr != end())
	  os << "\n";
      }

    return os.str();
  }

  std::tr1::array<GLfloat, 4>
  DataSet::getCursorPosition(uint32_t objID)
  {
    std::vector<GLfloat> pos = _positionSel->getValue(objID);
    pos.resize(3, 0);
    std::tr1::array<GLfloat, 4> vec = {{pos[0], pos[1], pos[2], 1.0}};
    return vec;
  }

  magnet::math::Vector 
  DataSet::getDimensions() const
  {
    std::vector<GLfloat> mins = _positionSel->getMin();
    std::vector<GLfloat> maxs = _positionSel->getMax();
    
    magnet::math::Vector retVal(0,0,0);
    for (size_t i(0); (i < maxs.size()) && (i < mins.size()); ++i)
      retVal[i] = maxs[i] - mins[i];
    
    return retVal;
  }

  magnet::math::Vector 
  DataSet::getCentre() const
  {
    std::vector<GLfloat> mins = _positionSel->getMin();
    std::vector<GLfloat> maxs = _positionSel->getMax();
    
    magnet::math::Vector retVal(0,0,0);
    for (size_t i(0); (i < maxs.size()) && (i < mins.size()) && (i < 3); ++i)
      retVal[i] = (maxs[i] + mins[i]) * 0.5;
    
    return retVal;
  }
}
