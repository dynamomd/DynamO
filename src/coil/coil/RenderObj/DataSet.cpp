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

extern const guint8 DataSet_Icon[];
extern const size_t DataSet_Icon_size;

namespace coil {  
  AttributeSelector::AttributeSelector(size_t attrnum, bool enableDataFiltering):
    _lastAttribute(NULL),
    _lastAttributeDataCount(-1),
    _lastComponentSelected(-1),
    _context(NULL),
    _components(0),
    _attrnum(attrnum),
    _enableDataFiltering(enableDataFiltering)
  {
    pack_start(_selectorRow, false, false, 5);
    _selectorRow.show();
    //Label
    _label.show();
    _selectorRow.pack_start(_label, false, false, 5);
    _context = &(magnet::GL::Context::getContext());
    //combo box
    _model = Gtk::ListStore::create(_modelColumns);
    _comboBox.set_model(_model);      
    _comboBox.pack_start(_modelColumns.m_name);
    _comboBox.show();
    _selectorRow.pack_start(_comboBox, false, false, 5);
      
    _selectorRow.pack_start(_componentSelect, false, false, 5);
      
    _singleValueLabel.show();
    _singleValueLabel.set_text("Value:");
    _singleValueLabel.set_alignment(1.0, 0.5);

    _selectorRow.pack_start(_singleValueLabel, true, true, 5);
    for (size_t i(0); i < 4; ++i)
      {
	_selectorRow.pack_start(_scalarvalues[i], false, false, 0);
	_scalarvalues[i].signal_changed()
	  .connect(sigc::bind<Gtk::Entry&>(&magnet::gtk::forceNumericEntry, _scalarvalues[i]));
	_scalarvalues[i].set_text("1.0");
	_scalarvalues[i].set_max_length(0);
	_scalarvalues[i].set_width_chars(5);	  
      }

    show();

    _comboBox.signal_changed()
      .connect(sigc::mem_fun(this, &AttributeSelector::updateGui));
  }

  AttributeColorSelector::AttributeColorSelector():
    AttributeSelector(magnet::GL::Context::vertexColorAttrIndex, true),
    _autoScaling("Autoscale to data range"),
    _lastColorMap(-1)
  {
    pack_start(_colorMapSelector, false, false, 5);
    pack_start(_autoScaling, false, false, 5);
    _autoScaling.show();
    
    _colorMapSelector.signal_changed()
      .connect(sigc::mem_fun(*this, &AttributeColorSelector::colorMapChanged));
    _autoScaling.signal_toggled()
      .connect(sigc::mem_fun(*this, &AttributeColorSelector::colorMapChanged));
  }


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

  void 
  DataSet::init(const std::tr1::shared_ptr<magnet::thread::TaskQueue>& systemQueue)
  {
    RenderObj::init(systemQueue); 
    initGtk(); 
    
    for (std::vector<std::tr1::shared_ptr<DataSetChild> >::iterator iPtr = _children.begin();
	 iPtr != _children.end(); ++iPtr)
      (*iPtr)->init(systemQueue);
    
    _overlay.init(600, 600);
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

    _gtkOptList->show();
    rebuildGui();
  }

  void 
  DataSet::addGlyphs()
  {
    std::tr1::shared_ptr<Glyphs> glyph(new Glyphs("Glyphs", *this));
    _children.push_back(glyph);
    
    if (!_context) M_throw() << "Cannot add glyphs before the Dataset is initialised";
    
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
    _overlay.deinit();
    RenderObj::deinit();
  }

  
  void 
  DataSet::glyphClicked(size_t id, Vector loc)
  {
    _selectedGlyph = id;
    _selectedGlyphLocation = loc;
  }

  void 
  DataSet::interfaceRender(const magnet::GL::Camera& camera)
  {
    if (_selectedGlyph >= 0)
      {
	_overlay.resize(camera.getWidth(), camera.getHeight());
	std::tr1::array<GLfloat, 4> vec = {{_selectedGlyphLocation[0],
					    _selectedGlyphLocation[1],
					    _selectedGlyphLocation[2],
					    1.0}};
	vec = camera.getViewMatrix() * vec;
	vec = camera.getProjectionMatrix() * vec;
	
	_overlay.setPosition((0.5 + 0.5 * vec[0] / vec[3]) * camera.getWidth(), 
			     (0.5 - 0.5 * vec[1] / vec[3]) * camera.getHeight());

	_overlay.clear();    
	for (const_iterator iPtr = begin(); iPtr != end();)
	  {
	    size_t comps = iPtr->second->components();
	    _overlay << iPtr->first << ":"; 
	
	    for (size_t i(0); i < comps; ++i)
	      _overlay << " " << (*(iPtr->second))[_selectedGlyph * comps + i];
	    if (++iPtr != end())
	      _overlay << "\n";
	  }

	_overlay.glRender();
      }
  }
}
