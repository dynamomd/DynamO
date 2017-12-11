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
#include <coil/images/images.hpp>

namespace coil {  
  Glib::RefPtr<Gdk::Pixbuf> 
  DataSet::getIcon()
  { return coil::images::DataSet_Icon(); }

  void 
  DataSetChild::request_delete()
  { _ds.deleteChild(this); }

  void 
  DataSet::deleteChildWorker(DataSetChild* childtodelete)
  {
    for (auto iPtr = _children.begin(); iPtr != _children.end(); ++iPtr)
      if (iPtr->get() == childtodelete)
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

  magnet::math::NVector<GLfloat, 4>
  DataSetChild::getCursorPosition(uint32_t objID)
  {
    return _ds.getCursorPosition(objID); 
  }

  std::string 
  DataSetChild::getCursorText(uint32_t objID)
  { return _ds.getCursorText(objID); }

  void 
  DataSet::init(const std::shared_ptr<magnet::thread::TaskQueue>& systemQueue)
  {
    RenderObj::init(systemQueue); 
    initGtk(); 
    
    for (auto& child : _children) child->init(systemQueue);
    
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
      Gtk::Frame* frame = Gtk::manage(new Gtk::Frame("Data Set Information")); frame->show();
      _gtkOptList->pack_start(*frame, false, true, 5);
      Gtk::VBox* vbox = Gtk::manage(new Gtk::VBox); vbox->show();
      frame->add(*vbox);

      _infolabel.reset(new Gtk::Label("Points: " + boost::lexical_cast<std::string>(_N))); 
      _infolabel->show();
      vbox->pack_start(*_infolabel, false, true, 5);	
    }

    //Glyph adding mechanism
    {
      Gtk::HBox* box = Gtk::manage(new Gtk::HBox); box->show();
      _gtkOptList->pack_start(*box, false, false, 5);

      _comboPointSet.reset(new Gtk::ComboBoxText); _comboPointSet->show();
      box->pack_start(*_comboPointSet, false, false, 5);

      //Check the combo box is correct
      _comboPointSet->get_model().clear();
      for (const auto& pointset: _pointSets)
	_comboPointSet->insert(-1, pointset.first);
      _comboPointSet->set_active(0);

      Gtk::Button* btn = Gtk::manage(new Gtk::Button("Add Glyphs"));
      btn->signal_clicked().connect(sigc::mem_fun(*this, &DataSet::addGlyphs));
      btn->show();
      box->pack_start(*btn, false, false, 5);

      _comboLinkSet.reset(new Gtk::ComboBoxText); _comboLinkSet->show();
      box->pack_start(*_comboLinkSet, false, false, 5);
      //Check the combo box is correct
      _comboLinkSet->get_model().clear();
      for (const auto& linkset: _linkSets)
	_comboLinkSet->insert(-1, linkset.first);
      _comboLinkSet->set_active(0);

      btn = Gtk::manage(new Gtk::Button("Add Links"));
      btn->signal_clicked().connect(sigc::mem_fun(*this, &DataSet::addLinkGlyphs));
      btn->show();
      box->pack_start(*btn, false, false, 5);
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
   
    std::shared_ptr<Glyphs> glyph(new Glyphs(_comboPointSet->get_active_text(), *this));
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
  DataSet::addLinkGlyphs()
  {
  }

  void 
  DataSet::addAttribute(std::string name, int type, size_t components)
  {
    if (_attributes.find(name) != _attributes.end())
      M_throw() << "Trying to add an Attribute with a existing name, " << name;

    //Spinlock to force that the Data set is initialised before the attribute is created
    for (;;) if (_context) break;
    
    std::shared_ptr<Attribute> ptr(new Attribute(_N, type, components, _context));
    _attributes.insert(std::make_pair(name, ptr));

    //If we're initialised, we should rebuild the view of attributes
    if (_context) _context->queueTask(std::bind(&DataSet::rebuildGui, this));
  }
  
  void
  DataSet::rebuildGui()
  {
    _attrtreestore->clear();

    _infolabel->set_text("Points: " + boost::lexical_cast<std::string>(_N) 
			 + " Point Sets: " +boost::lexical_cast<std::string>(_pointSets.size())
			 + " Link Sets: " +boost::lexical_cast<std::string>(_linkSets.size())); 

    for (const auto& attribute : _attributes)
      {
	Gtk::TreeModel::iterator iter = _attrtreestore->append();
	(*iter)[_attrcolumns->name] = attribute.first;
	(*iter)[_attrcolumns->components] = attribute.second->components();
	const std::vector<GLfloat>& mins = attribute.second->minVals();
	const std::vector<GLfloat>& maxs = attribute.second->maxVals();
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
  DataSet::addPointSet(std::string name, const std::vector<GLuint>& data, int datatype)
  {
    _context->queueTask(std::bind(&DataSet::addPointSetWorker, this, name, data, datatype));
  }

  void
  DataSet::addPointSetWorker(std::string name, std::vector<GLuint> data, int datatype)
  {
    _pointSets[name].init(data, 1);
    _pointSets[name].glyphType = datatype;
    _comboPointSet->get_model().clear();
    for (const auto& pointset: _pointSets)
      _comboPointSet->insert(-1, pointset.first);
    _comboPointSet->set_active(0);

    std::shared_ptr<Glyphs> glyph(new Glyphs(name, *this));
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
  DataSet::deinit()
  {
    _positionSel.reset();
    _gtkOptList.reset();
    _attrcolumns.reset();
    _attrview.reset();
    _attrtreestore.reset();
    for (auto& child : _children) child->deinit();
    for (auto& attribute :_attributes) attribute.second->deinit();
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
    for (const auto& attribute : _attributes)
      {
	size_t comps = attribute.second->components();
	os << attribute.first << " <"; 
	
	for (size_t i(0); i < comps-1; ++i)
	  os << (*(attribute.second))[objID * comps + i] << ", ";
	os << (*(attribute.second))[(objID + 1) * comps - 1] << ">";
	os << "\n";
      }

    return os.str();
  }

  magnet::math::NVector<GLfloat, 4>
  DataSet::getCursorPosition(uint32_t objID)
  {
    std::vector<GLfloat> pos = _positionSel->getValue(objID);
    pos.resize(3, 0);
    return magnet::math::NVector<GLfloat, 4>{pos[0], pos[1], pos[2], 1.0};
  }

  magnet::math::Vector 
  DataSet::getMinCoord() const
  {
    magnet::math::Vector min{HUGE_VAL, HUGE_VAL, HUGE_VAL}; 
    for (const auto& child : _children)
      {
	magnet::math::Vector child_min = child->getMinCoord();
	for (size_t i(0); i < 3; ++i)
	  min[i] = std::min(min[i], child_min[i]);
      }
    return min;
  }

  magnet::math::Vector 
  DataSet::getMaxCoord() const
  {
    magnet::math::Vector max{-HUGE_VAL, -HUGE_VAL, -HUGE_VAL}; 
    for (const auto& child : _children)
      {
	magnet::math::Vector child_max = child->getMaxCoord();
	for (size_t i(0); i < 3; ++i)
	  max[i] = std::max(max[i], child_max[i]);
      }
    return max;
  }
}
