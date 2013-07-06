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

#include <coil/RenderObj/AttributeSelectors.hpp>
#include <coil/RenderObj/AttributeColorSelector.hpp>

namespace coil {  
  AttributeColorSelector::AttributeColorSelector():
    AttributeSelector(true),
    _autoScaling("Autoscale to data range"),
    _lastColorMap(-1)
  {
    pack_start(_colorMapSelector, false, false, 5);
    pack_start(_autoScaling, false, false, 5);

    _autoScaling.set_active(true);
    _autoScaling.show();
    
    _colorMapSelector.signal_changed()
      .connect(sigc::mem_fun(*this, &AttributeColorSelector::colorMapChanged));
    _autoScaling.signal_toggled()
      .connect(sigc::mem_fun(*this, &AttributeColorSelector::colorMapChanged));

    _selectorRow.pack_end(_colorButton,false, false, 5);
    _colorButton.show();
    _colorButton.signal_color_set().connect(sigc::mem_fun(*this, &AttributeColorSelector::colorButtonUsed));

    _componentSelect.signal_changed()
      .connect(sigc::mem_fun(this, &AttributeColorSelector::updateComponent));

    //Change the handler for the single value boxes over to one which
    //will update the color button.
    for (size_t i(0); i < 4; ++i)
      _scalarvalues[i].signal_changed()
	.connect(sigc::mem_fun(*this, &AttributeColorSelector::colorValuesChanged));
    colorValuesChanged();
  }

  void 
  AttributeColorSelector::colorButtonUsed() {
    Gdk::Color color = _colorButton.get_color();
    _scalarvalues[0].set_text(boost::lexical_cast<std::string>(color.get_red()/65535.0));
    _scalarvalues[1].set_text(boost::lexical_cast<std::string>(color.get_green()/65535.0));
    _scalarvalues[2].set_text(boost::lexical_cast<std::string>(color.get_blue()/65535.0));
  }
  
  void  
  AttributeColorSelector::colorValuesChanged()
  {
    magnet::gtk::forceNumericEntry(_scalarvalues + 0);
    magnet::gtk::forceNumericEntry(_scalarvalues + 1);
    magnet::gtk::forceNumericEntry(_scalarvalues + 2);
    magnet::gtk::forceNumericEntry(_scalarvalues + 3);
    Gdk::Color color = _colorButton.get_color();

    double colors[3] = {1,1,1};
    for (size_t i(0); i < 3; ++i)
      try {
	colors[i] = boost::lexical_cast<double>(_scalarvalues[i].get_text());
      } catch (...) {}

    color.set_rgb_p(colors[0],colors[1], colors[2]);
    _colorButton.set_color(color);
  }
  

  void 
  AttributeColorSelector::updateGui()
  {
    AttributeSelector::updateGui();

    bool singlevalmode = singleValueMode();

    _colorButton.set_visible(_components && singlevalmode);

    for (size_t i(0); i < _components; ++i)
      _scalarvalues[i].set_visible(false);
      
    if (!singleValueMode() && _enableDataFiltering)
      {
	//Default to coloring using the raw values
	_componentSelect.set_active(1);
      }
    
    updateComponent();
  }

  AttributeSelector::AttributeSelector(bool enableDataFiltering):
    _lastAttribute(NULL),
    _lastAttributeDataCount(-1),
    _lastComponentSelected(-1),
    _components(0),
    _enableDataFiltering(enableDataFiltering)
  {
    pack_start(_selectorRow, false, false, 5);
    _selectorRow.show();
    //Label
    _label.show();
    _selectorRow.pack_start(_label, false, false, 5);
    _context = magnet::GL::Context::getContext();
    //combo box
    _model = Gtk::ListStore::create(_modelColumns);
    _comboBox.set_model(_model);      
    _comboBox.pack_start(_modelColumns.m_name);
    _comboBox.show();
    _selectorRow.pack_start(_comboBox, false, false, 5);
      
    _selectorRow.pack_start(_componentSelect, false, false, 5);


    {
      _singleValueLabel.show();
      _singleValueLabel.set_text("Value:");
      _singleValueLabel.set_alignment(1.0, 0.5);
      
      _selectorRow.pack_start(_singleValueLabel, true, true, 5);
      for (size_t i(0); i < 4; ++i)
	{
	  _selectorRow.pack_start(_scalarvalues[i], false, false, 0);
	  _scalarvalues[i].signal_changed()
	    .connect(sigc::bind(&magnet::gtk::forceNumericEntry, _scalarvalues + i));
	  _scalarvalues[i].set_text("1.0");
	  _scalarvalues[i].set_max_length(0);
	  _scalarvalues[i].set_width_chars(5);	  
	}
    }

    show();

    _comboBox.signal_changed()
      .connect(sigc::mem_fun(this, &AttributeSelector::updateGui));
  }

  void 
  AttributeSelector::buildEntries(std::string name, DataSet& ds, size_t minComponents, size_t maxComponents, 
				  int typeMask, size_t components, int defaultMask)
  {
    _components = components;
    _label.set_text(name);
    _model->clear();
      
    updateGui();

    if (_components)
      {
	Gtk::TreeModel::Row row = *(_model->append());
	row[_modelColumns.m_name] = "Single Value";
      }

    for (auto& data : ds.getAttributes())
      if (((data.second->getType()) & typeMask) && (data.second->components() >=  minComponents) && (data.second->components() <=  maxComponents))
	{
	  Gtk::TreeModel::Row row = *(_model->append());
	  row[_modelColumns.m_name] = data.first;
	  row[_modelColumns.m_ptr] = data.second;
	}
      
    typedef Gtk::TreeModel::Children::iterator iterator;
    iterator selected = _comboBox.get_model()->children().begin();

    for (auto iPtr = _comboBox.get_model()->children().begin(); iPtr != _comboBox.get_model()->children().end(); ++iPtr)
      {
	std::shared_ptr<Attribute> attr_ptr = (*iPtr)[_modelColumns.m_ptr];
	if ((attr_ptr) && (attr_ptr->getType() & defaultMask))
	  {
	    selected = iPtr;
	    break;
	  }
      }

    _comboBox.set_active(selected);
  }

  magnet::GL::Buffer<GLfloat>& 
  AttributeSelector::getBuffer()
  {
    if (singleValueMode())
      M_throw() << "Cannot get the attribute buffer when in single value mode.";

    Gtk::TreeModel::iterator iter = _comboBox.get_active();
    std::shared_ptr<Attribute> ptr = (*iter)[_modelColumns.m_ptr];
 
    if ((!_componentSelect.get_visible())
	|| (_componentSelect.get_active_row_number() == 0))
      return ptr->getBuffer();

    //Check if the data actually needs updating
    if ((_lastAttribute != ptr.get())
	|| (_lastAttributeDataCount != ptr->getUpdateCount())
	|| (_lastComponentSelected != _componentSelect.get_active_row_number())
	|| _filteredData.empty())
      {
	_lastAttribute = ptr.get();
	_lastAttributeDataCount = ptr->getUpdateCount();
	_lastComponentSelected = _componentSelect.get_active_row_number();
	  
	std::vector<GLfloat> scalardata;
	generateFilteredData(scalardata, ptr, _lastComponentSelected);
	_filteredData.init(scalardata, 1);
      }
      
    return _filteredData;
  }

  void 
  AttributeSelector::generateFilteredData(std::vector<GLfloat>& scalardata,
					  const std::shared_ptr<Attribute>& ptr,
					  int mode)
  {
    //Update the data according to what was selected
    scalardata.resize(ptr->num_elements());
    const size_t components = ptr->components();
    const std::vector<GLfloat>& attrdata = *ptr;
      
    if (mode <= 1)
      //Magnitude calculation
      {
	for (size_t i(0); i < scalardata.size(); ++i)
	  {
	    scalardata[i] = 0;
	    for (size_t j(0); j < components; ++j)
	      {
		GLfloat val = attrdata[i * components + j];
		scalardata[i] += val * val;
	      }
	    scalardata[i] = std::sqrt(scalardata[i]);
	  }
      }
    else
      {
	//Component wise selection
	size_t component = mode - 2;
#ifdef COIL_DEBUG
	if (component >= components)
	  M_throw() << "Trying to filter an invalid component";
#endif
	for (size_t i(0); i < scalardata.size(); ++i)
	  scalardata[i] = attrdata[i * components + component];
      }
  }

  void 
  AttributeSelector::setConstantAttribute(size_t attr)
  {
    _context->disableAttributeArray(attr);

    GLfloat val[4] = {1,1,1,1};

    for (size_t i(0); i < 4; ++i)
      try {
	val[i] = boost::lexical_cast<double>(_scalarvalues[i].get_text());
      } catch (...) {}
      
    _context->setAttribute(attr, val[0], val[1], val[2], val[3]);
  }


  std::vector<GLfloat>
  AttributeSelector::getMin()
  {
    if (singleValueMode())
      {
	std::vector<GLfloat> _retVal(4);
	for (size_t i(0); i < 4; ++i)
	  try {
	    _retVal[i] = boost::lexical_cast<double>(_scalarvalues[i].get_text());
	  } catch (...) {}
	
	return _retVal;
      }
    
    Gtk::TreeModel::iterator iter = _comboBox.get_active();
    if (!iter) return std::vector<GLfloat>();
    
    std::shared_ptr<Attribute> ptr = (*iter)[_modelColumns.m_ptr];
    if (!ptr) return std::vector<GLfloat>();
        
    return ptr->minVals();
  }

  std::vector<GLfloat>
  AttributeSelector::getMax()
  {
    if (singleValueMode())
      return getMin();

    Gtk::TreeModel::iterator iter = _comboBox.get_active();
    if (!iter) return std::vector<GLfloat>();
    
    std::shared_ptr<Attribute> ptr = (*iter)[_modelColumns.m_ptr];
    if (!ptr) return std::vector<GLfloat>();
        
    return ptr->maxVals();
  }

  void 
  AttributeSelector::updateGui()
  {
    _singleValueLabel.set_visible(false);
    for (size_t i(0); i < 4; ++i)
      _scalarvalues[i].set_visible(false);

    bool singlevalmode = singleValueMode();

    if (_components && singlevalmode)
      {
	_singleValueLabel.set_visible(true);
	for (size_t i(0); i < _components; ++i)
	  _scalarvalues[i].set_visible(true);
      }

    _componentSelect.clear_items();
    if (singlevalmode || !_enableDataFiltering)
      _componentSelect.set_visible(false);
    else
      {
	_componentSelect.set_visible(true);

	Gtk::TreeModel::iterator iter = _comboBox.get_active();
	std::shared_ptr<Attribute> ptr = (*iter)[_modelColumns.m_ptr];

	_componentSelect.append_text("Raw Data");
	_componentSelect.append_text("Magnitude");
	_componentSelect.append_text("X");
	    
	if (ptr->components() > 1)
	  _componentSelect.append_text("Y");
	if (ptr->components() > 2)
	  _componentSelect.append_text("Z");
	if (ptr->components() > 3)
	  _componentSelect.append_text("W");
	  
	//Default to coloring using the raw values
	_componentSelect.set_active(0);
      }
  }
}
