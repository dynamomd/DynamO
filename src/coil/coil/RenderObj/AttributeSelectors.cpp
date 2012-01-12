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

#include <coil/RenderObj/AttributeSelectors.hpp>

namespace coil {  
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

    for (DataSet::iterator iPtr = ds.begin();
	 iPtr != ds.end(); ++iPtr)
      if (((iPtr->second->getType()) & typeMask)
	  && (iPtr->second->components() >=  minComponents)
	  && (iPtr->second->components() <=  maxComponents))
	{
	  Gtk::TreeModel::Row row = *(_model->append());
	  row[_modelColumns.m_name] = iPtr->first;
	  row[_modelColumns.m_ptr] = iPtr->second;
	}
      
    typedef Gtk::TreeModel::Children::iterator iterator;
    iterator selected = _comboBox.get_model()->children().begin();

    for (iterator iPtr = _comboBox.get_model()->children().begin();
	 iPtr != _comboBox.get_model()->children().end(); ++iPtr)
      {
	std::tr1::shared_ptr<Attribute> attr_ptr = (*iPtr)[_modelColumns.m_ptr];
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
    std::tr1::shared_ptr<Attribute> ptr = (*iter)[_modelColumns.m_ptr];
 
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
	_filteredData = scalardata;
      }
      
    return _filteredData;
  }

  void 
  AttributeSelector::generateFilteredData(std::vector<GLfloat>& scalardata,
					  const std::tr1::shared_ptr<Attribute>& ptr,
					  size_t mode)
  {
    //Update the data according to what was selected
    scalardata.resize(ptr->num_elements());
    const size_t components = ptr->components();
    const std::vector<GLfloat>& attrdata = *ptr;
      
    if (mode == 1)
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

  void 
  AttributeSelector::updateGui()
  {
    _singleValueLabel.set_visible(false);
    for (size_t i(0); i < 4; ++i)
      _scalarvalues[i].hide();

    bool singlevalmode = singleValueMode();

    if (_components && singlevalmode)
      {
	_singleValueLabel.set_visible(true);
	for (size_t i(0); i < _components; ++i)
	  _scalarvalues[i].show();
      }

    _componentSelect.clear_items();
    if (singlevalmode || !_enableDataFiltering)
      _componentSelect.set_visible(false);
    else
      {
	_componentSelect.set_visible(true);

	Gtk::TreeModel::iterator iter = _comboBox.get_active();
	std::tr1::shared_ptr<Attribute> ptr = (*iter)[_modelColumns.m_ptr];

	_componentSelect.append_text("Raw Data");
	_componentSelect.append_text("Magnitude");
	_componentSelect.append_text("X");
	    
	if (ptr->components() > 1)
	  _componentSelect.append_text("Y");
	if (ptr->components() > 2)
	  _componentSelect.append_text("Z");
	if (ptr->components() > 3)
	  _componentSelect.append_text("W");
	  
	//Default to coloring using the magnitude
	_componentSelect.set_active(1);
      }

    for (size_t i(0); i < _components; ++i)
      _scalarvalues[i].set_sensitive(singlevalmode);
  }
}
