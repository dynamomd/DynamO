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

#pragma once
#include <coil/RenderObj/AttributeSelectors.hpp>

namespace coil {
  class AttributeColorSelector : public AttributeSelector
  {
  public:
    AttributeColorSelector();
    
    virtual void bindAttribute(size_t attrnum, size_t divisor = 1)
    {
      Gtk::TreeModel::iterator iter = _comboBox.get_active();

      if (singleValueMode())
	setConstantAttribute(attrnum);
      else
	{
	  std::tr1::shared_ptr<Attribute> ptr = (*iter)[_modelColumns.m_ptr];
	  //We have an attribute, check the mode the ComboBox is in,
	  //and determine if we have to do something with the data!

	  //Detect if it is in simple pass-through mode
	  if ((_componentSelect.get_visible())
	      && (_componentSelect.get_active_row_number() == 0))
	    {
	      ptr->bindAttribute(attrnum, false, divisor);
	      return;
	    }
	    
	  //Check if the data actually needs updating
	  if ((_lastAttribute != ptr.get())
	      || (_lastAttributeDataCount != ptr->getUpdateCount())
	      || (_lastComponentSelected != _componentSelect.get_active_row_number())
	      || (_lastColorMap != _colorMapSelector.getMode())
	      || _filteredData.empty())
	    {
	      std::vector<GLfloat> scalardata;
	      generateFilteredData(scalardata, ptr, _lastComponentSelected);
	      
	      if (_autoScaling.get_active() && !scalardata.empty())
		{
		  GLfloat min(scalardata[0]), max(scalardata[0]);
		  
		  for (std::vector<GLfloat>::const_iterator iPtr = scalardata.begin();
		       iPtr != scalardata.end();
		       ++iPtr)
		    {
		      min = std::min(min, *iPtr);
		      max = std::max(max, *iPtr);
		    }
		  _colorMapSelector.setRange(min, max);
		}

	      //Now convert to HSV or whatever
	      _filteredData.init(4 * scalardata.size());
	      GLfloat* data_ptr = _filteredData.map();
	      for (size_t i(0); i < scalardata.size(); ++i)
		_colorMapSelector.map(data_ptr + 4 * i, scalardata[i]);
	      
	      _filteredData.unmap();
	      
	      _lastAttribute = ptr.get();
	      _lastAttributeDataCount = ptr->getUpdateCount();
	      _lastComponentSelected = _componentSelect.get_active_row_number();
	      _lastColorMap = _colorMapSelector.getMode();
	    }

	  _filteredData.attachToAttribute(attrnum, 4, divisor);
	}
    }

  protected:
    void colorMapChanged() { _lastColorMap = -2; }

    inline virtual void updateGui()
    {
      AttributeSelector::updateGui();

      if (!singleValueMode() && _enableDataFiltering)
	{
	  //Default to coloring using the raw values
	  _componentSelect.set_active(1);
	}

      updateComponent();
    }

    void updateComponent()
    {
      if (singleValueMode() || (_componentSelect.get_active_row_number() == 0))
	{
	  _colorMapSelector.hide();
	  _autoScaling.hide();
	}
      else
	{
	  _colorMapSelector.show();
	  _autoScaling.show();
	}      
    }

    magnet::gtk::ColorMapSelector _colorMapSelector;
    Gtk::CheckButton _autoScaling;
    int _lastColorMap;
  };
}

