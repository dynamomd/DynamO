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
#include <magnet/math/quaternion.hpp>

namespace coil {
  class AttributeOrientationSelector : public AttributeSelector
  {
  public:
    AttributeOrientationSelector():
      AttributeSelector(false)
    {
      for (size_t i(0); i < 3; ++i)
	_scalarvalues[i].set_text("0.0");
      _scalarvalues[3].set_text("1.0");
    }

    virtual void bindAttribute(size_t attrnum, size_t divisor = 1)
    {
      Gtk::TreeModel::iterator iter = _comboBox.get_active();

      if (singleValueMode())
	{
	  setConstantAttribute(attrnum);
	  return;
	}
      
      std::shared_ptr<Attribute> ptr = (*iter)[_modelColumns.m_ptr];	  
      if (ptr->components() == 4)
	{
	  ptr->bindAttribute(attrnum, false, divisor);
	  return;
	}
      
      if (ptr->components() != 3)
	M_throw() << "Cannot create orientation from anything other than a 3 or 4 component Attribute";


      if ((_lastAttribute != ptr.get())
	  || (_lastAttributeDataCount != ptr->getUpdateCount())
	  || _filteredData.empty())
	{
	  _lastAttribute = ptr.get();
	  _lastAttributeDataCount = ptr->getUpdateCount();
      
	  const size_t elements = ptr->num_elements();
	  _filteredData.init(4 * elements, 4);
	  const std::vector<GLfloat>& attrdata = *ptr;
	  GLfloat* glptr = _filteredData.map();
	  
	  for (size_t i(0); i < elements; ++i)
	    {
	      //Convert the data to a quaternion
	      Vector data{attrdata[3 * i + 0], attrdata[3 * i + 1], attrdata[3 * i + 2]};
	      magnet::math::Quaternion q = magnet::math::Quaternion::fromToVector(data.normal(), Vector{0,0,1});
	      for (size_t j(0); j < 3; ++j)
		glptr[4 * i + j] = q.imaginary()[j];
	      glptr[4 * i + 3] = q.real();
	    }
	  
	  _filteredData.unmap();
	}

      _filteredData.attachToAttribute(attrnum, divisor);
    }
  };
}
