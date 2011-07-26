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
#pragma once
#include <coil/RenderObj/RenderObj.hpp>
#include <coil/RenderObj/Attribute.hpp>

namespace coil {
  class DataSetChild: public RenderObj
  {
  public:
    
  protected:
  };

  /*! \brief A container class for a collection of \ref Attribute
   * instances forming a dataset, and any active filters/glyphs or any
   * other type derived from \ref DataSetChild.
   */
  class DataSet: public RenderObj
  {
  public:
    DataSet(std::string name): RenderObj(name) {}
    
    /** @name The host code interface. */
    /**@{*/

    /*! \brief Method to add attributes to the DataSet.
     */
    void addAttribute(const Attribute& attr)
    { _attributes.push_back(attr); }
    
    /**@}*/
    
    /*! \brief Returns the list of attributes.
     */
    std::vector<Attribute>& getAttributes()
    { return _attributes; }

    virtual void clTick(const magnet::GL::Camera& cam) 
    {
      for (std::vector<DataSetChild>::iterator iPtr = _children.begin();
	   iPtr != _children.end(); ++iPtr)
	iPtr->clTick(cam);
    }

    virtual void glRender(magnet::GL::FBO& fbo, const magnet::GL::Camera& cam)
    {
      for (std::vector<DataSetChild>::iterator iPtr = _children.begin();
	   iPtr != _children.end(); ++iPtr)
	iPtr->glRender(fbo, cam);
    }

  protected:
    std::vector<Attribute> _attributes;
    std::vector<DataSetChild> _children;
  };
}
