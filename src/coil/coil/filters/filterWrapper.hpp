/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include <gtkmm.h>

namespace coil 
{
  template <class T, bool reqNormalDepth>
  class magnetFilterWrapper: public filter
  {
  public:
    magnetFilterWrapper() { _filter.build(); }

    inline virtual size_t type_id() { return detail::filterEnum<magnetFilterWrapper<T,reqNormalDepth> >::val; }    
    inline virtual bool isEditable() { return false; }
    inline virtual void invoke(GLuint colorTextureUnit, GLuint depthTextureUnit, size_t width, size_t height) 
    { _filter.invoke(colorTextureUnit, depthTextureUnit, width, height); }

    inline virtual bool needsNormalDepth()  { return reqNormalDepth; }

  protected:
    T _filter;
  };

}
