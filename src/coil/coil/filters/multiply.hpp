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
#include "filter.hpp"
#include <magnet/GL/multiplyTexture.hpp>

namespace coil 
{
  class MultiplyFilter: public filter
  {
  public:
    MultiplyFilter() { _filter.build(); }

    inline virtual size_t type_id() { return detail::filterEnum<MultiplyFilter>::val; }
    inline virtual bool isEditable() { return false; }
    inline virtual void invoke(GLuint colorTextureUnit, size_t width, size_t height)
    { _filter.invoke(colorTextureUnit, 0, width, height); }

    inline virtual bool needsNormalDepth()  { return false; }
  protected:
    magnet::GL::MultiplyTexture _filter;
  };
}
