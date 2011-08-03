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
#include "filter.hpp"
#include <magnet/GL/shader/multiplyTexture.hpp>

namespace coil 
{
  class MultiplyFilter: public Filter
  {
  public:
    MultiplyFilter() { _filter.build(); }

    inline virtual size_t type_id() { return detail::filterEnum<MultiplyFilter>::val; }
    inline virtual bool isEditable() { return false; }
    inline virtual void invoke(GLint colorTextureUnit, size_t width, size_t height,
			       const magnet::GL::Camera& vp)
    { 
      _filter.attach();
      _filter["u_Texture0"] = colorTextureUnit;
      _filter["u_Texture1"] = 0;
      _filter.invoke(); 
      _filter.detach();
    }

    inline virtual bool needsNormalDepth()  { return false; }
  protected:
    magnet::GL::shader::MultiplyTexture _filter;
  };
}
