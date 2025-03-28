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

#ifdef DYNAMO_visualizer
#include <coil/RenderObj/RenderObj.hpp>
#include <dynamo/base.hpp>
#include <magnet/GL/context.hpp>
#endif

namespace dynamo {
struct CoilRenderObj {
#ifdef DYNAMO_visualizer
  virtual shared_ptr<coil::RenderObj> getCoilRenderObj() const = 0;
  virtual void initRenderData(magnet::GL::Context::ContextPtr) const {}
  virtual void updateRenderData() const = 0;
#endif
};
} // namespace dynamo
