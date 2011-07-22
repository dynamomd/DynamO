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
#ifdef DYNAMO_visualizer
# include <gtkmm.h>
# include <coil/RenderObj/Arrows.hpp>
# include <memory>
# include <magnet/gtk/colorMapSelector.hpp>

class LineParticleRenderer: public coil::RArrows
{
public:
  LineParticleRenderer(size_t N, std::string name):
    RArrows(N, name)
  {
    _particleData.resize(N * 6);
  }
  
  std::vector<cl_float> _particleData;

  void sendRenderData(magnet::GL::Context& context)
  {
    context.getCLCommandQueue().enqueueWriteBuffer
      (getPointData(), false, 0, 3 * _N * sizeof(cl_float), &_particleData[0]);
    
    context.getCLCommandQueue().enqueueWriteBuffer
      (getDirectionData(), false, 0, 3 * _N * sizeof(cl_float), &_particleData[3 * _N]);
  }

protected:
};
#endif
