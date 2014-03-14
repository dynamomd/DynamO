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
#include <dynamo/globals/cells.hpp>
#include <dynamo/particle.hpp>
#include <dynamo/coilRenderObj.hpp>
#include <magnet/containers/ordering.hpp>
#ifdef DYNAMO_visualizer
# include <coil/RenderObj/Volume.hpp>
#endif 
#include <unordered_map>
#include <vector>
#include <array>

namespace dynamo {
  /*! \brief An implementation of volumetric potentials.
   */
  class GVolumetricPotential: public GCells, public CoilRenderObj
  {
  public:
  public:
    GVolumetricPotential(const magnet::xml::Node& XML, dynamo::Simulation* ptrSim):
      GCells(ptrSim, "VolumetricPotential")
    { operator<<(XML); }
    
    virtual void runEvent(Particle& p, const double dt);

    virtual void initialise(size_t id);

    virtual void operator<<(const magnet::xml::Node&);

#ifdef DYNAMO_visualizer
    virtual shared_ptr<coil::RenderObj> getCoilRenderObj() const;
    virtual void initRenderData(magnet::GL::Context::ContextPtr) const;
    virtual void updateRenderData() const;
#endif

  protected:
    virtual void outputXML(magnet::xml::XmlStream&) const;

#ifdef DYNAMO_visualizer
    mutable shared_ptr<coil::RVolume> _renderObj;
#endif

    std::string _fileName;
    size_t _sampleBytes;
    std::vector<unsigned char> _volumeData;
    std::array<size_t, 3> _imageDimensions;
    std::array<size_t, 3> _offset;
  };
}
