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

#include <dynamo/globals/volumetric_potential.hpp>
#include <dynamo/globals/globEvent.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/ranges/IDRangeAll.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/BC/LEBC.hpp>
#include <dynamo/ranges/IDRangeList.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/exception.hpp>
#include <fstream>
#include <cstdio>
#include <set>
#include <algorithm>

namespace dynamo {
  GlobalEvent 
  GVolumetricPotential::getEvent(const Particle &p) const {
    return GlobalEvent(p, HUGE_VAL, CELL, *this);
  }
  
  void 
  GVolumetricPotential::runEvent(Particle& p, const double dt) const {
    M_throw() << "Not implemented!";
  }

  void 
  GVolumetricPotential::outputXML(magnet::xml::XmlStream& XML) const {
    XML << magnet::xml::tag("Global")
	<< magnet::xml::attr("Type") << "VolumetricPotential"
	<< magnet::xml::attr("Name") << globName
	<< magnet::xml::attr("RawFile") << _fileName
	<< magnet::xml::attr("SampleBytes") << _sampleBytes
	<< magnet::xml::tag("Dimensions")
	<< magnet::xml::attr("x") << _ordering.getDimensions()[0]
	<< magnet::xml::attr("y") << _ordering.getDimensions()[1]
	<< magnet::xml::attr("z") << _ordering.getDimensions()[2]
	<< magnet::xml::endtag("Dimensions")
	<< magnet::xml::endtag("Global");
  }
  
  void 
  GVolumetricPotential::operator<<(const magnet::xml::Node& XML) {
    globName = XML.getAttribute("Name");
    _fileName = XML.getAttribute("RawFile");
    _sampleBytes = XML.getAttribute("SampleBytes").as<size_t>();
    auto XMLdim = XML.getNode("Dimensions");
    _ordering = Ordering(std::array<size_t, 3>{{XMLdim.getAttribute("x").as<size_t>(), XMLdim.getAttribute("y").as<size_t>(), XMLdim.getAttribute("z").as<size_t>()}});


    //Load the file data in
    std::vector<unsigned char> fileData(_ordering.size() * _sampleBytes);
    dout << "Opening " << _fileName << std::endl;
    std::ifstream file(_fileName.c_str(), std::ifstream::binary);

    if (!file.good())
      M_throw() << "Failed open the file " << _fileName;

    dout << "Loading " << _ordering.size() * _sampleBytes << " bytes" <<  std::endl;
    file.read(reinterpret_cast<char*>(fileData.data()), _ordering.size() * _sampleBytes);
    
    if (!file)
      M_throw() << "Failed reading volumetric data (read " << file.gcount() << " bytes of an expected " << _ordering.size() * _sampleBytes << "  from " << _fileName << ")";
    file.close();

    //Just pull the most significant byte in for now
    _volumeData.resize(_ordering.size());
    for (auto index : _ordering)
      _volumeData[index] = fileData[index * _sampleBytes];
  }

#ifdef DYNAMO_visualizer

  shared_ptr<coil::RenderObj>
  GVolumetricPotential::getCoilRenderObj() const
  {
    if (!_renderObj)
      {
	_renderObj.reset(new coil::RVolume(getName()));
      }
    return std::static_pointer_cast<coil::RenderObj>(_renderObj);
  }

  void 
  GVolumetricPotential::initRenderData(magnet::GL::Context::ContextPtr context) const
  {
    if (!_renderObj)
      M_throw() << "Initialising before the render object has been created";

    context->queueTask(std::bind(&coil::RVolume::loadData, _renderObj.get(), _volumeData, _ordering.getDimensions(), Vector(Sim->primaryCellSize / Sim->units.unitLength())));
  }

  void
  GVolumetricPotential::updateRenderData() const
  {
    if (!_renderObj)
      M_throw() << "Updating before the render object has been created";
  }
#endif
}
