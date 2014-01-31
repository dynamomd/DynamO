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
#include <dynamo/outputplugins/outputplugin.hpp>
#include <dynamo/ranges/IDRangeList.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/exception.hpp>
#include <fstream>
#include <cstdio>
#include <set>
#include <algorithm>

namespace dynamo {
  void 
  GVolumetricPotential::initialise(size_t id) {
    Global::initialise(id);
    GNeighbourList::reinitialise();      

    for (size_t iDim = 0; iDim < NDIM; iDim++)
      _cellLatticeWidth[iDim] = Sim->primaryCellSize[iDim] / _ordering.getDimensions()[iDim];
    _cellDimension = _cellLatticeWidth;
    _cellOffset = Vector(0,0,0);

    buildCells();
  }

  
  void 
  GVolumetricPotential::runEvent(Particle& part, const double dt) const {
    //Despite the system not being streamed this must be done.  This is
    //because the scheduler and all interactions, locals and systems
    //expect the particle to be up to date.
    Sim->dynamics->updateParticle(part);

    const size_t oldCellIndex = _cellData.getCellID(part.getID());
    const int cellDirectionInt(Sim->dynamics->getSquareCellCollision3(part, calcPosition(oldCellIndex, part), _cellDimension));
    const size_t cellDirection = abs(cellDirectionInt) - 1;

    GlobalEvent iEvent(getEvent(part));

#ifdef DYNAMO_DEBUG 
    if (std::isnan(iEvent.getdt()))
      M_throw() << "A NAN Interaction collision time has been found when recalculating this global"
		<< iEvent.stringData(Sim);
#endif

    Sim->systemTime += iEvent.getdt();
    Sim->ptrScheduler->stream(iEvent.getdt());  
    Sim->stream(iEvent.getdt());

    //Calculate which cell the particle might end up in
    const auto oldCellCoord = _ordering.toCoord(oldCellIndex);
    auto newCellCoord = oldCellCoord;
    newCellCoord[cellDirection] += _ordering.getDimensions()[cellDirection] + ((cellDirectionInt > 0) ? 1 : -1);
    newCellCoord[cellDirection] %= _ordering.getDimensions()[cellDirection];
    const size_t newCellIndex = _ordering.toIndex(newCellCoord);


    Vector vNorm(0,0,0);
    vNorm[cellDirection] = (cellDirectionInt > 0) ? -1 : 1;

    NEventData EDat;

    //Run the collision and catch the data
    Sim->dynamics->updateParticle(part);    
    Vector pos(part.getPosition()), vel(part.getVelocity());
    Sim->BCs->applyBC(pos, vel);
    double potEnergyChange = 0.5 * (double(_volumeData[newCellIndex]) - double(_volumeData[oldCellIndex]));
    double arg =  vel[cellDirection] * vel[cellDirection] - 2 * potEnergyChange / Sim->species(part)->getMass(part);
    if (arg > 0)
      {
	EDat = ParticleEventData(part, *Sim->species(part), WALL);
	part.getVelocity()[cellDirection] *= std::sqrt(arg) / std::abs(part.getVelocity()[cellDirection]);
	_cellData.moveTo(oldCellIndex, newCellIndex, part.getID());
      }
    else
      EDat = Sim->dynamics->runPlaneEvent(part, vNorm, 1.0, 0.0);
    
    //Now we're past the event update the scheduler and plugins
    Sim->_sigParticleUpdate(EDat);
    Sim->ptrScheduler->fullUpdate(part);  
    for (shared_ptr<OutputPlugin> & Ptr : Sim->outputPlugins)
      Ptr->eventUpdate(iEvent, EDat);
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
