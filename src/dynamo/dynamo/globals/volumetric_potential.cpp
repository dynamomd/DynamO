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

#include <algorithm>
#include <cstdio>
#include <dynamo/BC/LEBC.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/globals/volumetric_potential.hpp>
#include <dynamo/outputplugins/outputplugin.hpp>
#include <dynamo/ranges/IDRangeAll.hpp>
#include <dynamo/ranges/IDRangeList.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/units/units.hpp>
#include <fstream>
#include <magnet/exception.hpp>
#include <magnet/xmlreader.hpp>
#include <magnet/xmlwriter.hpp>
#include <set>

namespace dynamo {
void GVolumetricPotential::initialise(size_t id) {
  Global::initialise(id);
  GNeighbourList::reinitialise();

  for (size_t iDim = 0; iDim < NDIM; iDim++)
    _cellLatticeWidth[iDim] =
        Sim->primaryCellSize[iDim] / _ordering.getDimensions()[iDim];
  _cellDimension = _cellLatticeWidth;
  _cellOffset = Vector{0, 0, 0};

  buildCells();
}

void GVolumetricPotential::runEvent(Particle &part, const double dt) {
  // Despite the system not being streamed this must be done.  This is
  // because the scheduler and all interactions, locals and systems
  // expect the particle to be up to date.
  Sim->dynamics->updateParticle(part);

  const size_t oldCellIndex = _cellData.getCellID(part.getID());
  const int cellDirectionInt(Sim->dynamics->getSquareCellCollision3(
      part, calcPosition(oldCellIndex, part), _cellDimension));
  const size_t cellDirection = abs(cellDirectionInt) - 1;

  Event iEvent = getEvent(part);

#ifdef DYNAMO_DEBUG
  if (std::isnan(iEvent._dt))
    M_throw() << "A NAN Interaction collision time has been found when "
                 "recalculating this global"
              << iEvent;
#endif

  Sim->systemTime += iEvent._dt;
  Sim->scheduler->stream(iEvent._dt);
  Sim->stream(iEvent._dt);

  // Calculate which cell the particle might end up in
  const auto oldCellCoord = _ordering.toCoord(oldCellIndex);
  auto newCellCoord = oldCellCoord;
  newCellCoord[cellDirection] += _ordering.getDimensions()[cellDirection] +
                                 ((cellDirectionInt > 0) ? 1 : -1);
  newCellCoord[cellDirection] %= _ordering.getDimensions()[cellDirection];
  const size_t newCellIndex = _ordering.toIndex(newCellCoord);

  Vector vNorm{0, 0, 0};
  vNorm[cellDirection] = (cellDirectionInt > 0) ? -1 : 1;

  NEventData EDat;

  // Run the collision and catch the data
  Sim->dynamics->updateParticle(part);
  Vector pos(part.getPosition()), vel(part.getVelocity());
  Sim->BCs->applyBC(pos, vel);
  double potEnergyChange = 0.5 * (double(_volumeData[newCellIndex]) -
                                  double(_volumeData[oldCellIndex]));
  double arg = vel[cellDirection] * vel[cellDirection] -
               2 * potEnergyChange / Sim->species[part]->getMass(part);
  if (arg > 0) {
    EDat = ParticleEventData(part, *Sim->species[part], WALL);
    part.getVelocity()[cellDirection] *=
        std::sqrt(arg) / std::abs(part.getVelocity()[cellDirection]);
    _cellData.moveTo(oldCellIndex, newCellIndex, part.getID());
  } else
    EDat = Sim->dynamics->runPlaneEvent(part, vNorm, 1.0, 0.0);

  // Now we're past the event update the scheduler and plugins
  Sim->_sigParticleUpdate(EDat);
  Sim->scheduler->fullUpdate(part);
  for (shared_ptr<OutputPlugin> &Ptr : Sim->outputPlugins)
    Ptr->eventUpdate(iEvent, EDat);
}

void GVolumetricPotential::outputXML(magnet::xml::XmlStream &XML) const {
  XML << magnet::xml::tag("Global") << magnet::xml::attr("Type")
      << "VolumetricPotential" << magnet::xml::attr("Name") << globName
      << magnet::xml::attr("RawFile") << _fileName
      << magnet::xml::attr("SampleBytes") << _sampleBytes
      << magnet::xml::tag("Dimensions") << magnet::xml::attr("x")
      << _imageDimensions[0] << magnet::xml::attr("y") << _imageDimensions[1]
      << magnet::xml::attr("z") << _imageDimensions[2]
      << magnet::xml::endtag("Dimensions");

  if (_offset != std::array<size_t, 3>{{0, 0, 0}})
    XML << magnet::xml::tag("Offset") << magnet::xml::attr("x") << _offset[0]
        << magnet::xml::attr("y") << _offset[1] << magnet::xml::attr("z")
        << _offset[2] << magnet::xml::endtag("Offset");

  if (_ordering.getDimensions() != _imageDimensions)
    XML << magnet::xml::tag("SampleDimensions") << magnet::xml::attr("x")
        << _ordering.getDimensions()[0] << magnet::xml::attr("y")
        << _ordering.getDimensions()[1] << magnet::xml::attr("z")
        << _ordering.getDimensions()[2]
        << magnet::xml::endtag("SampleDimensions");

  XML << magnet::xml::endtag("Global");
}

void GVolumetricPotential::operator<<(const magnet::xml::Node &XML) {
  globName = XML.getAttribute("Name");
  _fileName = XML.getAttribute("RawFile");
  _sampleBytes = XML.getAttribute("SampleBytes").as<size_t>();

  // Load the dimensions of the data set (and its subset of data if
  // only processing a smaller section)
  auto XMLdim = XML.getNode("Dimensions");
  _imageDimensions =
      std::array<size_t, 3>{{XMLdim.getAttribute("x").as<size_t>(),
                             XMLdim.getAttribute("y").as<size_t>(),
                             XMLdim.getAttribute("z").as<size_t>()}};

  _offset = std::array<size_t, 3>{{0, 0, 0}};
  if (XML.hasNode("Offset")) {
    auto XMLdim = XML.getNode("Offset");
    _offset = std::array<size_t, 3>{{XMLdim.getAttribute("x").as<size_t>(),
                                     XMLdim.getAttribute("y").as<size_t>(),
                                     XMLdim.getAttribute("z").as<size_t>()}};
  }

  std::array<size_t, 3> sampleDimensions = _imageDimensions;
  if (XML.hasNode("SampleDimensions")) {
    auto XMLdim = XML.getNode("SampleDimensions");
    sampleDimensions =
        std::array<size_t, 3>{{XMLdim.getAttribute("x").as<size_t>(),
                               XMLdim.getAttribute("y").as<size_t>(),
                               XMLdim.getAttribute("z").as<size_t>()}};
  }

  Ordering fileOrdering(_imageDimensions);
  std::vector<unsigned char> fileData(fileOrdering.size() * _sampleBytes);
  dout << "Opening " << _fileName << std::endl;
  std::ifstream file(_fileName.c_str(), std::ifstream::binary);

  if (!file.good())
    M_throw() << "Failed open the file " << _fileName;
  dout << "Reading " << fileOrdering.size() * _sampleBytes
       << " bytes of data into memory" << std::endl;
  file.read(reinterpret_cast<char *>(fileData.data()),
            fileOrdering.size() * _sampleBytes);

  if (!file)
    M_throw() << "Failed reading volumetric data (read " << file.gcount()
              << " bytes of an expected " << fileOrdering.size() * _sampleBytes
              << "  from " << _fileName << ")";
  file.close();

  _ordering = Ordering(sampleDimensions);
  _volumeData.resize(_ordering.size());

  dout << "Resampling " << _ordering.size()
       << " bytes of data from the file into the simulation" << std::endl;
  if (_sampleBytes == 1) {
    if (sampleDimensions == _imageDimensions)
      std::swap(_volumeData, fileData);
    else
      for (size_t z = 0; z < sampleDimensions[2]; ++z)
        for (size_t y = 0; y < sampleDimensions[1]; ++y) {
          size_t startindex = fileOrdering.toIndex(std::array<size_t, 3>{
              {_offset[0], y + _offset[1], z + _offset[2]}});
          std::copy(fileData.begin() + startindex,
                    fileData.begin() + startindex + sampleDimensions[0],
                    _volumeData.begin() +
                        _ordering.toIndex(std::array<size_t, 3>{{0, y, z}}));
        }
  } else
    M_throw() << "Do not have an optimised loader for resampling data yet";
  dout << "Loading complete" << std::endl;
}

#ifdef DYNAMO_visualizer

shared_ptr<coil::RenderObj> GVolumetricPotential::getCoilRenderObj() const {
  if (!_renderObj) {
    _renderObj.reset(new coil::RVolume(getName()));
  }
  return std::static_pointer_cast<coil::RenderObj>(_renderObj);
}

void GVolumetricPotential::initRenderData(
    magnet::GL::Context::ContextPtr context) const {
  if (!_renderObj)
    M_throw() << "Initialising before the render object has been created";

  context->queueTask(
      std::bind(&coil::RVolume::loadData, _renderObj.get(), _volumeData,
                _ordering.getDimensions(),
                Vector{Sim->primaryCellSize / Sim->units.unitLength()}));
}

void GVolumetricPotential::updateRenderData() const {
  if (!_renderObj)
    M_throw() << "Updating before the render object has been created";
}
#endif
} // namespace dynamo
