/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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


#include "cubecomponents.hpp"
#include "../../dynamics/include.hpp"
#include <boost/foreach.hpp>
#include <boost/array.hpp>
#include "../../base/is_simdata.hpp"
#include "../0partproperty/collMatrix.hpp"

OPCubeComp::OPCubeComp(const DYNAMO::SimData* tmp, const XMLNode&):
  OutputPlugin(tmp, "CubeComponents")
{}

void 
OPCubeComp::initialise()
{}

void 
OPCubeComp::eventUpdate(const CIntEvent& iEvent, const C2ParticleData& pDat)
{
  mapdata& ref = angles[mapKey(iEvent.getType(), getClassKey(iEvent))];

  std::vector<Iflt> vals(NDIM, 0);

  for (size_t i(0); i < NDIM; ++i)
    vals[i] = pDat.rij[i] * pDat.rij[i] / Sim->dynamics.units().unitArea();

  std::sort(vals.begin(), vals.end());
  
  for (size_t i(0); i < NDIM; ++i)
    ref.angles[i].addVal(vals[i]);
}

void 
OPCubeComp::eventUpdate(const CGlobEvent& globEvent, const CNParticleData& SDat)
{
  BOOST_FOREACH(const C2ParticleData& pDat, SDat.L2partChanges)
    {
      mapdata& ref = angles[mapKey(globEvent.getType(), 
				   getClassKey(globEvent))];
      std::vector<Iflt> vals(NDIM, 0);
      
      for (size_t i(0); i < NDIM; ++i)
	vals[i] = pDat.rij[i] * pDat.rij[i] / Sim->dynamics.units().unitArea();
      
      std::sort(vals.begin(), vals.end());
      
      for (size_t i(0); i < NDIM; ++i)
	ref.angles[i].addVal(vals[i]);
    }
}

void 
OPCubeComp::eventUpdate(const CLocalEvent& localEvent, const CNParticleData& SDat)
{
  BOOST_FOREACH(const C2ParticleData& pDat, SDat.L2partChanges)
    {
      mapdata& ref = angles[mapKey(localEvent.getType(), 
				   getClassKey(localEvent))];
      
      std::vector<Iflt> vals(NDIM, 0);
      
      for (size_t i(0); i < NDIM; ++i)
	vals[i] = pDat.rij[i] * pDat.rij[i] / Sim->dynamics.units().unitArea();
      
      std::sort(vals.begin(), vals.end());
      
      for (size_t i(0); i < NDIM; ++i)
	ref.angles[i].addVal(vals[i]);
    }
}

void
OPCubeComp::eventUpdate(const CSystem& sysEvent, const CNParticleData& SDat, const Iflt&)
{
  BOOST_FOREACH(const C2ParticleData& pDat, SDat.L2partChanges)
    {
      mapdata& ref = angles[mapKey(sysEvent.getType(), 
				     getClassKey(sysEvent))];
      
      std::vector<Iflt> vals(NDIM, 0);
      
      for (size_t i(0); i < NDIM; ++i)
	vals[i] = pDat.rij[i] * pDat.rij[i] / Sim->dynamics.units().unitArea();
      
      std::sort(vals.begin(), vals.end());

      for (size_t i(0); i < NDIM; ++i)
	ref.angles[i].addVal(vals[i]);
    } 
}

void
OPCubeComp::output(xmlw::XmlStream &XML)
{
  XML << xmlw::tag("CubeComponents");
  
  typedef std::pair<const mapKey, mapdata> mappair;

  BOOST_FOREACH(const mappair& pair1, angles)
    {
      XML << xmlw::tag("Element")
	  << xmlw::attr("Type") 
	  << CIntEvent::getCollEnumName(pair1.first.first)
	  << xmlw::attr("EventName") 
	  << getName(pair1.first.second, Sim);
      
      for (size_t i(0); i < NDIM; ++i)
	pair1.second.angles[i].outputHistogram(XML, 1.0);
      
      XML << xmlw::endtag("Element");
    }

    XML << xmlw::endtag("CubeComponents");
}
