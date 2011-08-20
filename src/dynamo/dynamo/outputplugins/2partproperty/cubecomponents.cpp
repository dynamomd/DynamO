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


#include <dynamo/outputplugins/1partproperty/cubecomponents.hpp>
#include <dynamo/dynamics/include.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/outputplugins/0partproperty/collMatrix.hpp>
#include <boost/foreach.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
  OPCubeComp::OPCubeComp(const dynamo::SimData* tmp, const magnet::xml::Node&):
    OutputPlugin(tmp, "CubeComponents")
  {}

  void 
  OPCubeComp::initialise()
  {}

  void 
  OPCubeComp::eventUpdate(const IntEvent& iEvent, const PairEventData& pDat)
  {
    mapdata& ref = angles[mapKey(iEvent.getType(), getClassKey(iEvent))];

    std::vector<double> vals(NDIM, 0);

    for (size_t i(0); i < NDIM; ++i)
      vals[i] = pDat.rij[i] * pDat.rij[i] / Sim->dynamics.units().unitArea();

    std::sort(vals.begin(), vals.end());
  
    for (size_t i(0); i < NDIM; ++i)
      ref.angles[i].addVal(vals[i]);
  }

  void 
  OPCubeComp::eventUpdate(const GlobalEvent& globEvent, const NEventData& SDat)
  {
    BOOST_FOREACH(const PairEventData& pDat, SDat.L2partChanges)
      {
	mapdata& ref = angles[mapKey(globEvent.getType(), 
				     getClassKey(globEvent))];
	std::vector<double> vals(NDIM, 0);
      
	for (size_t i(0); i < NDIM; ++i)
	  vals[i] = pDat.rij[i] * pDat.rij[i] / Sim->dynamics.units().unitArea();
      
	std::sort(vals.begin(), vals.end());
      
	for (size_t i(0); i < NDIM; ++i)
	  ref.angles[i].addVal(vals[i]);
      }
  }

  void 
  OPCubeComp::eventUpdate(const LocalEvent& localEvent, const NEventData& SDat)
  {
    BOOST_FOREACH(const PairEventData& pDat, SDat.L2partChanges)
      {
	mapdata& ref = angles[mapKey(localEvent.getType(), 
				     getClassKey(localEvent))];
      
	std::vector<double> vals(NDIM, 0);
      
	for (size_t i(0); i < NDIM; ++i)
	  vals[i] = pDat.rij[i] * pDat.rij[i] / Sim->dynamics.units().unitArea();
      
	std::sort(vals.begin(), vals.end());
      
	for (size_t i(0); i < NDIM; ++i)
	  ref.angles[i].addVal(vals[i]);
      }
  }

  void
  OPCubeComp::eventUpdate(const System& sysEvent, const NEventData& SDat, const double&)
  {
    BOOST_FOREACH(const PairEventData& pDat, SDat.L2partChanges)
      {
	mapdata& ref = angles[mapKey(sysEvent.getType(), 
				     getClassKey(sysEvent))];
      
	std::vector<double> vals(NDIM, 0);
      
	for (size_t i(0); i < NDIM; ++i)
	  vals[i] = pDat.rij[i] * pDat.rij[i] / Sim->dynamics.units().unitArea();
      
	std::sort(vals.begin(), vals.end());

	for (size_t i(0); i < NDIM; ++i)
	  ref.angles[i].addVal(vals[i]);
      } 
  }

  void
  OPCubeComp::output(magnet::xml::XmlStream &XML)
  {
    XML << magnet::xml::tag("CubeComponents");
  
    typedef std::pair<const mapKey, mapdata> mappair;

    BOOST_FOREACH(const mappair& pair1, angles)
      {
	XML << magnet::xml::tag("Element")
	    << magnet::xml::attr("Type") 
	    << pair1.first.first
	    << magnet::xml::attr("EventName") 
	    << getName(pair1.first.second, Sim);
      
	for (size_t i(0); i < NDIM; ++i)
	  pair1.second.angles[i].outputHistogram(XML, 1.0);
      
	XML << magnet::xml::endtag("Element");
      }

    XML << magnet::xml::endtag("CubeComponents");
  }
}
