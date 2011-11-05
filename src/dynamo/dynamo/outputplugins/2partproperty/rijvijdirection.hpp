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
#include <dynamo/outputplugins/outputplugin.hpp>
#include <magnet/math/histogram.hpp>
#include <dynamo/outputplugins/eventtypetracking.hpp>
#include <map>

namespace dynamo {
  using namespace EventTypeTracking;

  class OPRijVij: public OutputPlugin
  {
  public:
    OPRijVij(const dynamo::SimData*, const magnet::xml::Node&);

    virtual void initialise();
  
    virtual void eventUpdate(const IntEvent&, const PairEventData&);

    virtual void eventUpdate(const GlobalEvent&, const NEventData&);

    virtual void eventUpdate(const LocalEvent&, const NEventData&);

    virtual void eventUpdate(const System&, const NEventData&, const double&);

    void output(magnet::xml::XmlStream &);

    virtual void changeSystem(OutputPlugin* plug) 
    { std::swap(Sim, static_cast<OPRijVij*>(plug)->Sim); }
  
  protected:
    struct mapdata
    {
      mapdata(): anglemapcount(0)
      {
	for (size_t iDim(0); iDim < NDIM; ++iDim)
	  {
	    rij[iDim] = magnet::math::Histogram(0.001);
	    vij[iDim] = magnet::math::Histogram(0.001);

	    rijcostheta[iDim].resize(2000, std::pair<size_t,double>(0, 0));
	    costhetarij[iDim].resize(1000, std::pair<size_t,double>(0, 0));
	    anglemap[iDim].resize(200, std::vector<size_t>(100, 0));
	  }
      }
    
      magnet::math::Histogram rij[NDIM];
      magnet::math::Histogram vij[NDIM];
    
      std::vector<std::pair<size_t,double> > rijcostheta[NDIM];
      std::vector<std::pair<size_t,double> > costhetarij[NDIM];
      std::vector<std::vector<size_t> > anglemap[NDIM];
      size_t anglemapcount;
    };
  
    typedef std::pair<EEventType, classKey> mapKey;

    std::map<mapKey, mapdata> rvdotacc;

    void process2PED(mapdata&, const PairEventData&);
  };
}
