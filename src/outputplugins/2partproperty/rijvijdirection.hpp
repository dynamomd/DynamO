/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef OPRijVij_H
#define OPRijVij_H

#include "../outputplugin.hpp"
#include "../../datatypes/histogram.hpp"
#include "../eventtypetracking.hpp"
#include <map>

using namespace EventTypeTracking;

class OPRijVij: public OutputPlugin
{
 public:
  OPRijVij(const DYNAMO::SimData*, const XMLNode&);

  virtual void initialise();
  
  virtual void eventUpdate(const IntEvent&, const C2ParticleData&);

  virtual void eventUpdate(const CGlobEvent&, const CNParticleData&);

  virtual void eventUpdate(const CLocalEvent&, const CNParticleData&);

  virtual void eventUpdate(const CSystem&, const CNParticleData&, const Iflt&);

  void output(xmlw::XmlStream &);

  virtual void changeSystem(OutputPlugin* plug) 
  { std::swap(Sim, static_cast<OPRijVij*>(plug)->Sim); }
  
  virtual OutputPlugin *Clone() const { return new OPRijVij(*this); };
  
 protected:
  struct mapdata
  {
    mapdata(): anglemapcount(0)
    {
      for (size_t iDim(0); iDim < NDIM; ++iDim)
	{
	  rij[iDim] = C1DHistogram(0.001);
	  vij[iDim] = C1DHistogram(0.001);

	  rijcostheta[iDim].resize(2000, std::pair<size_t,Iflt>(0, 0));
	  costhetarij[iDim].resize(1000, std::pair<size_t,Iflt>(0, 0));
	  anglemap[iDim].resize(200, std::vector<size_t>(100, 0));
	}
    }
    
    C1DHistogram rij[NDIM];
    C1DHistogram vij[NDIM];
    
    std::vector<std::pair<size_t,Iflt> > rijcostheta[NDIM];
    std::vector<std::pair<size_t,Iflt> > costhetarij[NDIM];
    std::vector<std::vector<size_t> > anglemap[NDIM];
    size_t anglemapcount;
  };
  
  typedef std::pair<EEventType, classKey> mapKey;

  std::map<mapKey, mapdata> rvdotacc;

  void process2PED(mapdata&, const C2ParticleData&);
};

#endif

