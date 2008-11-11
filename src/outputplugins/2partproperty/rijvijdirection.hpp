/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef COPRijVij_H
#define COPRijVij_H

#include "../outputplugin.hpp"
#include "../../datatypes/histogram.hpp"
#include "../eventtypetracking.hpp"
#include <map>

using namespace EventTypeTracking;

class COPRijVij: public COutputPlugin
{
 public:
  COPRijVij(const DYNAMO::SimData*, const XMLNode&);

  virtual void initialise();
  
  virtual void eventUpdate(const CIntEvent&, const C2ParticleData&);

  virtual void eventUpdate(const CGlobEvent&, const CNParticleData&);

  virtual void eventUpdate(const CLocalEvent&, const CNParticleData&);

  virtual void eventUpdate(const CSystem&, const CNParticleData&, const Iflt&);

  void output(xmlw::XmlStream &);

  virtual void changeSystem(COutputPlugin* plug) 
  { std::swap(Sim, static_cast<COPRijVij*>(plug)->Sim); }
  
  virtual COutputPlugin *Clone() const { return new COPRijVij(*this); };
  
 protected:
  struct mapdata
  {
    mapdata()
    {
      for (size_t iDim(0); iDim < NDIM; ++iDim)
	{
	  rij[iDim] = C1DHistogram(0.001);
	  vij[iDim] = C1DHistogram(0.001);

	  rijcostheta[iDim].resize(2000, std::pair<size_t,Iflt>(0, 0));
	  costhetarij[iDim].resize(2000, std::pair<size_t,Iflt>(0, 0));
	}
    }
    
    C1DHistogram rij[NDIM];
    C1DHistogram vij[NDIM];
    
    std::vector<std::pair<size_t,Iflt> > rijcostheta[NDIM];
    std::vector<std::pair<size_t,Iflt> > costhetarij[NDIM];
  };
  
  typedef std::pair<EEventType, classKey> mapKey;

  std::map<mapKey, mapdata> rvdotacc;

  void process2PED(mapdata&, const C2ParticleData&);
};

#endif

