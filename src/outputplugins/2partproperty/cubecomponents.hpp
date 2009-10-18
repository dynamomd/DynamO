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

#ifndef OPCubeComp_H
#define OPCubeComp_H

#include "../outputplugin.hpp"
#include "../../datatypes/histogram.hpp"
#include "../eventtypetracking.hpp"
#include <map>

using namespace EventTypeTracking;

class OPCubeComp: public OutputPlugin
{
 public:
  OPCubeComp(const DYNAMO::SimData*, const XMLNode&);

  virtual void initialise();
  
  virtual void eventUpdate(const CIntEvent&, const C2ParticleData&);

  virtual void eventUpdate(const CGlobEvent&, const CNParticleData&);

  virtual void eventUpdate(const CLocalEvent&, const CNParticleData&);

  virtual void eventUpdate(const CSystem&, const CNParticleData&, const Iflt&);

  void output(xmlw::XmlStream &);

  virtual void changeSystem(OutputPlugin* plug) { std::swap(Sim, static_cast<OPCubeComp*>(plug)->Sim); }
  
  virtual OutputPlugin *Clone() const { return new OPCubeComp(*this); };
  
 protected:
  struct mapdata
  {
    mapdata(): angles(C1DHistogram(0.01)) {}
    
    CVector<C1DHistogram> angles;
  };
  
  typedef std::pair<EEventType, classKey> mapKey;

  std::map<mapKey, mapdata> angles;
};

#endif

