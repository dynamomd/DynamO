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

#ifndef OPMSD_H
#define OPMSD_H

#include "../outputplugin.hpp"
#include "../../datatypes/vector.hpp"
#include <vector>

class OPMSD: public OutputPlugin
{
 public:
  OPMSD(const DYNAMO::SimData*, const XMLNode&);
  ~OPMSD();

  virtual void initialise();

  virtual void eventUpdate(const IntEvent&, const C2ParticleData&) {}

  virtual void eventUpdate(const CGlobEvent&, const CNParticleData&) {}

  virtual void eventUpdate(const CLocalEvent&, const CNParticleData&) {}

  virtual void eventUpdate(const CSystem&, const CNParticleData&, const Iflt&) {}

  void output(xmlw::XmlStream &); 

  virtual OutputPlugin *Clone() const { return new OPMSD(*this); };

  Iflt calcMSD() const;

  Iflt calcStructMSD(const CTopology&) const;
  
 protected:
  
  std::vector<Vector  > initPos;
};

#endif
