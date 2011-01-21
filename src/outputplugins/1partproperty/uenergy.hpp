/*  DYNAMO:- Event driven molecular dynamics simulator 
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

#ifndef OPUEnergy_H
#define OPUEnergy_H

#include "1partproperty.hpp"
#include "../../datatypes/vector.hpp"

class OPUEnergy: public OP1PP
{
 public:
  OPUEnergy(const DYNAMO::SimData*, const XMLNode&);

  void A1ParticleChange(const ParticleEventData&);

  void A2ParticleChange(const PairEventData&);

  void stream(const double&);  

  void output(xml::XmlStream &); 

  void periodicOutput();

  virtual void initialise();

  virtual OutputPlugin *Clone() const { return new OPUEnergy(*this); }

  double getAvgU() const;

  double getAvgSqU() const;

  double getSimU() const { return intECurrent; }

  void changeSystem(OutputPlugin*);

 protected:

  double intECurrent;
  double intEsqAcc;
  double intEAcc;
};

#endif
