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

#ifndef OPKEnergy_H
#define OPKEnergy_H

#include "1partproperty.hpp"
#include "../../datatypes/vector.hpp"

class OPKEnergy: public OP1PP
{
 public:
  OPKEnergy(const DYNAMO::SimData*, const XMLNode&);

  void A1ParticleChange(const ParticleEventData&);

  void A2ParticleChange(const PairEventData&);

  void stream(const Iflt&);  

  void output(xml::XmlStream &); 

  void periodicOutput();

  virtual void initialise();

  virtual OutputPlugin *Clone() const { return new OPKEnergy(*this); }

  Iflt getAvgkT() const;

  Iflt getAvgTheta() const;
     
  Iflt getAvgSqTheta() const;
  
  void changeSystem(OutputPlugin*);

  void temperatureRescale(const Iflt&);
  
  const Iflt& getCurrentkT() const { return KECurrent; }

 protected:

  Iflt InitialKE;
  Iflt KEacc;
  Iflt KEsqAcc;
  Iflt KECurrent;  
};

#endif
