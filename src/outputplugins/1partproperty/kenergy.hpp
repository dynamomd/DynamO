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

#ifndef COPKEnergy_H
#define COPKEnergy_H

#include "1partproperty.hpp"
#include "../../datatypes/vector.hpp"

class COPKEnergy: public COP1PP
{
 public:
  COPKEnergy(const DYNAMO::SimData*, const XMLNode&);

  void A1ParticleChange(const C1ParticleData&);

  void A2ParticleChange(const C2ParticleData&);

  void stream(const Iflt&);  

  void output(xmlw::XmlStream &); 

  void periodicOutput();

  virtual void initialise();

  virtual COutputPlugin *Clone() const { return new COPKEnergy(*this); }

  Iflt getAvgkT() const;

  Iflt getAvgTheta() const;
     
  Iflt getAvgTheta(long) const;

  Iflt getAvgSqTheta() const;
     
  Iflt getAvgSqTheta(long) const;

  void changeSystem(COutputPlugin*);

  void temperatureRescale(const Iflt&);
  
 protected:

  Iflt sumPowerLoss;
  Iflt powerLossAcc;
  CVector<> kTacc;
  CVector<> kTsqAcc;
  CVector<> kTCurrent;
};

#endif
