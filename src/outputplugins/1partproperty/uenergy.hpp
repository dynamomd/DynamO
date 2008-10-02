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

#ifndef COPUEnergy_H
#define COPUEnergy_H

#include "1partproperty.hpp"
#include "../../datatypes/vector.hpp"

class COPUEnergy: public COP1PP
{
 public:
  COPUEnergy(const DYNAMO::SimData*);

  void A1ParticleChange(const C1ParticleData&);

  void A2ParticleChange(const C2ParticleData&);

  void stream(Iflt);  

  void output(xmlw::XmlStream &); 

  void periodicOutput();

  virtual void initialise();

  virtual COutputPlugin *Clone() const { return new COPUEnergy(*this); }

  Iflt getAvgU() const;

  Iflt getAvgSqU() const;

  Iflt getSimU() const { return intECurrent; }

  void changeSystem(COutputPlugin*);

  void temperatureRescale(const double&) {}
  
 protected:

  Iflt intECurrent;
  Iflt intEsqAcc;
  Iflt intEAcc;
};

#endif
