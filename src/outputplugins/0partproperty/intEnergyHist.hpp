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

#ifndef OPIntEnergyHist_H
#define OPIntEnergyHist_H

#include "collticker.hpp"
#include "../../datatypes/histogram.hpp"

class OPUEnergy;

class OPIntEnergyHist: public OPCollTicker
{
 public:
  OPIntEnergyHist(const DYNAMO::SimData*, const XMLNode&);

  virtual OutputPlugin *Clone() const
  { return new OPIntEnergyHist(*this); }

  virtual void initialise();

  virtual void stream(Iflt);

  virtual void ticker();

  virtual void output(xml::XmlStream&);

  virtual void changeSystem(OutputPlugin*);
  
  void operator<<(const XMLNode&);

 protected:

  C1DWeightHistogram intEnergyHist;
  const OPUEnergy* ptrOPEnergy;
  Iflt weight;
  Iflt binwidth;

};

#endif
