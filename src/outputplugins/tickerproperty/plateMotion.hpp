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

#ifndef OPPlateMotion_H
#define OPPlateMotion_H

#include "ticker.hpp"
#include <fstream>
#include <vector>

class OPPlateMotion: public OPTicker
{
 public:
  OPPlateMotion(const DYNAMO::SimData*, const XMLNode&);

  OPPlateMotion(const OPPlateMotion&);
  
  virtual OutputPlugin *Clone() const
  { return new OPPlateMotion(*this); }

  virtual void initialise();

  virtual void stream(Iflt) {}

  virtual void ticker();

  virtual void operator<<(const XMLNode&);
  
  virtual void eventUpdate(const LocalEvent&, const NEventData&);

  virtual void eventUpdate(const IntEvent&, const PairEventData&);

  virtual void output(xmlw::XmlStream&);

 protected:
  mutable std::ofstream logfile;
  size_t plateID;
  std::string plateName;
  typedef std::pair<Iflt,std::vector<Iflt> > localEntry;
  std::vector<localEntry> localEnergyFlux;
  std::vector<localEntry> localEnergyLoss;
  Iflt partpartEnergyLoss;
  Iflt oldPlateEnergy;

};

#endif
