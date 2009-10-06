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

#ifndef COPPlateMotion_H
#define COPPlateMotion_H

#include "ticker.hpp"
#include <fstream>
#include <vector>

class COPPlateMotion: public COPTicker
{
 public:
  COPPlateMotion(const DYNAMO::SimData*, const XMLNode&);

  COPPlateMotion(const COPPlateMotion&);
  
  virtual COutputPlugin *Clone() const
  { return new COPPlateMotion(*this); }

  virtual void initialise();

  virtual void stream(Iflt) {}

  virtual void ticker();

  virtual void operator<<(const XMLNode&);
  
  virtual void eventUpdate(const CLocalEvent&, const CNParticleData&);

  virtual void eventUpdate(const CIntEvent&, const C2ParticleData&);

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
