/*  dynamo:- Event driven molecular dynamics simulator 
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

#pragma once
#include <dynamo/outputplugins/tickerproperty/ticker.hpp>
#include <fstream>
#include <vector>

namespace dynamo {
  class OPPlateMotion: public OPTicker
  {
  public:
    OPPlateMotion(const dynamo::SimData*, const magnet::xml::Node&);

    OPPlateMotion(const OPPlateMotion&);
  
    virtual OutputPlugin *Clone() const
    { return new OPPlateMotion(*this); }

    virtual void initialise();

    virtual void stream(double) {}

    virtual void ticker();

    virtual void operator<<(const magnet::xml::Node&);
  
    virtual void eventUpdate(const LocalEvent&, const NEventData&);

    virtual void eventUpdate(const IntEvent&, const PairEventData&);

    virtual void output(magnet::xml::XmlStream&);

  protected:
    mutable std::ofstream logfile;
    size_t plateID;
    std::string plateName;
    typedef std::pair<double,std::vector<double> > localEntry;
    std::vector<localEntry> localEnergyFlux;
    std::vector<localEntry> localEnergyLoss;
    double partpartEnergyLoss;
    double oldPlateEnergy;
    Vector momentumChange;
  };
}
