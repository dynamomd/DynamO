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

#ifndef OPTrajectory_HPP
#define OPTrajectory_HPP

#include "../outputplugin.hpp"
#include <fstream>

class OPTrajectory: public OutputPlugin
{
public:
  OPTrajectory(const DYNAMO::SimData*, const XMLNode&);

  OPTrajectory(const OPTrajectory&);
  
  ~OPTrajectory() {}

  void eventUpdate(const IntEvent&, const PairEventData&);

  void eventUpdate(const GlobalEvent&, const NEventData&);

  void eventUpdate(const LocalEvent&, const NEventData&);
  
  void eventUpdate(const CSystem&, const NEventData&, const Iflt&);

  OutputPlugin *Clone() const { return new OPTrajectory(*this); }

  virtual void changeSystem(OutputPlugin*) {}

  virtual void initialise();

  virtual void output(xmlw::XmlStream&);

private:

  void printData(const size_t&,
		 const size_t&) const;

  mutable std::ofstream logfile;
};

#endif
