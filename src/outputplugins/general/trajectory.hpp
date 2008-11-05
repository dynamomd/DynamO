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

#ifndef COPTrajectory_HPP
#define COPTrajectory_HPP

#include "../outputplugin.hpp"
#include <fstream>

class COPTrajectory: public COutputPlugin
{
public:
  COPTrajectory(const DYNAMO::SimData*, const XMLNode&);

  COPTrajectory(const COPTrajectory&);
  
  ~COPTrajectory() {}

  void eventUpdate(const CIntEvent&, const C2ParticleData&);

  void eventUpdate(const CGlobEvent&, const CNParticleData&);

  void eventUpdate(const CLocalEvent&, const CNParticleData&);
  
  void eventUpdate(const CSystem&, const CNParticleData&, const Iflt&);

  COutputPlugin *Clone() const { return new COPTrajectory(*this); }

  virtual void changeSystem(COutputPlugin*) {}

  virtual void initialise();

  virtual void output(xmlw::XmlStream&);

private:

  void printData(const CParticle&,
		 const CParticle&) const;

  mutable std::ofstream logfile;
};

#endif
