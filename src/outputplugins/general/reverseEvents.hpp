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

#ifndef COPReverseEventsCheck_HPP
#define COPReverseEventsCheck_HPP

#include "../outputplugin.hpp"

class COPCollMatrix;

class COPReverseEventsCheck: public COutputPlugin
{
public:
  COPReverseEventsCheck(const DYNAMO::SimData*, const XMLNode&);

  ~COPReverseEventsCheck() {}

  void eventUpdate(const CIntEvent&, const C2ParticleData&);

  void eventUpdate(const CGlobEvent&, const CNParticleData&);
  
  void eventUpdate(const CSystem&, const CNParticleData&, const Iflt&);

  COutputPlugin *Clone() const { return new COPReverseEventsCheck(*this); }

  virtual void changeSystem(COutputPlugin*) {}

  virtual void initialise();

  virtual void output(xmlw::XmlStream&);

private:

  unsigned long lReverseEvents;
  Iflt localeps;
};

#endif
