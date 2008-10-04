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

#ifndef COPCollTicker_HPP
#define COPCollTicker_HPP

#include "../outputplugin.hpp"

class COPCollTicker: public COutputPlugin
{
public:
  COPCollTicker(const DYNAMO::SimData*, const char*, unsigned char order=100);

  virtual void eventUpdate(const CIntEvent&, const C2ParticleData&);

  virtual void eventUpdate(const CGlobEvent&, const CNParticleData&);

  virtual void eventUpdate(const CSystem&, const CNParticleData&, const Iflt&);

  virtual void output(xmlw::XmlStream&) {}

  virtual void changeSystem(const DYNAMO::SimData*) 
  { D_throw() << "This plugin hasn't been prepared for changes of system\nPlugin " <<  name; }

private:
  virtual void stream(Iflt) = 0;  
  virtual void ticker() = 0;
};

#endif
