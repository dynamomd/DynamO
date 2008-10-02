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

#ifndef COPTicker_HPP
#define COPTicker_HPP

#include "../outputplugin.hpp"

class COPTicker: public COutputPlugin
{
public:
  COPTicker(const DYNAMO::SimData*, const char*);

  //Not virtual to stop you using them! ticker() is called by the sys event
  void eventUpdate(const CIntEvent&, const C2ParticleData&) {}
  void eventUpdate(const CGlobEvent&, const CNParticleData&) {}
  void eventUpdate(const CSystem&, const CNParticleData&, const Iflt&) {}

  virtual void output(xmlw::XmlStream&) {}

  virtual void ticker() = 0;
  
  Iflt getTickerdt() const { return tickerdt; } 

  virtual void changeSystem(const DYNAMO::SimData*) 
  { I_throw() << "This plugin hasn't been prepared for changes of system\nPlugin " <<  name; }

  virtual void periodicOutput() {}

private:
  //The stream member funtion no longer exists as ticker is called on
  //events to be processed before the typical eventUpdate
  //virtual void stream(Iflt) {}  

  Iflt tickerdt;
};

#endif
