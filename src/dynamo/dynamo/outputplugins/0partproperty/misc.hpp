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
#include "../outputplugin.hpp"
#include <ctime>
#include <time.h>

class OPMisc: public OutputPlugin
{
public:
  OPMisc(const dynamo::SimData*, const magnet::xml::Node&);
  
  virtual void initialise();
  
  virtual void eventUpdate(const IntEvent&, const PairEventData&);
  
  virtual void eventUpdate(const GlobalEvent&, const NEventData&);

  virtual void eventUpdate(const LocalEvent&, const NEventData&);
  
  virtual void eventUpdate(const System&, const NEventData&, 
			   const double&);
  
  void output(magnet::xml::XmlStream &); 
  
  void periodicOutput();
  
  virtual OutputPlugin *Clone() const 
  { return new OPMisc(*this); }
  
  double getMFT() const;
  
  void changeSystem(OutputPlugin*);

protected:
  std::time_t tstartTime;
  timespec acc_tstartTime;

  double oldSysTime;

  unsigned long dualEvents;  
  unsigned long singleEvents;
  unsigned long oldcoll;
};
