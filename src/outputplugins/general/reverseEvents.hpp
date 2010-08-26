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

#ifndef OPReverseEventsCheck_HPP
#define OPReverseEventsCheck_HPP

#include "../outputplugin.hpp"

class OPCollMatrix;

class OPReverseEventsCheck: public OutputPlugin
{
public:
  OPReverseEventsCheck(const DYNAMO::SimData*, const XMLNode&);

  ~OPReverseEventsCheck() {}

  void eventUpdate(const IntEvent&, const PairEventData&);

  void eventUpdate(const GlobalEvent&, const NEventData&);

  void eventUpdate(const LocalEvent&, const NEventData&);
  
  void eventUpdate(const System&, const NEventData&, const Iflt&);

  OutputPlugin *Clone() const { return new OPReverseEventsCheck(*this); }

  virtual void changeSystem(OutputPlugin*) {}

  virtual void initialise();

  virtual void output(xml::XmlStream&);

private:

  unsigned long lReverseEvents;
  Iflt localeps;
};

#endif
