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
#include "reverseEvents.hpp"
#include "../../dynamics/dynamics.hpp"
#include "../../dynamics/units/units.hpp"
#include "../../dynamics/globals/globEvent.hpp"
#include "../../dynamics/interactions/intEvent.hpp"
#include "../../dynamics/locals/localEvent.hpp"
#include <magnet/xmlwriter.hpp>

OPReverseEventsCheck::OPReverseEventsCheck(const dynamo::SimData* t1, const magnet::xml::Node&):
  OutputPlugin(t1,"ReverseEventsChecker"),
  lReverseEvents(0)
{
}

void
OPReverseEventsCheck::initialise()
{
}

void 
OPReverseEventsCheck::eventUpdate(const IntEvent& eevent, 
				   const PairEventData&)
{
  if (eevent.getdt() < 0) 
  	++lReverseEvents;
}

void 
OPReverseEventsCheck::eventUpdate(const GlobalEvent& eevent, 
				   const NEventData&)
{
  if (eevent.getdt() < 0) ++lReverseEvents;
}

void 
OPReverseEventsCheck::eventUpdate(const LocalEvent& eevent, 
				   const NEventData&)
{
  if (eevent.getdt() < 0) ++lReverseEvents;
}

void 
OPReverseEventsCheck::eventUpdate(const System&, const NEventData&, 
				   const double& dt)
{
  if (dt < 0) ++lReverseEvents;  
}

void 
OPReverseEventsCheck::output(magnet::xml::XmlStream& XML)
{
  dout << "Reverse Event Count " << lReverseEvents << std::endl;

  XML << magnet::xml::tag("ReverseEvents")
      << magnet::xml::attr("Count") << lReverseEvents
      << magnet::xml::endtag("ReverseEvents");

}
