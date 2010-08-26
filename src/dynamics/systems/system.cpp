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

#include "system.hpp"
#include "ghost.hpp"
#include "rescale.hpp"
#include "DSMCspheres.hpp"
#include "RingDSMC.hpp"
#include "umbrella.hpp"
#include "../../extcode/xmlParser.h"
#include "../../extcode/xmlwriter.hpp"
#include "../../simulation/particle.hpp"
#include "../interactions/intEvent.hpp"
#include "../globals/globEvent.hpp"
#include "../ranges/1RAll.hpp"

bool 
System::operator<(const IntEvent& iEvent) const
{
  return dt < iEvent.getdt();
}

bool 
System::operator<(const GlobalEvent& gEvent) const
{
  return dt < gEvent.getdt();
}

bool 
System::operator<(const System& sEvent) const
{
  return dt < sEvent.dt;
}


System::System(DYNAMO::SimData* tmp):
  SimBase(tmp, "SystemInteraction", IC_blue),
  dt(HUGE_VAL)
{
  type = NONE;
}

xml::XmlStream& operator<<(xml::XmlStream& XML, 
			    const System& g)
{
  g.outputXML(XML);
  return XML;
}

System* 
System::getClass(const XMLNode& XML, DYNAMO::SimData* Sim)
{
  if (!strcmp(XML.getAttribute("Type"),"Andersen"))
    return new CSysGhost(XML,Sim);
  else if (!strcmp(XML.getAttribute("Type"), "DSMCSpheres"))
    return new CSDSMCSpheres(XML, Sim);
  else if (!strcmp(XML.getAttribute("Type"), "Rescale"))
    return new CSysRescale(XML, Sim);
  else if (!strcmp(XML.getAttribute("Type"), "RingDSMC"))
    return new CSRingDSMC(XML, Sim);
  else if (!strcmp(XML.getAttribute("Type"), "Umbrella"))
    return new CSUmbrella(XML, Sim);
  else
    D_throw() << XML.getAttribute("Type")
	      << ", Unknown type of System Interaction encountered";
}
