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

#include "scheduler.hpp"
#include "../dynamics/interactions/intEvent.hpp"
#include "../dynamics/globals/global.hpp"
#include "../dynamics/globals/globEvent.hpp"
#include "../base/is_simdata.hpp"
#include "../extcode/xmlwriter.hpp"
#include "../extcode/xmlParser.h"
#include <boost/foreach.hpp>
#include "include.hpp"

//A declaration of one member of a virtual class MUST exist in an
//object file somewhere! This is the one for the scheduler class
CScheduler::~CScheduler()
{}

CGlobEvent
CScheduler::getGlobEvent(const CParticle& particle) const
{
  I_cout() << "This is buggy given global cellular transitions";

  CGlobEvent col1, col2;
  BOOST_FOREACH(const smrtPlugPtr<CGlobal>& ptrGlob, Sim->Dynamics.getGlobals())
    if (ptrGlob->isInteraction(particle))
      {
	col1 = ptrGlob->getEvent(particle);
	if (col1 < col2)
	  col2 = col1;
      }
  
  return col2;
}

CScheduler* 
CScheduler::getClass(const XMLNode& XML, const DYNAMO::SimData* Sim)
{
  if ((!strcmp(XML.getAttribute("Type"),"MultList"))
      || (!strcmp(XML.getAttribute("Type"),"Cellular")))
    return new CSMultList(XML,Sim);
  else if (!strcmp(XML.getAttribute("Type"), "MultListSpecial"))
    return new CSMultListSpecial(XML,Sim);
  else if (!strcmp(XML.getAttribute("Type"), "MultListShear"))
    return new CSMultListShear(XML,Sim);
  else if (!strcmp(XML.getAttribute("Type"),"FastSingle"))
    return new CSFastSingle(XML,Sim);
  else if (!strcmp(XML.getAttribute("Type"),"GlobalCellular"))
    return new CSGlobCellular(XML,Sim);
  else 
    D_throw() << "Unknown type of Scheduler encountered";
}

xmlw::XmlStream& operator<<(xmlw::XmlStream& XML, 
			    const CScheduler& g)
{
  g.outputXML(XML);
  return XML;
}
