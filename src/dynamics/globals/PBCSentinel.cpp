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

#include "PBCSentinel.hpp"
#include "globEvent.hpp"
#include "../../base/is_simdata.hpp"

CGPBCSentinel::CGPBCSentinel(DYNAMO::SimData* nSim, const std::string& name):
  CGlobal(nSim, "PBCSentinel"),
  maxintdist(0)
{
  globName = name;
  I_cout() << "PBCSentinel Loaded";
}

CGPBCSentinel::CGPBCSentinel(const XMLNode &XML, DYNAMO::SimData* ptrSim):
  CGlobal(ptrSim, "PBCSentinel"),
  maxintdist(0)
{
  operator<<(XML);

  I_cout() << "PBCSentinel Loaded";
}


CGlobEvent 
CGPBCSentinel::getEvent(const CParticle& part) const
{

}

void 
CGPBCSentinel::runEvent(const CParticle& part) const
{

}

void 
CGPBCSentinel::initialise(size_t nID)
{
  ID=nID;
  
  maxintdist = Sim->Dynamics.getLongestInteraction();
}

void 
CGPBCSentinel::operator<<(const XMLNode& XML)
{
  try {
    globName = XML.getAttribute("Name");	
  }
  catch(...)
    {
      D_throw() << "Error loading CGPBCSentinel";
    }
}

void 
CGPBCSentinel::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "PBCSentinel"
      << xmlw::attr("Name") << globName;
}
