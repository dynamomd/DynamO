/*  DYNAMO:- Event driven molecular dynamics simulator 
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

#include "nblistentry.hpp"
#include "../../dynamics/globals/neighbourList.hpp"
#include "../../base/is_simdata.hpp"
#include "../../extcode/xmlwriter.hpp"

CSCENBList::CSCENBList(const XMLNode& XML, DYNAMO::SimData* const nSim):
  CSCEntry(nSim, "ComplexNBlistEntry"),
  nblistID(-1)
{
  operator<<(XML);
}

void 
CSCENBList::operator<<(const XMLNode& XML)
{
  range.set_ptr(CRange::loadClass(XML, Sim));
  
  try 
    {
      if (strcmp(XML.getAttribute("Type"),"NeighbourList"))
	M_throw() << "Attempting to load NeighbourList from "
		  << XML.getAttribute("Type") << " entry";
  
      name = XML.getAttribute("NBListName");
    }
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in CSCENBList";
    }
}

void 
CSCENBList::initialise()
{
  try {
    nblistID = Sim->dynamics.getGlobal(name)->getID();
  }
  catch (std::exception& cep)
    {
      M_throw() << "Failed to find the global named " 
		<< name << " for the CSCENBList entry."
		<< "\n" << cep.what();
    }

  if (dynamic_cast<const CGNeighbourList* const>(Sim->dynamics.getGlobals()[nblistID].get_ptr()) == NULL)
    M_throw() << "Global named " << name << " is not a CGNeighbourList";
  
  static_cast<CGNeighbourList&>(*Sim->dynamics.getGlobals()[nblistID])
    .markAsUsedInScheduler();				 
}

void 
CSCENBList::getParticleNeighbourhood(const Particle& part, 
				     const CGNeighbourList::nbHoodFunc& func) const
{
#ifdef DYNAMO_DEBUG
  if (!isApplicable(part))
    M_throw() << "This complexNBlist entry ("
	      << name << ") is not valid for this particle (" 
	      << part.getID() << ") yet it is being used anyway!";
#endif

  static_cast<const CGNeighbourList&>(*Sim->dynamics.getGlobals()[nblistID])
    .getParticleNeighbourhood(part, func);
}

void 
CSCENBList::getParticleLocalNeighbourhood(const Particle& part, 
					  const CGNeighbourList::nbHoodFunc& func) const
{
#ifdef DYNAMO_DEBUG
  if (!isApplicable(part))
    M_throw() << "This complexNBlist entry ("
	      << name << ") is not valid for this particle (" 
	      << part.getID() << ") yet it is being used anyway!";
#endif

  static_cast<const CGNeighbourList&>(*Sim->dynamics.getGlobals()[nblistID])
    .getParticleLocalNeighbourhood(part, func);
}


void 
CSCENBList::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Type") << "NeighbourList"
      << xml::attr("NBListName")
      << Sim->dynamics.getGlobals()[nblistID]->getName()
      << range
    ;
}
