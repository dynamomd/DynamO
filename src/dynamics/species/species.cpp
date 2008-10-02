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

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include "species.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../ranges/1range.hpp"
#include "../ranges/1RAll.hpp"
#include "../../simulation/particle.hpp"
#include "../../base/is_simdata.hpp"

CSpecies::CSpecies(DYNAMO::SimData* tmp, CRange* nr, Iflt nMass, 
		   std::string nName, unsigned int nID, std::string nIName):
  SimBase_const(tmp,"Species", IC_blue),
  mass(nMass),range(nr),spName(nName),intName(nIName),IntPtr(NULL),
  ID(nID)
{}

CSpecies::CSpecies(const XMLNode& XML, DYNAMO::SimData* tmp, unsigned int nID):
  SimBase_const(tmp,"Species", IC_blue),
  mass(1.0),range(NULL),IntPtr(NULL),
  ID(nID)
{ operator<<(XML); }

const CInteraction* 
CSpecies::getIntPtr() const
{ 
#ifdef DYNAMO_DEBUG
  if (IntPtr == NULL)
    I_throw() << "Fetching an unset interaction pointer for a species";
#endif

  return IntPtr; 
}

void
CSpecies::setIntPtr(CInteraction* nPtr)
{ IntPtr = nPtr; }

void
CSpecies::initialise()
{ 
  if (IntPtr == NULL)
    I_throw() << "Species missing a matching interaction";
}

xmlw::XmlStream& operator<<(xmlw::XmlStream& XML, const CSpecies& g)
{
  g.outputXML(XML);
  return XML;
}

void 
CSpecies::operator<<(const XMLNode& XML)
{
  range.set_ptr(CRange::loadClass(XML,Sim));
    
    try {
      mass = boost::lexical_cast<Iflt>(XML.getAttribute("Mass"));
      spName = XML.getAttribute("Name");
      intName = XML.getAttribute("IntName");
    } 
    catch (boost::bad_lexical_cast &)
      {
	I_throw() << "Failed a lexical cast in CSpecies";
      }

}

bool 
CSpecies::isSpecies(const CParticle &p1) const
{ return range->isInRange(p1); }

void 
CSpecies::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Mass") << mass
      << xmlw::attr("Name") << spName
      << xmlw::attr("IntName") << intName
      << range;
}

unsigned long 
CSpecies::getCount() const
{
  return range->size();
}
