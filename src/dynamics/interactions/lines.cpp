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

#include "lines.hpp"
#include <boost/lexical_cast.hpp>
#include <cmath>
#include <iomanip>
#include "../../base/is_exception.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../../dynamics/interactions/intEvent.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../units/units.hpp"
#include "../../base/is_simdata.hpp"
#include "../2particleEventData.hpp"
#include "../BC/BC.hpp"

CILines::CILines(const DYNAMO::SimData* tmp, Iflt nd, 
		 Iflt ne, C2Range* nR):
  CInteraction(tmp, nR),
  length(nd), e(ne) {}

CILines::CILines(const XMLNode& XML, const DYNAMO::SimData* tmp):
  CInteraction(tmp,NULL)
{
  operator<<(XML);
}

void 
CILines::initialise(size_t nID)
{ ID=nID; }

void 
CILines::operator<<(const XMLNode& XML)
{ 
  if (strcmp(XML.getAttribute("Type"),"Lines"))
    D_throw() << "Attempting to load Lines from non Lines entry";
  
  range.set_ptr(C2Range::loadClass(XML,Sim));
  
  try 
    {
      length = Sim->Dynamics.units().unitLength() * 
	boost::lexical_cast<Iflt>(XML.getAttribute("Length"));
      
      e = boost::lexical_cast<Iflt>(XML.getAttribute("Elasticity"));
      
      intName = XML.getAttribute("Name");
    }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in CILines";
    }
}

Iflt 
CILines::maxIntDist() const 
{ return length; }

Iflt 
CILines::hardCoreDiam() const 
{ return 0.0; }

void 
CILines::rescaleLengths(Iflt scale) 
{ 
  length += scale * length;
}

CInteraction* 
CILines::Clone() const 
{ return new CILines(*this); }
  
CIntEvent 
CILines::getCollision(const CParticle &p1, const CParticle &p2) const 
{ 
  Sim->Dynamics.Liouvillean().updateParticlePair(p1, p2);

  CPDData colldat(*Sim, p1, p2);

  //getLineLineCollision();
  /*if (Sim->Dynamics.Liouvillean().SphereSphereInRoot(colldat, d2))
    {
#ifdef DYNAMO_OverlapTesting
      if (Sim->Dynamics.Liouvillean().sphereOverlap(colldat, d2))
	D_throw() << "Overlapping particles found" 
		  << ", particle1 " << p1.getID() << ", particle2 " 
		  << p2.getID() << "\nOverlap = " << (sqrt(colldat.r2) - sqrt(d2))/Sim->Dynamics.units().unitLength();
#endif

      return CIntEvent(p1, p2, colldat.dt, CORE, *this);
    }
  */

  return CIntEvent(p1,p2,HUGE_VAL, NONE, *this);  
}

C2ParticleData 
CILines::runCollision(const CIntEvent &event) const
{ 
  return Sim->Dynamics.Liouvillean().runLineLineCollision(); 
}
   
void 
CILines::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "Lines"
      << xmlw::attr("Length") << length / Sim->Dynamics.units().unitLength()
      << xmlw::attr("Elasticity") << e
      << xmlw::attr("Name") << intName
      << range;
}

void
CILines::checkOverlaps(const CParticle& part1, const CParticle& part2) const
{}
