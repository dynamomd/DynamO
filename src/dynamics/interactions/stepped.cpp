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

#include "stepped.hpp"
#include <boost/lexical_cast.hpp>
#include <cmath>
#include "../../base/is_exception.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../BC/BC.hpp"
#include "../dynamics.hpp"
#include "../units/units.hpp"
#include "../globals/global.hpp"
#include "../../simulation/particle.hpp"
#include "../interactions/intEvent.hpp"
#include "../species/species.hpp"
#include "../2particleEventData.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../../base/is_simdata.hpp"
#include "../../schedulers/scheduler.hpp"
#include "../NparticleEventData.hpp"
#include <iomanip>
#include <boost/math/special_functions/pow.hpp>

CIStepped::CIStepped(DYNAMO::SimData* tmp, 
		     const std::vector<steppair>& vec, C2Range* nR):
  CIMultiCapture(tmp,nR),
  steps(vec)
{}

CIStepped::CIStepped(const XMLNode& XML, DYNAMO::SimData* tmp):
  CIMultiCapture(tmp, NULL) //A temporary value!
{
  operator<<(XML);
}

void 
CIStepped::operator<<(const XMLNode& XML)
{
  if (strcmp(XML.getAttribute("Type"),"Stepped"))
    D_throw() << "Attempting to load Stepped from non Stepped entry";
  
  range.set_ptr(C2Range::loadClass(XML,Sim));
  
  try {
    intName = XML.getAttribute("Name");

    if (XML.nChildNode("Step"))
      {
	int xml_iter = 0;
	long counter = XML.nChildNode("Step");
	for (long i = 0; i < counter; ++i)
	  {
	    XMLNode browseNode = XML.getChildNode("Step",&xml_iter);
	    steps.push_back(steppair(boost::lexical_cast<Iflt>
				     (browseNode.getAttribute("R"))
				     * Sim->Dynamics.units().unitLength(),
				     boost::lexical_cast<Iflt>
				     (browseNode.getAttribute("E"))
				     * Sim->Dynamics.units().unitEnergy()
				     ));
	    
	  }
      }
    else
      D_throw() << "No steppings defined for stepped potential " 
		<< intName;
    
    std::sort(steps.rbegin(), steps.rend());

    CIMultiCapture::loadCaptureMap(XML);   
  }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in CIStepped";
    }
}

CInteraction* 
CIStepped::Clone() const 
{ return new CIStepped(*this); }

Iflt 
CIStepped::hardCoreDiam() const 
{ return steps.back().first; }

Iflt 
CIStepped::maxIntDist() const 
{ return steps.front().first; }

void 
CIStepped::rescaleLengths(Iflt scale) 
{ 
  BOOST_FOREACH(steppair& p, steps)
    p.first += scale * p.first;
}

void 
CIStepped::initialise(size_t nID)
{
  ID = nID;
  CIMultiCapture::initCaptureMap();
}

int 
CIStepped::captureTest(const CParticle& p1, const CParticle& p2) const
{
  Vector  rij = p1.getPosition() - p2.getPosition();
  Sim->Dynamics.BCs().setPBC(rij);
  
  Iflt r = rij.nrm();

  for (size_t i(0); i < steps.size(); ++i)
    if (r > steps[i].first) return i;

  return steps.size() - 1;
}

Iflt 
CIStepped::getInternalEnergy() const 
{ 
  //Once the capture maps are loaded just iterate through that determining energies
  Iflt Energy = 0.0;

  typedef std::pair<const std::pair<size_t, size_t>, int> locpair;

  BOOST_FOREACH(const locpair& IDs, captureMap)
    Energy += steps[IDs.second - 1].second;
  
  return Energy; 
}

CIntEvent
CIStepped::getEvent(const CParticle &p1, 
		    const CParticle &p2) const 
{
  
#ifdef DYNAMO_DEBUG
  if (!Sim->Dynamics.Liouvillean().isUpToDate(p1))
    D_throw() << "Particle 1 is not up to date";
  
  if (!Sim->Dynamics.Liouvillean().isUpToDate(p2))
    D_throw() << "Particle 2 is not up to date";

  if (p1 == p2)
    D_throw() << "You shouldn't pass p1==p2 events to the interactions!";
#endif 

  CPDData colldat(*Sim, p1, p2);
  
  const_cmap_it capstat = getCMap_it(p1,p2);

  if (capstat == captureMap.end())
    {
      //Not captured, test for capture
      if (Sim->Dynamics.Liouvillean().SphereSphereInRoot
	  (colldat, boost::math::pow<2>(steps.front().first)))
	{
#ifdef DYNAMO_OverlapTesting
	  //Check that there is no overlap 
	  if (Sim->Dynamics.Liouvillean().sphereOverlap
	      (colldat, boost::math::pow<2>(steps.front().first)))
	    D_throw() << "Overlapping particles found" 
		      << ", particle1 " << p1.getID() 
		      << ", particle2 " 
		      << p2.getID() << "\nOverlap = " 
		      << (sqrt(colldat.r2) - steps.front().first)
	      /Sim->Dynamics.units().unitLength();
#endif
	  
	  return CIntEvent(p1, p2, colldat.dt, CAPTURE, *this);
	}
    }
  else
    {
      //Within the potential, look for further capture or release
      if (Sim->Dynamics.Liouvillean().SphereSphereInRoot
	  (colldat, boost::math::pow<2>(steps[capstat->second].first)))
	{
#ifdef DYNAMO_OverlapTesting
	  //Check that there is no overlap 
	  if (Sim->Dynamics.Liouvillean().sphereOverlap
	      (colldat, boost::math::pow<2>(steps[capstat->second].first)))
	    D_throw() << "Overlapping particles found" 
		      << ", particle1 " << p1.getID() 
		      << ", particle2 " 
		      << p2.getID() << "\nOverlap = " 
		      << (sqrt(colldat.r2) - steps[capstat->second].first)
	      /Sim->Dynamics.units().unitLength();
#endif
	  
	  return CIntEvent(p1, p2, colldat.dt, 
			   (capstat->second + 1 == static_cast<int>(steps.size())) 
			   ? CORE : WELL_IN , *this);
	}
      else if (Sim->Dynamics.Liouvillean().SphereSphereOutRoot
	       (colldat, boost::math::pow<2>(steps[capstat->second-1].first)))
	return CIntEvent(p1, p2, colldat.dt, 
			 (capstat->second == 1) ? RELEASE : WELL_OUT, *this);
    }

  return CIntEvent(p1, p2, HUGE_VAL, NONE, *this);
}

void
CIStepped::runEvent(const CParticle& p1, 
		    const CParticle& p2,
		    const CIntEvent& iEvent) const
{
  ++Sim->lNColl;

  switch (iEvent.getType())
    {
    case CORE:
      {
	C2ParticleData retVal(Sim->Dynamics.Liouvillean().SmoothSpheresColl
			      (iEvent, 1.0, boost::math::pow<2>(steps.back().first), CORE));
	Sim->signalParticleUpdate(retVal);
	
	Sim->ptrScheduler->fullUpdate(p1, p2);
	
	BOOST_FOREACH(smrtPlugPtr<COutputPlugin> & Ptr, Sim->outputPlugins)
	  Ptr->eventUpdate(iEvent, retVal);

	break;
      }
    case CAPTURE:
      {
	C2ParticleData retVal(Sim->Dynamics.Liouvillean()
			      .SphereWellEvent(iEvent, -steps.front().second, 
					       boost::math::pow<2>(steps.front().first)));
	
	if (retVal.getType() != BOUNCE)
	  addToCaptureMap(p1,p2);
	  
	Sim->ptrScheduler->fullUpdate(p1, p2);
	Sim->signalParticleUpdate(retVal);
	
	BOOST_FOREACH(smrtPlugPtr<COutputPlugin> & Ptr, Sim->outputPlugins)
	  Ptr->eventUpdate(iEvent, retVal);

	break;
      }
    case RELEASE:
      {
	C2ParticleData retVal(Sim->Dynamics.Liouvillean()
			      .SphereWellEvent(iEvent, steps.front().second, 
					       boost::math::pow<2>(steps.front().first)));
	
	if (retVal.getType() != BOUNCE)
	  delFromCaptureMap(p1,p2);
	  
	Sim->ptrScheduler->fullUpdate(p1, p2);
	Sim->signalParticleUpdate(retVal);
	
	BOOST_FOREACH(smrtPlugPtr<COutputPlugin> & Ptr, Sim->outputPlugins)
	  Ptr->eventUpdate(iEvent, retVal);

	break;
      }
    case WELL_OUT:
      {
	cmap_it capstat = getCMap_it(p1,p2);
	
	C2ParticleData retVal(Sim->Dynamics.Liouvillean().SphereWellEvent
			      (iEvent, -(steps[capstat->second - 2].second 
					 - steps[capstat->second - 1].second), 
			       boost::math::pow<2>(steps[capstat->second -1].first)));
	
	if (retVal.getType() != BOUNCE)
	  --(capstat->second);

	Sim->signalParticleUpdate(retVal);

	Sim->ptrScheduler->fullUpdate(p1, p2);
	
	BOOST_FOREACH(smrtPlugPtr<COutputPlugin> & Ptr, Sim->outputPlugins)
	  Ptr->eventUpdate(iEvent, retVal);
	break;
      }
    case WELL_IN:
      {
	cmap_it capstat = getCMap_it(p1,p2);
	
	C2ParticleData retVal(Sim->Dynamics.Liouvillean().SphereWellEvent
			      (iEvent, -(steps[capstat->second].second 
					 - steps[capstat->second - 1].second), 
			       boost::math::pow<2>(steps[capstat->second].first)));
	
	if (retVal.getType() != BOUNCE)
	  ++(capstat->second);

	Sim->signalParticleUpdate(retVal);

	Sim->ptrScheduler->fullUpdate(p1, p2);
	
	BOOST_FOREACH(smrtPlugPtr<COutputPlugin> & Ptr, Sim->outputPlugins)
	  Ptr->eventUpdate(iEvent, retVal);
	break;
      }
    default:
      D_throw() << "Unknown collision type";
    } 
}

void
CIStepped::checkOverlaps(const CParticle& part1, const CParticle& part2) const
{
  const_cmap_it capstat = getCMap_it(part1,part2);

  if (captureTest(part1,part2) != capstat->second)
    I_cerr() << "Particle " << part1.getID() << " and Particle " << part2.getID()
	     << "\nFailing as captureTest gives " << captureTest(part1,part2)
	     << "\nAnd recorded value is " << capstat->second;
}
  
void 
CIStepped::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "Stepped"
      << xmlw::attr("Name") << intName
      << range;

  BOOST_FOREACH(const steppair& s, steps)
    XML << xmlw::tag("Step")
	<< xmlw::attr("R") 
	<< s.first / Sim->Dynamics.units().unitLength()
	<< xmlw::attr("E")
	<< s.second / Sim->Dynamics.units().unitEnergy()
	<< xmlw::endtag("Step");
  
  CIMultiCapture::outputCaptureMap(XML);  
}

void 
CIStepped::write_povray_desc(const DYNAMO::RGB& rgb, 
				const size_t& specID, 
				std::ostream& os) const
{
  /*
  os << "#declare intrep" << ID << "center = " 
     << "sphere {\n <0,0,0> " 
     << diameter / 2.0
     << "\n texture { pigment { color rgb<" << rgb.R << "," << rgb.G 
     << "," << rgb.B << "> }}\nfinish { phong 0.9 phong_size 60 }\n}\n"
     << "#declare intrep" << ID << "well = sphere {\n <0,0,0> " << diameter * lambda * 0.5
     << "\n texture { pigment { color rgbt <1,1,1,0.9> }}\n}\n";

  BOOST_FOREACH(const size_t& part, *(Sim->Dynamics.getSpecies()[specID]->getRange()))
    {
      Vector  pos(Sim->vParticleList[part].getPosition());
      Sim->Dynamics.BCs().setPBC(pos);
      
      os << "object {\n intrep" << ID << "center\n translate < "
	 << pos[0] << ", " << pos[1] << ", " << pos[2] << ">\n}\n";
    }
  
  os << "merge {\n";
  BOOST_FOREACH(const size_t& part, *(Sim->Dynamics.getSpecies()[specID]->getRange()))
    {
      Vector  pos(Sim->vParticleList[part].getPosition());
      Sim->Dynamics.BCs().setPBC(pos);

      os << "object {\n intrep" << ID << "well\n translate < "
	 << pos[0] << ", " << pos[1] << ", " << pos[2] << ">\n}\n";
    }
  
  os << "}\n";
  */
}
