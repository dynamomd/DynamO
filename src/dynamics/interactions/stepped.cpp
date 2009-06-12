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
  
//  if (isCaptured(p1, p2)) 
//    {
//      if (Sim->Dynamics.Liouvillean().SphereSphereInRoot(colldat, d2)) 
//	{
//#ifdef DYNAMO_OverlapTesting
//	  //Check that there is no overlap 
//	  if (Sim->Dynamics.Liouvillean().sphereOverlap(colldat, d2))
//	    D_throw() << "Overlapping particles found" 
//		      << ", particle1 " << p1.getID() 
//		      << ", particle2 " 
//		      << p2.getID() << "\nOverlap = " << (sqrt(colldat.r2) - sqrt(d2))/Sim->Dynamics.units().unitLength();
//#endif
//	  
//	  return CIntEvent(p1, p2, colldat.dt, CORE, *this);
//	}
//      else
//	if (Sim->Dynamics.Liouvillean().SphereSphereOutRoot(colldat, ld2))
//	  {  
//	    return CIntEvent(p1, p2, colldat.dt, WELL_OUT, *this);
//	  }
//    }
//  else if (Sim->Dynamics.Liouvillean().SphereSphereInRoot(colldat, ld2)) 
//    {
//#ifdef DYNAMO_OverlapTesting
//      if (Sim->Dynamics.Liouvillean().sphereOverlap(colldat,ld2))
//	{
//	  if (Sim->Dynamics.Liouvillean().sphereOverlap(colldat,d2))
//	    D_throw() << "Overlapping cores (but not registerd as captured) particles found in square well" 
//		      << "\nparticle1 " << p1.getID() << ", particle2 " 
//		      << p2.getID() << "\nOverlap = " 
//		      << (sqrt(colldat.r2) - sqrt(d2)) 
//	      / Sim->Dynamics.units().unitLength();
//	  else
//	    D_throw() << "Overlapping wells (but not registerd as captured) particles found" 
//		      << "\nparticle1 " << p1.getID() << ", particle2 " 
//		      << p2.getID() << "\nOverlap = " 
//		      << (sqrt(colldat.r2) - sqrt(ld2)) 
//	      / Sim->Dynamics.units().unitLength();
//	  
//	}
//#endif
//
//      return CIntEvent(p1, p2, colldat.dt, WELL_IN, *this);
//    }

  return CIntEvent(p1, p2, HUGE_VAL, NONE, *this);
}

void
CIStepped::runEvent(const CParticle& p1, 
		    const CParticle& p2,
		    const CIntEvent& iEvent) const
{
}

void
CIStepped::checkOverlaps(const CParticle& part1, const CParticle& part2) const
{
  /*
  Vector  rij = part1.getPosition() - part2.getPosition();
  Sim->Dynamics.BCs().setPBC(rij);
  Iflt r2 = rij.nrm2();

  if (isCaptured(part1, part2))
    {
      if (r2 < d2)
	I_cerr() << std::setprecision(std::numeric_limits<float>::digits10)
		 << "Possible captured overlap occured in diagnostics\n ID1=" << part1.getID() 
		 << ", ID2=" << part2.getID() << "\nR_ij^2=" 
		 << r2 / pow(Sim->Dynamics.units().unitLength(),2)
		 << "\nd^2=" 
		 << d2 / pow(Sim->Dynamics.units().unitLength(),2);

      if (r2 > lambda * lambda * d2)
	I_cerr() << std::setprecision(std::numeric_limits<float>::digits10)
		 << "Possible escaped captured pair in diagnostics\n ID1=" << part1.getID() 
		 << ", ID2=" << part2.getID() << "\nR_ij^2=" 
		 << r2 / pow(Sim->Dynamics.units().unitLength(),2)
		 << "\n(lambda * d)^2=" 
		 << ld2 / pow(Sim->Dynamics.units().unitLength(),2);
    }
  else 
    if (r2 < ld2)
      I_cerr() << std::setprecision(std::numeric_limits<float>::digits10)
	       << "Possible missed captured pair in diagnostics\n ID1=" << part1.getID() 
	       << ", ID2=" << part2.getID() << "\nR_ij^2=" 
	       << r2 / pow(Sim->Dynamics.units().unitLength(),2)
	       << "\n(lambda * d)^2=" 
	       << ld2 / pow(Sim->Dynamics.units().unitLength(),2);
  */
}
  
void 
CIStepped::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "SquareWell"
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
