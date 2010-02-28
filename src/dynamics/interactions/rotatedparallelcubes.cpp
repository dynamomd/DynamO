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

#include "rotatedparallelcubes.hpp"
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
#include <sstream>
#include "../ranges/1range.hpp"
#include "../../schedulers/scheduler.hpp"
#include "../NparticleEventData.hpp"
#include "../../datatypes/vector.xml.hpp"

CIRotatedParallelCubes::CIRotatedParallelCubes(DYNAMO::SimData* tmp, Iflt nd, 
					       Iflt ne, const Matrix& rot,
					       C2Range* nR):
  CInteraction(tmp, nR),
  Rotation(rot),
  diameter(nd), e(ne)
{}

CIRotatedParallelCubes::CIRotatedParallelCubes(const XMLNode& XML, DYNAMO::SimData* tmp):
  CInteraction(tmp,NULL)
{
  operator<<(XML);
}

void 
CIRotatedParallelCubes::initialise(size_t nID)
{ 
  ID=nID; 
}

void 
CIRotatedParallelCubes::operator<<(const XMLNode& XML)
{ 
  if (strcmp(XML.getAttribute("Type"),"RotatedParallelCubes"))
    D_throw() << "Attempting to load RotatedParallelCubes from " 
	      << XML.getAttribute("Type") << " entry";
  
  range.set_ptr(C2Range::loadClass(XML,Sim));
  
  try 
    {
      diameter = Sim->dynamics.units().unitLength() * 
	boost::lexical_cast<Iflt>(XML.getAttribute("Diameter"));
      
      e = boost::lexical_cast<Iflt>(XML.getAttribute("Elasticity"));
            
      intName = XML.getAttribute("Name");

      ::operator<<(Rotation, XML.getChildNode("Rotation"));
    }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in CIRotatedParallelCubes";
    }
}

Iflt 
CIRotatedParallelCubes::maxIntDist() const 
{ return std::sqrt(NDIM) * diameter; }

Iflt 
CIRotatedParallelCubes::hardCoreDiam() const 
{ return diameter; }

void 
CIRotatedParallelCubes::rescaleLengths(Iflt scale) 
{ 
  diameter += scale*diameter;
}

CInteraction* 
CIRotatedParallelCubes::Clone() const 
{ return new CIRotatedParallelCubes(*this); }
  
CIntEvent 
CIRotatedParallelCubes::getEvent(const CParticle &p1, const CParticle &p2) const 
{ 
#ifdef DYNAMO_DEBUG
  if (!Sim->dynamics.getLiouvillean().isUpToDate(p1))
    D_throw() << "Particle 1 is not up to date";
  
  if (!Sim->dynamics.getLiouvillean().isUpToDate(p2))
    D_throw() << "Particle 2 is not up to date";
#endif

#ifdef DYNAMO_DEBUG
  if (p1 == p2)
    D_throw() << "You shouldn't pass p1==p2 events to the interactions!";
#endif 

  CPDData colldat(*Sim, p1, p2);

  if (Sim->dynamics.getLiouvillean().CubeCubeInRoot(colldat, diameter, Rotation))
    {
#ifdef DYNAMO_OverlapTesting
      if (Sim->dynamics.getLiouvillean().cubeOverlap(colldat, diameter, Rotation))
	D_throw() << "Overlapping particles found" 
		  << ", particle1 " << p1.getID() << ", particle2 " 
		  << p2.getID() << "\nOverlap = " 
		  << (sqrt(colldat.r2) - diameter) 
	  / Sim->dynamics.units().unitLength();
#endif

      return CIntEvent(p1, p2, colldat.dt, CORE, *this);
    }
  
  return CIntEvent(p1,p2,HUGE_VAL, NONE, *this);
}

void
CIRotatedParallelCubes::runEvent(const CParticle& p1,
				 const CParticle& p2,
				 const CIntEvent& iEvent) const
{

  ++Sim->lNColl;
    
  //Run the collision and catch the data
  C2ParticleData EDat
    (Sim->dynamics.getLiouvillean().parallelCubeColl(iEvent, e, diameter, Rotation)); 

  Sim->signalParticleUpdate(EDat);

  //Now we're past the event, update the scheduler and plugins
  Sim->ptrScheduler->fullUpdate(p1, p2);
  
  BOOST_FOREACH(smrtPlugPtr<OutputPlugin> & Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(iEvent,EDat);
}
   
void 
CIRotatedParallelCubes::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "RotatedParallelCubes"
      << xmlw::attr("Diameter") << diameter / Sim->dynamics.units().unitLength()
      << xmlw::attr("Elasticity") << e
      << xmlw::attr("Name") << intName
      << range
      << xmlw::tag("Rotation")
      << Rotation
      << xmlw::endtag("Rotation");
}

void
CIRotatedParallelCubes::checkOverlaps(const CParticle& part1, const CParticle& part2) const
{
  Vector  rij = part1.getPosition() - part2.getPosition();  
  Sim->dynamics.BCs().applyBC(rij); 
  
  if ((rij | rij) < diameter*diameter)
    I_cerr() << std::setprecision(std::numeric_limits<float>::digits10)
	     << "Possible overlap occured in diagnostics\n ID1=" << part1.getID() 
	     << ", ID2=" << part2.getID() << "\nR_ij^2=" 
	     << (rij | rij) / pow(Sim->dynamics.units().unitLength(),2)
	     << "\nd^2=" 
	     << diameter * diameter / pow(Sim->dynamics.units().unitLength(),2);
}

void 
CIRotatedParallelCubes::write_povray_desc(const DYNAMO::RGB& rgb, const size_t& specID, 
				std::ostream& os) const
{ 
  os << "#declare intrep" << ID << " = "
     << "object {\n"
     << " box {\n <" << -diameter / 2.0 << "," << -diameter / 2.0 << "," 
     << -diameter / 2.0 << ">, " 
     << " <" << diameter / 2.0 << "," << diameter / 2.0 << "," << diameter / 2.0 << "> "
     << "\n  texture { pigment { color rgb<" << rgb.R << "," << rgb.G 
     << "," << rgb.B << "> }}\n  finish { phong 0.9 phong_size 60 }\n}\n"
     << " matrix < " << Rotation(0,0) << "," << Rotation(0,1) << "," << Rotation(0,2)
     << ","<< Rotation(1,0) << "," << Rotation(1,1) << "," << Rotation(1,2)
     << ","<< Rotation(2,0) << "," << Rotation(2,1) << "," << Rotation(2,2) 
     << ",0,0,0>"
     << "\n}\n";
  
  BOOST_FOREACH(const size_t& pid, *(Sim->dynamics.getSpecies()[specID]->getRange()))
    {
      Vector  pos(Sim->vParticleList[pid].getPosition());
      Sim->dynamics.BCs().applyBC(pos);

      os << "object {\n intrep" << ID << "\n translate <"
	 << pos[0];

      for (size_t iDim(1); iDim < NDIM; ++iDim)
	os << "," << pos[iDim];
      
       os << ">\n}\n";
    }
}

