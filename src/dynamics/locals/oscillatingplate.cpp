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

#include "oscillatingplate.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "localEvent.hpp"
#include "../NparticleEventData.hpp"
#include "../overlapFunc/CubePlane.hpp"
#include "../units/units.hpp"
#include "../../datatypes/vector.xml.hpp"
#include "../../schedulers/scheduler.hpp"


CLOscillatingPlate::CLOscillatingPlate(DYNAMO::SimData* nSim, 
				       Iflt nx0, Iflt nxi, 
				       Iflt nomega0, Iflt nsigma, Iflt ne,
				       std::string nname, CRange* nRange):
  CLocal(nRange, nSim, "LocalWall"),
  x0(nx0), xi(nxi), omega0(nomega0), sigma(nsigma), e(ne)
{
  localName = nname;
}

CLOscillatingPlate::CLOscillatingPlate(const XMLNode& XML, DYNAMO::SimData* tmp):
  CLocal(tmp, "LocalWall")
{
  operator<<(XML);
}

CLocalEvent 
CLOscillatingPlate::getEvent(const CParticle& part) const
{
#ifdef ISSS_DEBUG
  if (!Sim->Dynamics.Liouvillean().isUpToDate(part))
    D_throw() << "Particle is not up to date";
#endif

  Iflt rvec = part.getPosition()[0] - x0;
  Iflt maxlength = sigma + xi;
  Iflt tmax = 0;

  //Simple test to see if they are approaching
  if (((rvec > maxlength) && (part.getVelocity()[0] > 0))
      || ((rvec < -maxlength) && (part.getVelocity()[0] < 0)))
      return CLocalEvent(part, HUGE_VAL, NONE, *this);
  
  if (part.getVelocity()[0] > 0)
    tmax = (maxlength - rvec) / part.getVelocity()[0];
  else
    tmax = (-maxlength - rvec) / part.getVelocity()[0];

  
  
  I_cout() << "True";

  return CLocalEvent(part, HUGE_VAL, NONE, *this);
}



void
CLOscillatingPlate::runEvent(const CParticle& part, const CLocalEvent& iEvent) const
{
//  ++Sim->lNColl;
//
//  //Run the collision and catch the data
//  CNParticleData EDat(Sim->Dynamics.Liouvillean().runWallCollision
//		      (part, vNorm, e));
//
//  Sim->signalParticleUpdate(EDat);
//
//  //Now we're past the event update the scheduler and plugins
//  Sim->ptrScheduler->fullUpdate(part);
//  
//  BOOST_FOREACH(smrtPlugPtr<COutputPlugin> & Ptr, Sim->outputPlugins)
//    Ptr->eventUpdate(iEvent, EDat);
}

bool 
CLOscillatingPlate::isInCell(const Vector & Origin, const Vector& CellDim) const
{
  return true;
}

void 
CLOscillatingPlate::initialise(size_t nID)
{
  ID = nID;
}

void 
CLOscillatingPlate::operator<<(const XMLNode& XML)
{
  range.set_ptr(CRange::loadClass(XML,Sim));
  
  try {
    e = boost::lexical_cast<Iflt>(XML.getAttribute("Elasticity"));

    x0 =  boost::lexical_cast<Iflt>(XML.getAttribute("X0"))
      * Sim->Dynamics.units().unitLength();

    xi =  boost::lexical_cast<Iflt>(XML.getAttribute("Xi"))
      * Sim->Dynamics.units().unitLength();

    omega0 =  boost::lexical_cast<Iflt>(XML.getAttribute("Omega0"))
      / Sim->Dynamics.units().unitTime();

    sigma = boost::lexical_cast<Iflt>(XML.getAttribute("Sigma"))
      * Sim->Dynamics.units().unitLength();

    localName = XML.getAttribute("Name");
  } 
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in CLOscillatingPlate";
    }
}

void 
CLOscillatingPlate::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "OscillatingPlate" 
      << xmlw::attr("Name") << localName
      << xmlw::attr("Elasticity") << e
      << xmlw::attr("X0") << x0 / Sim->Dynamics.units().unitLength()
      << xmlw::attr("Xi") << xi / Sim->Dynamics.units().unitLength()
      << xmlw::attr("Omega0") << omega0 * Sim->Dynamics.units().unitTime()
      << xmlw::attr("Sigma") << sigma / Sim->Dynamics.units().unitLength();

}

double 
CLOscillatingPlate::getPosition() const
{
  return xi * std::cos(omega0 * Sim->dSysTime) + x0;
}

void 
CLOscillatingPlate::write_povray_info(std::ostream& os) const
{
  os << "object {\n union { plane {\n  <" << -1 << ", " << 0
     << ", " << 0 << ">, 0 texture{pigment { color rgb<0.5,0.5,0.5>}}\n translate <" << -sigma 
     << ",0,0> } \n plane {\n  <" << 1 << ", " << 0
     << ", " << 0 << ">, 0 texture{pigment { color rgb<0.5,0.5,0.5>}}\n translate <" << +sigma 
     << ",0,0> } }\n clipped_by{box {\n  <" << -Sim->aspectRatio[0]/2 
     << ", " << -Sim->aspectRatio[1]/2 << ", " << -Sim->aspectRatio[2]/2 
     << ">, <" << Sim->aspectRatio[0]/2 << ", " << Sim->aspectRatio[1]/2 
     << ", " << Sim->aspectRatio[2]/2 << "> }\n}\n translate <" << getPosition()<< 
    ","<< 0 << "," << 0 << ">\n}\n";
}
