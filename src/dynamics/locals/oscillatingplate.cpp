/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
				       Vector nrw0, Vector nnhat,
				       Iflt nomega0, Iflt nsigma, Iflt ne,
				       Iflt ndelta, Iflt nmass, std::string nname, 
				       CRange* nRange, Iflt timeshift):
  CLocal(nRange, nSim, "OscillatingPlate"),
  rw0(nrw0), nhat(nnhat), omega0(nomega0), sigma(nsigma), 
  e(ne), delta(ndelta), mass(nmass), timeshift(0), lastID(-1), lastdSysTime(HUGE_VAL)
{
  localName = nname;
}

CLOscillatingPlate::CLOscillatingPlate(const XMLNode& XML, DYNAMO::SimData* tmp):
  CLocal(tmp, "OscillatingPlate"),
  lastID(-1), lastdSysTime(HUGE_VAL)
{
  operator<<(XML);
}

CLocalEvent 
CLOscillatingPlate::getEvent(const CParticle& part) const
{
#ifdef ISSS_DEBUG
  if (!Sim->dynamics.getLiouvillean().isUpToDate(part))
    D_throw() << "Particle is not up to date";
#endif

  bool caution = ((part.getID() == lastID) && (lastdSysTime == Sim->dSysTime));

  Iflt reducedt = Sim->dSysTime 
    - 2.0 * M_PIl * int(Sim->dSysTime * omega0 / (2.0*M_PIl)) / omega0;

  Iflt dt = Sim->dynamics.getLiouvillean().getPointPlateCollision
    (part, rw0, nhat, delta, omega0, sigma, reducedt + timeshift, 
     caution);

  if (dt != HUGE_VAL)
    return CLocalEvent(part, dt, WALL, *this);
  else
    return CLocalEvent(part, HUGE_VAL, NONE, *this);
}



void
CLOscillatingPlate::runEvent(const CParticle& part, const CLocalEvent& iEvent) const
{
  ++Sim->lNColl;
  
  //Run the collision and catch the data
  CNParticleData EDat(Sim->dynamics.getLiouvillean().runOscilatingPlate
		      (part, rw0, nhat, delta, omega0, sigma, mass, 
		       e, timeshift));

  lastdSysTime = Sim->dSysTime;
  lastID = part.getID();

  Sim->signalParticleUpdate(EDat);

  BOOST_FOREACH(smrtPlugPtr<OutputPlugin> & Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(iEvent, EDat);

  //Now we're past the event update the scheduler and plugins
  //Sim->ptrScheduler->fullUpdate(part);
  Sim->ptrScheduler->rebuildList();
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

    XMLNode xBrowseNode = XML.getChildNode("Norm");
    nhat << xBrowseNode;
    nhat /= nhat.nrm();

    xBrowseNode = XML.getChildNode("Origin");
    rw0 << xBrowseNode;
    rw0 *= Sim->dynamics.units().unitLength();

    omega0 =  boost::lexical_cast<Iflt>(XML.getAttribute("Omega0"))
      / Sim->dynamics.units().unitTime();

    sigma = boost::lexical_cast<Iflt>(XML.getAttribute("Sigma"))
      * Sim->dynamics.units().unitLength();

    delta = boost::lexical_cast<Iflt>(XML.getAttribute("Delta"))
      * Sim->dynamics.units().unitLength();

    mass = boost::lexical_cast<Iflt>(XML.getAttribute("Mass"))
      * Sim->dynamics.units().unitMass();

    timeshift = boost::lexical_cast<Iflt>(XML.getAttribute("TimeShift"))
      * Sim->dynamics.units().unitTime();

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
  Iflt tmp = Sim->dSysTime + timeshift;

  tmp -= 2.0 * PI * int(tmp * omega0 / (2.0 * PI) ) / omega0;

  XML << xmlw::attr("Type") << "OscillatingPlate" 
      << xmlw::attr("Name") << localName
      << xmlw::attr("Elasticity") << e
      << xmlw::attr("Omega0") << omega0 * Sim->dynamics.units().unitTime()
      << xmlw::attr("Sigma") << sigma / Sim->dynamics.units().unitLength()
      << xmlw::attr("Delta") << delta / Sim->dynamics.units().unitLength()
      << xmlw::attr("Mass") << mass / Sim->dynamics.units().unitMass()
      << xmlw::attr("TimeShift") << tmp / Sim->dynamics.units().unitTime()
      << range
      << xmlw::tag("Norm")
      << nhat
      << xmlw::endtag("Norm")
      << xmlw::tag("Origin")
      << rw0 / Sim->dynamics.units().unitLength()
      << xmlw::endtag("Origin");

}

Vector
CLOscillatingPlate::getPosition() const
{
  return nhat * (delta * std::cos(omega0 * (Sim->dSysTime + timeshift))) + rw0;
}

Vector
CLOscillatingPlate::getVelocity() const
{
  return - nhat * (delta * omega0 * std::sin(omega0 * (Sim->dSysTime + timeshift)));
}

Iflt 
CLOscillatingPlate::getPlateEnergy() const
{
  return 0.5 * mass 
    * (std::pow(omega0 * delta * std::cos(omega0 * (Sim->dSysTime + timeshift)), 2)
       +  std::pow(omega0 * delta * std::sin(omega0 * (Sim->dSysTime + timeshift)), 2));
}

void 
CLOscillatingPlate::write_povray_info(std::ostream& os) const
{
  Vector pos = getPosition();
  //The walls are \pm0.25 thick and are set 0.5 away from the COM surface
  //to give the appearance of proper walls, not COM walls
  Vector WallLoc1 = pos + nhat * (sigma + 0.75 * Sim->dynamics.units().unitLength());
  Vector WallLoc2 = pos - nhat * (sigma + 0.75 * Sim->dynamics.units().unitLength());

  Sim->dynamics.BCs().applyBC(WallLoc1);
  Sim->dynamics.BCs().applyBC(WallLoc2);

  os << "object { intersection { object { box { <-0.5, " << -0.25 * Sim->dynamics.units().unitLength() 
     << ", -0.5>, <0.5, " << +0.25 * Sim->dynamics.units().unitLength() 
     << ", 0.5> } Point_At_Trans(<"
     << nhat[0] << "," << nhat[1] << "," << nhat[2] << ">) translate <"
     <<  WallLoc1[0] << "," <<  WallLoc1[1] << "," <<  WallLoc1[2] 
     << "> }"

     << "\nbox { <" 
     << -Sim->aspectRatio[0]/2 - Sim->dynamics.units().unitLength() 
     << "," << -Sim->aspectRatio[1]/2 - Sim->dynamics.units().unitLength()  
     << "," << -Sim->aspectRatio[2]/2 - Sim->dynamics.units().unitLength() 
     << ">,"
     << "<" << Sim->aspectRatio[0]/2 + Sim->dynamics.units().unitLength()
     << "," << Sim->aspectRatio[1]/2 + Sim->dynamics.units().unitLength()
     << "," << Sim->aspectRatio[2]/2 + Sim->dynamics.units().unitLength()
     << "> }\n"
     << "} pigment { Col_Glass_Bluish } finish { F_Glass5 } }\n";

  os << "object { intersection { object { box { <-0.5, " << -0.25 * Sim->dynamics.units().unitLength()
     << ", -0.5>, <0.5, " << 0.25 * Sim->dynamics.units().unitLength() 
     << ", 0.5> } Point_At_Trans(<"
     << -nhat[0] << "," << -nhat[1] << "," << -nhat[2] << ">) translate <"
     <<  WallLoc2[0] << "," <<  WallLoc2[1] << "," <<  WallLoc2[2] 
     << "> }"

     << "\nbox { <" 
     << -Sim->aspectRatio[0]/2 - Sim->dynamics.units().unitLength() 
     << "," << -Sim->aspectRatio[1]/2 - Sim->dynamics.units().unitLength()  
     << "," << -Sim->aspectRatio[2]/2 - Sim->dynamics.units().unitLength() 
     << ">,"
     << "<" << Sim->aspectRatio[0]/2 + Sim->dynamics.units().unitLength()
     << "," << Sim->aspectRatio[1]/2 + Sim->dynamics.units().unitLength()
     << "," << Sim->aspectRatio[2]/2 + Sim->dynamics.units().unitLength()
     << "> }\n"
     << "} pigment { Col_Glass_Bluish } finish { F_Glass5 } }\n";
}
