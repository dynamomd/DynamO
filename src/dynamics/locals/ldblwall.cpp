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

#include "ldblwall.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "localEvent.hpp"
#include "../NparticleEventData.hpp"
#include "../overlapFunc/CubePlane.hpp"
#include "../units/units.hpp"
#include "../../datatypes/vector.xml.hpp"
#include "../../schedulers/scheduler.hpp"


CLDblWall::CLDblWall(DYNAMO::SimData* nSim, Iflt ne, Vector  nnorm, 
	       Vector  norigin, std::string nname, CRange* nRange):
  Local(nRange, nSim, "LocalWall"),
  vNorm(nnorm),
  vPosition(norigin),
  e(ne)
{
  localName = nname;
}

CLDblWall::CLDblWall(const XMLNode& XML, DYNAMO::SimData* tmp):
  Local(tmp, "LocalDoubleWall")
{
  operator<<(XML);
}

LocalEvent 
CLDblWall::getEvent(const Particle& part) const
{
#ifdef ISSS_DEBUG
  if (!Sim->dynamics.getLiouvillean().isUpToDate(part))
    M_throw() << "Particle is not up to date";
#endif

  if (part.getID() == lastID) return LocalEvent(part, HUGE_VAL, NONE, *this);
  
  Vector rij = part.getPosition() - vPosition;
  Sim->dynamics.BCs().applyBC(rij);

  Vector norm(vNorm);
  if ((norm | rij) < 0)
    norm *= -1;

  return LocalEvent(part, Sim->dynamics.getLiouvillean().getWallCollision
		     (part, vPosition, norm), WALL, *this);
}

void
CLDblWall::runEvent(const Particle& part, const LocalEvent& iEvent) const
{
  ++Sim->eventCount;

  Vector norm = vNorm;

  Vector rij = part.getPosition() - vPosition;
  Sim->dynamics.BCs().applyBC(rij);
  
  if ((norm | rij) < 0)
    norm *= -1;

  //Run the collision and catch the data
  NEventData EDat(Sim->dynamics.getLiouvillean().runWallCollision
		      (part, norm, e));

  Sim->signalParticleUpdate(EDat);

  //Must do this after the signal is run
  lastID = part.getID();

  //Now we're past the event update the scheduler and plugins
  Sim->ptrScheduler->fullUpdate(part);
  
  BOOST_FOREACH(ClonePtr<OutputPlugin> & Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(iEvent, EDat);
}

bool 
CLDblWall::isInCell(const Vector & Origin, const Vector& CellDim) const
{
  return DYNAMO::OverlapFunctions::CubePlane
    (Origin, CellDim, vPosition, vNorm);
}

void 
CLDblWall::initialise(size_t nID)
{
  ID = nID;
  lastID = std::numeric_limits<size_t>::max();

  Sim->registerParticleUpdateFunc
    (fastdelegate::MakeDelegate(this, &CLDblWall::particleUpdate));

}

void
CLDblWall::particleUpdate(const NEventData& PDat) const
{
  if (lastID == std::numeric_limits<size_t>::max()) return;

  BOOST_FOREACH(const ParticleEventData& pdat, PDat.L1partChanges)
    if (pdat.getParticle().getID() == lastID)
      {
	lastID = std::numeric_limits<size_t>::max();
	return;
      }
  
  BOOST_FOREACH(const PairEventData& pdat, PDat.L2partChanges)
    if (pdat.particle1_.getParticle().getID() == lastID 
	|| pdat.particle2_.getParticle().getID() == lastID)
      {
	lastID = std::numeric_limits<size_t>::max();
	return;
      }
}

void 
CLDblWall::operator<<(const XMLNode& XML)
{
  range.set_ptr(CRange::loadClass(XML,Sim));
  
  try {
    e = boost::lexical_cast<Iflt>(XML.getAttribute("Elasticity"));
    XMLNode xBrowseNode = XML.getChildNode("Norm");
    localName = XML.getAttribute("Name");
    vNorm << xBrowseNode;
    vNorm /= vNorm.nrm();
    xBrowseNode = XML.getChildNode("Origin");
    vPosition << xBrowseNode;
    vPosition *= Sim->dynamics.units().unitLength();
  } 
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in CLDblWall";
    }
}

void 
CLDblWall::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Type") << "DoubleWall" 
      << xml::attr("Name") << localName
      << xml::attr("Elasticity") << e
      << range
      << xml::tag("Norm")
      << vNorm
      << xml::endtag("Norm")
      << xml::tag("Origin")
      << vPosition / Sim->dynamics.units().unitLength()
      << xml::endtag("Origin");
}

void 
CLDblWall::write_povray_info(std::ostream& os) const
{
  os << "object {\n plane {\n  <" << vNorm[0] << ", " << vNorm[1] 
     << ", " << vNorm[2] << ">, 0 texture{pigment { color rgb<0.5,0.5,0.5>}}}\n clipped_by{box {\n  <" << -Sim->aspectRatio[0]/2 
     << ", " << -Sim->aspectRatio[1]/2 << ", " << -Sim->aspectRatio[2]/2 
     << ">, <" << Sim->aspectRatio[0]/2 << ", " << Sim->aspectRatio[1]/2 
     << ", " << Sim->aspectRatio[2]/2 << "> }\n}\n translate <" << vPosition[0] << 
    ","<< vPosition[1] << "," << vPosition[2] << ">\n}\n";
}
