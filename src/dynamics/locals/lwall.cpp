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

#include "lwall.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "localEvent.hpp"
#include "../NparticleEventData.hpp"
#include "../overlapFunc/CubePlane.hpp"
#include "../units/units.hpp"
#include "../../datatypes/vector.xml.hpp"
#include "../../schedulers/scheduler.hpp"


CLWall::CLWall(DYNAMO::SimData* nSim, Iflt ne, CVector<> nnorm, 
	       CVector<> norigin, std::string nname, CRange* nRange):
  CLocal(nRange, nSim, "LocalWall"),
  vNorm(nnorm),
  vPosition(norigin),
  e(ne)
{
  localName = nname;
}

CLWall::CLWall(const XMLNode& XML, DYNAMO::SimData* tmp):
  CLocal(tmp, "LocalWall")
{
  operator<<(XML);
}

CLocalEvent 
CLWall::getEvent(const CParticle& part) const
{
#ifdef ISSS_DEBUG
  if (!Sim->Dynamics.Liouvillean().isUpToDate(part))
    D_throw() << "Particle is not up to date";
#endif

  return CLocalEvent(part, Sim->Dynamics.Liouvillean().getWallCollision
		     (part, vPosition, vNorm), WALL, *this);
}

void
CLWall::runEvent(const CParticle& part) const
{
  ++Sim->lNColl;

  Sim->Dynamics.Liouvillean().updateParticle(part);
  CLocalEvent iEvent(getEvent(part));
  
  if (iEvent.getType() == NONE)
    D_throw() << "No global collision found\n"
	      << iEvent.stringData(Sim);
  
#ifdef DYNAMO_DEBUG 
  if (isnan(iEvent.getdt()))
    D_throw() << "A NAN Global collision time has been found\n"
	      << iEvent.stringData(Sim);
  
  if (iEvent.getdt() == HUGE_VAL)
    D_throw() << "An infinite (not marked as NONE) Global collision time has been found\n"
	      << iEvent.stringData(Sim);
#endif
  
  Sim->dSysTime += iEvent.getdt();
  
  Sim->ptrScheduler->stream(iEvent.getdt());
  
  Sim->Dynamics.stream(iEvent.getdt());
    
  //Run the collision and catch the data
  CNParticleData EDat(Sim->Dynamics.Liouvillean().runWallCollision
		      (part, vNorm, e));

  //Now we're past the event update the scheduler and plugins
  Sim->ptrScheduler->fullUpdate(part);
  
  BOOST_FOREACH(smrtPlugPtr<COutputPlugin> & Ptr, Sim->outputPlugins)
    Ptr->eventUpdate(iEvent, EDat);

  Sim->signalParticleUpdate(EDat.L1partChanges.front());
}

bool 
CLWall::isInCell(const CVector<>& Origin, const CVector<>& CellDim) const
{
  return DYNAMO::OverlapFunctions::CubePlane
    (Origin, CellDim, vPosition, vNorm);
}

void 
CLWall::initialise(size_t nID)
{
  ID = nID;
}

void 
CLWall::operator<<(const XMLNode& XML)
{
  range.set_ptr(CRange::loadClass(XML,Sim));
  
  try {
    e = boost::lexical_cast<Iflt>(XML.getAttribute("Elasticity"));
    XMLNode xBrowseNode = XML.getChildNode("Norm");
    localName = XML.getAttribute("Name");
    vNorm << xBrowseNode;
    vNorm = vNorm.unitVector();
    xBrowseNode = XML.getChildNode("Origin");
    vPosition << xBrowseNode;
    vPosition *= Sim->Dynamics.units().unitLength();
  } 
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in CLWall";
    }
}

void 
CLWall::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "Wall" 
      << xmlw::attr("Name") << localName
      << xmlw::attr("Elasticity") << e
      << range
      << xmlw::tag("Norm")
      << vNorm
      << xmlw::endtag("Norm")
      << xmlw::tag("Origin")
      << vPosition / Sim->Dynamics.units().unitLength()
      << xmlw::endtag("Origin");
}
