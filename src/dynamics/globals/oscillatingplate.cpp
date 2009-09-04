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
#include "globEvent.hpp"
#include "../NparticleEventData.hpp"
#include "../../extcode/xmlParser.h"
#include "../../extcode/xmlwriter.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../units/units.hpp"
#include "../ranges/1RAll.hpp"
#include "../../schedulers/scheduler.hpp"
#include "../locals/local.hpp"
#include "../BC/LEBC.hpp"
#include <boost/static_assert.hpp>
#include <boost/math/special_functions/pow.hpp>


CGOscillatingPlate::CGOscillatingPlate(DYNAMO::SimData* nSim, Iflt nx0, Iflt nxi, 
				       Iflt nomega0, Iflt nsigma, const std::string& name):
  CGlobal(nSim, "OscillatingPlate"),
  x0(nx0), xi(nxi), omega0(nomega0), sigma(nsigma)
{
  globName = name;
  I_cout() << "Oscillating Plate loaded";
}

CGOscillatingPlate::CGOscillatingPlate(const XMLNode &XML, DYNAMO::SimData* ptrSim):
  CGlobal(ptrSim, "OscillatingPlate")
{
  operator<<(XML);

  I_cout() << "Single occupancy cells loaded";
}

void 
CGOscillatingPlate::operator<<(const XMLNode& XML)
{
  try {
    globName = XML.getAttribute("Name");	
    x0 =  boost::lexical_cast<Iflt>(XML.getAttribute("X0"));
    xi =  boost::lexical_cast<Iflt>(XML.getAttribute("Xi"));
    omega0 =  boost::lexical_cast<Iflt>(XML.getAttribute("Omega0"));
    sigma = boost::lexical_cast<Iflt>(XML.getAttribute("Sigma"));
  }
  catch(...)
    {
      D_throw() << "Error loading CGOscillatingPlate";
    }
}

CGlobEvent 
CGOscillatingPlate::getEvent(const CParticle& part) const
{
//#ifdef ISSS_DEBUG
//  if (!Sim->Dynamics.Liouvillean().isUpToDate(part))
//    D_throw() << "Particle is not up to date";
//#endif
//
//  //This 
//  //Sim->Dynamics.Liouvillean().updateParticle(part);
//  //is not required as we compensate for the delay using 
//  //Sim->Dynamics.Liouvillean().getParticleDelay(part)
//
//  Vector CellOrigin;
//  size_t ID(part.getID());
//
//  for (size_t iDim(0); iDim < NDIM; ++iDim)
//    {
//      CellOrigin[iDim] = (ID % cuberootN) * cellDimension[iDim] - 0.5*Sim->aspectRatio[iDim];
//      ID /= cuberootN;
//    }
//
//  return CGlobEvent(part,
//		    Sim->Dynamics.Liouvillean().
//		    getSquareCellCollision2
//		    (part, CellOrigin,
//		     cellDimension)
//		    -Sim->Dynamics.Liouvillean().getParticleDelay(part),
//		    CELL, *this);
  D_throw() << "Not implemented";
}

void
CGOscillatingPlate::runEvent(const CParticle& part) const
{
//  Sim->Dynamics.Liouvillean().updateParticle(part);
//
//  Vector CellOrigin;
//  size_t ID(part.getID());
//
//  for (size_t iDim(0); iDim < NDIM; ++iDim)
//    {
//      CellOrigin[iDim] = (ID % cuberootN) * cellDimension[iDim] - 0.5*Sim->aspectRatio[iDim];
//      ID /= cuberootN;
//    }
//  
//  //Determine the cell transition direction, its saved
//  size_t cellDirection(Sim->Dynamics.Liouvillean().
//		       getSquareCellCollision3
//		       (part, CellOrigin, 
//			cellDimension));
//
//  CGlobEvent iEvent(getEvent(part));
//
//#ifdef DYNAMO_DEBUG 
//  if (isnan(iEvent.getdt()))
//    D_throw() << "A NAN Interaction collision time has been found"
//	      << iEvent.stringData(Sim);
//  
//  if (iEvent.getdt() == HUGE_VAL)
//    D_throw() << "An infinite Interaction (not marked as NONE) collision time has been found\n"
//	      << iEvent.stringData(Sim);
//#endif
//
//  Sim->dSysTime += iEvent.getdt();
//    
//  Sim->ptrScheduler->stream(iEvent.getdt());
//  
//  Sim->Dynamics.stream(iEvent.getdt());
//
//  Vector vNorm(0,0,0);
//
//  Vector pos(part.getPosition()), vel(part.getVelocity());
//
//  Sim->Dynamics.BCs().setPBC(pos, vel);
//
//  vNorm[cellDirection] = (vel[cellDirection] > 0) ? -1 : +1; 
//    
//  //Run the collision and catch the data
//  CNParticleData EDat(Sim->Dynamics.Liouvillean().runWallCollision
//		      (part, vNorm, 1.0));
//
//  Sim->signalParticleUpdate(EDat);
//
//  //Now we're past the event update the scheduler and plugins
//  Sim->ptrScheduler->fullUpdate(part);
//  
//  BOOST_FOREACH(smrtPlugPtr<COutputPlugin> & Ptr, Sim->outputPlugins)
//    Ptr->eventUpdate(iEvent, EDat);
//
}

void 
CGOscillatingPlate::initialise(size_t nID)
{
  ID=nID;
}

void
CGOscillatingPlate::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "OscillatingPlate"
      << xmlw::attr("Name") << globName
      << xmlw::attr("X0") << x0
      << xmlw::attr("Xi") << xi
      << xmlw::attr("Omega0") << omega0
      << xmlw::attr("Sigma") << sigma;
}
