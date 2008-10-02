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

#include "multlist_special.hpp"
#include "../dynamics/interactions/intEvent.hpp"
#include "../extcode/threadpool.hpp"
#include "../simulation/particle.hpp"
#include "../dynamics/dynamics.hpp"
#include "../dynamics/liouvillean/liouvillean.hpp"
#include "../dynamics/BC/BC.hpp"
#include "../base/is_simdata.hpp"
#include "../base/is_base.hpp"
#include "../dynamics/globals/globEvent.hpp"
#include "../dynamics/systems/system.hpp"
#include <cmath> //for huge val
#include <boost/lexical_cast.hpp>
#include "../extcode/xmlParser.h"
#include "../dynamics/globals/global.hpp"
#include "../dynamics/globals/globEvent.hpp"
#include "../dynamics/ranges/1range.hpp"
#include "../dynamics/ranges/2RSingle.hpp"

CSMultListSpecial::CSMultListSpecial(const XMLNode& XML, const DYNAMO::SimData* Sim):
  CSMultList(Sim, "MultListSpecial"),
  specialParticles(0,0)
{
  I_cout() << "Multi List Cellular Algorithm with special particles";
  operator<<(XML);
}

CSMultListSpecial::CSMultListSpecial(const DYNAMO::SimData* Sim):
  CSMultList(Sim, "MultListSpecial"),
  specialParticles(0,0)
{ I_cout() << "Multi List Cellular Algorithmn with special particles"; }

void 
CSMultListSpecial::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "MultListSpecial";
  CSCells::outputXML(XML);
}

void 
CSMultListSpecial::operator<<(const XMLNode& XML)
{
  CSCells::operator<<(XML);
}

void
CSMultListSpecial::initialise()
{ 
  Iflt secondMaxDiam = 0.0;
  Iflt maxdiam = 0.0;
  const CInteraction* biggest = NULL;
  
  if (Sim->Dynamics.getInteractions().size() < 2)
    I_throw() << "This scheduler doesn't work unless you have more than 1 interaction";

  //Find the largest interaction
  BOOST_FOREACH(const smrtPlugPtr<CInteraction>& intPtr, Sim->Dynamics.getInteractions())
    if (intPtr->maxIntDist() > maxdiam)
      {
	maxdiam = intPtr->maxIntDist();
	biggest = intPtr.get_ptr();
      }

  //Find the second biggest
  BOOST_FOREACH(const smrtPlugPtr<CInteraction>& intPtr, Sim->Dynamics.getInteractions())
    if (intPtr->maxIntDist() > secondMaxDiam)
      if (intPtr.get_ptr() != biggest)
	secondMaxDiam = intPtr->maxIntDist();

  if (dynamic_cast<const C2RSingle*>(biggest->getRange().get_ptr()) == NULL)
    I_throw() << "For the MultListSpecial scheduler to work, the largest interaction must use 2Single to adapt a 1range to a 2range";

  if (static_cast<const CRRange*>(static_cast<const C2RSingle*>(biggest->getRange().get_ptr())->getRange().get_ptr()) == NULL)
    I_throw() << "I'm being a pain I know but if the largest interaction was using 2Single and a 1Range it would be quicker";

  //We have the biggest interaction, steal its range pointer
  specialParticles = *(static_cast<const CRRange*>(static_cast<const C2RSingle*>(biggest->getRange().get_ptr())->getRange().get_ptr()));

  I_cout() << "Found that interaction \"" << biggest->getName() <<  "\" had the longest interaction range"
	   << "\nUsing its range to speed calculations";

  reinitialise(secondMaxDiam);
}

void 
CSMultListSpecial::update(const CParticle& part)
{
  //Invalidate previous entries
  ++eventCount[part.getID()];
  eventHeap[part.getID()].clear();
  addNewEvents(part);
  addSpecialEvents(part);
  eventHeap.update(part.getID());
}

void
CSMultListSpecial::reinitialise(Iflt maxdiam)
{
  I_cout() << "Reinitialising on collision " << Sim->lNColl;
  eventHeap.clear();
  eventHeap.resize(Sim->lN);
  eventCount.clear();
  eventCount.resize(Sim->lN, 0);

  //Create the cells
  addCells(maxdiam, false);
  
  //Now initialise the interactions
  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    {
      addNewEvents_init(part);
      addSpecialEvents_init(part);
    }
  
  eventHeap.init();

#ifndef CBT
  I_cout() << "BPQ: Number of lists " << eventHeap.NLists();
  I_cout() << "BPQ: Scale Factor " << eventHeap.scaleFactor();
#endif
}

void
CSMultListSpecial::addSpecialEvents(const CParticle& part)
{
  if (specialParticles.isInRange(part))
    BOOST_FOREACH(const unsigned long& ID, specialParticles)
      if (part.getID() != ID)
	eventHeap.push(intPart(Sim->Dynamics.getEvent(part, Sim->vParticleList[ID]), 
			       eventCount[ID]), part.getID());
}

void
CSMultListSpecial::addSpecialEvents_init(const CParticle& part)
{
  if (specialParticles.isInRange(part))
    BOOST_FOREACH(const unsigned long& ID, specialParticles)
      if (part.getID() < ID)
	eventHeap.push(intPart(Sim->Dynamics.getEvent(part, Sim->vParticleList[ID]), 
			       eventCount[ID]), part.getID());
}
