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

#include "gListAndCell.hpp"
#include "globEvent.hpp"
#include "../NparticleEventData.hpp"
#include "../../extcode/xmlParser.h"
#include "../../extcode/xmlwriter.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../units/units.hpp"
#include "../ranges/1RAll.hpp"
#include "../ranges/2RSingle.hpp"
#include "../../schedulers/scheduler.hpp"
#include "../locals/local.hpp"

CGListAndCell::CGListAndCell(const DYNAMO::SimData* nSim, 
			     const std::string& name):
  CGCells(nSim, "ListAndCellNBList", NULL),
  largestParticles(NULL)
{
  globName = name;
  I_cout() << "Cells Loaded";
}

CGListAndCell::CGListAndCell(const XMLNode &XML, const DYNAMO::SimData* ptrSim):
  CGCells(ptrSim, "ListAndCellNBList"),
  largestParticles(NULL)
{
  operator<<(XML);

  I_cout() << "Cells Loaded";
}

void 
CGListAndCell::operator<<(const XMLNode& XML)
{
  try {
    if (XML.isAttributeSet("lambda"))
      lambda = boost::lexical_cast<Iflt>
	(XML.getAttribute("lambda"));
    
    globName = XML.getAttribute("Name");	
  }
  catch(...)
    {
      D_throw() << "Error loading CGListAndCell";
    }
  
  if (lambda < 0.0 || lambda > 1.0)
    D_throw() << "Lambda out of bounds [0,1), lambda = " << lambda;
}

void 
CGListAndCell::initialise(size_t nID)
{
  ID=nID;

  Iflt secondMaxDiam = 0.0;
  Iflt maxdiam = 0.0;
  const CInteraction* biggest = NULL;
  
  if (Sim->Dynamics.getInteractions().size() < 2)
    D_throw() << "This scheduler doesn't work unless you have more than 1 interaction";

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
    D_throw() << "For the MultListSpecial scheduler to work, the largest interaction"
	      << " must use C2RSingle to adapt a CRange to a C2Range";

  //We have the biggest interaction, steal its range pointer
  largestParticles.set_ptr(static_cast<const C2RSingle*>
			   (biggest->getRange().get_ptr())->getRange()->Clone());
  
  I_cout() << "Found that interaction \"" << biggest->getName() <<  "\" had the longest interaction range"
	   << "\nUsing its range to increase number of cells";
  
  reinitialise(secondMaxDiam);
}

void
CGListAndCell::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "ListAndCell"
      << xmlw::attr("Lambda") << lambda
      << xmlw::attr("Name") << globName;
}

void 
CGListAndCell::getParticleNeighbourhood(const CParticle& part,
				  const nbhoodFunc& func) const
{
  CGCells::getParticleNeighbourhood(part, func);
  
  if (largestParticles->isInRange(part))
    //This is a large particle so tell it to compare against all other large particles
    BOOST_FOREACH(const size_t& ID, *largestParticles)
      if (ID != part.getID())
	func(part, ID);
}
