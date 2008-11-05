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
#include "../../schedulers/scheduler.hpp"
#include "../locals/local.hpp"

CGListAndCell::CGListAndCell(const DYNAMO::SimData* nSim, 
			     const std::string& name):
  CGCells(nSim, "ListAndCellNBList", NULL)
{
  globName = name;
  I_cout() << "Cells Loaded";
}

CGListAndCell::CGListAndCell(const XMLNode &XML, const DYNAMO::SimData* ptrSim):
  CGCells(ptrSim, "GlobalCellularEvent")
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

CGlobEvent 
CGListAndCell::getEvent(const CParticle& part) const
{
  //Update the wall collision

  return CGlobEvent(part, Sim->Dynamics.Liouvillean().
		    getSquareCellCollision2
		    (part, cells[partCellData[part.getID()].cell].origin, 
		     cellDimension), VIRTUAL, *this);
}

CNParticleData
CGListAndCell::runEvent(const CGlobEvent& event) const
{
}

void 
CGListAndCell::initialise(size_t nID)
{
  ID=nID;
  reinitialise(Sim->Dynamics.getLongestInteraction());
}

void
CGListAndCell::reinitialise(const Iflt& maxdiam)
{
  I_cout() << "Reinitialising on collision " << Sim->lNColl;

  //Create the cells
  addCells(maxdiam, false);

  addLocalEvents();

  ReInitNotify();
}

void
CGListAndCell::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "ListAndCells"
      << xmlw::attr("Name") << globName;
}

void 
CGListAndCell::getParticleNeighbourhood(const CParticle& part,
				  const nbhoodFunc& func) const
{
  BOOST_FOREACH(const int& nb, cells[partCellData[part.getID()].cell].neighbours)
    for (int next = cells[nb].list;
	 next != -1; next = partCellData[next].next)
      if (next != static_cast<int>(part.getID()))
	func(part, next);
}

void 
CGListAndCell::getParticleLocalNeighbourhood(const CParticle& part, 
				       const nbhoodFunc& func) const
{
  BOOST_FOREACH(const size_t& id, cells[partCellData[part.getID()].cell].locals)
    func(part, id);
}
