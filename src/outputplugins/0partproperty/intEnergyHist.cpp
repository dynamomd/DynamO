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

#include "intEnergyHist.hpp"
#include <fstream>
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_exception.hpp"
#include "../../base/is_simdata.hpp"
#include "../1partproperty/uenergy.hpp"

COPIntEnergyHist::COPIntEnergyHist(const DYNAMO::SimData* tmp, const XMLNode& XML):
  COPCollTicker(tmp,"InternalEnergyHistogram", 10),//Before COPEnergy
  intEnergyHist(1.0),
  ptrCOPEnergy(NULL),
  weight(0.0),
  binwidth(1.0)
{
  operator<<(XML);
}

void 
COPIntEnergyHist::operator<<(const XMLNode& XML)
{  
  try 
    {
      if (XML.isAttributeSet("BinWidth"))
	binwidth = boost::lexical_cast<Iflt>
	   (XML.getAttribute("BinWidth"));
    }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in COPIntEnergyHist";
    }  
}

void 
COPIntEnergyHist::initialise() 
{
  ptrCOPEnergy = Sim->getOutputPlugin<COPUEnergy>();
  intEnergyHist = C1DWeightHistogram(binwidth * Sim->dynamics.units().unitEnergy());
}

void 
COPIntEnergyHist::changeSystem(COutputPlugin* EHist2)
{
  //Add the current data
  intEnergyHist.addVal(ptrCOPEnergy->getSimU(), weight);
  //Same for the other histogram
  static_cast<COPIntEnergyHist*>(EHist2)->intEnergyHist.addVal
    (static_cast<COPIntEnergyHist*>(EHist2)->ptrCOPEnergy->getSimU(), 
     static_cast<COPIntEnergyHist*>(EHist2)->weight);

  //Now swap over the data
  std::swap(Sim, static_cast<COPIntEnergyHist*>(EHist2)->Sim);

  //NEVER SWAP THE PLUGIN POINTERS! they don't change
  //std::swap(ptrCOPEnergy, static_cast<COPIntEnergyHist*>(EHist2)->ptrCOPEnergy);

  //Reset the weighting
  weight = 0.0;
  static_cast<COPIntEnergyHist*>(EHist2)->weight = 0.0;
}

void 
COPIntEnergyHist::stream(Iflt dt)
{
  weight += dt;
}

void 
COPIntEnergyHist::ticker()
{
  intEnergyHist.addVal(ptrCOPEnergy->getSimU(), weight);
  weight = 0.0;
}

void 
COPIntEnergyHist::output(xmlw::XmlStream& XML)
{
  XML << xmlw::tag("EnergyHist");
  intEnergyHist.outputClearHistogram(XML, Sim->dynamics.units().unitEnergy());
  XML << xmlw::endtag("EnergyHist");
}
