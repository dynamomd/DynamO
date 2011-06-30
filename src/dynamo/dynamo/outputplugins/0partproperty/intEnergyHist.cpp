/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
#include "../../dynamics/include.hpp"
#include "../../dynamics/liouvillean/NewtonMCL.hpp"
#include "../../base/is_simdata.hpp"
#include "../1partproperty/uenergy.hpp"
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>
#include <fstream>

OPIntEnergyHist::OPIntEnergyHist(const dynamo::SimData* tmp, const magnet::xml::Node& XML):
  OPCollTicker(tmp,"InternalEnergyHistogram", 10),//Before OPEnergy
  intEnergyHist(1.0),
  ptrOPEnergy(NULL),
  weight(0.0),
  binwidth(1.0)
{
  operator<<(XML);
}

void 
OPIntEnergyHist::operator<<(const magnet::xml::Node& XML)
{  
  try 
    {
      binwidth = XML.getAttribute("BinWidth").as<double>(1.0);
    }
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in OPIntEnergyHist";
    }  
}

void 
OPIntEnergyHist::initialise() 
{
  ptrOPEnergy = Sim->getOutputPlugin<OPUEnergy>();
  intEnergyHist = C1DWeightHistogram(binwidth * Sim->dynamics.units().unitEnergy());
}

void 
OPIntEnergyHist::changeSystem(OutputPlugin* EHist2)
{
  //Add the current data
  intEnergyHist.addVal(ptrOPEnergy->getSimU(), weight);
  //Same for the other histogram
  static_cast<OPIntEnergyHist*>(EHist2)->intEnergyHist.addVal
    (static_cast<OPIntEnergyHist*>(EHist2)->ptrOPEnergy->getSimU(), 
     static_cast<OPIntEnergyHist*>(EHist2)->weight);

  //Now swap over the data
  std::swap(Sim, static_cast<OPIntEnergyHist*>(EHist2)->Sim);

  //NEVER SWAP THE PLUGIN POINTERS! they don't change
  //std::swap(ptrOPEnergy, static_cast<OPIntEnergyHist*>(EHist2)->ptrOPEnergy);

  //Reset the weighting
  weight = 0.0;
  static_cast<OPIntEnergyHist*>(EHist2)->weight = 0.0;
}

void 
OPIntEnergyHist::stream(double dt)
{
  weight += dt;
}

void 
OPIntEnergyHist::ticker()
{
  intEnergyHist.addVal(ptrOPEnergy->getSimU(), weight);
  weight = 0.0;
}

boost::unordered_map<int, double> 
OPIntEnergyHist::getImprovedW() const
{
  boost::unordered_map<int, double> retval;

  bool isMC(Sim->dynamics.liouvilleanTypeTest<LNewtonianMC>());

  typedef std::pair<const long, double> lv1pair;
  BOOST_FOREACH(const lv1pair &p1, intEnergyHist.data.data)
    {
      double E = p1.first * intEnergyHist.data.binWidth;
      
      //Fetch the current W value
      double W = 0;
      
      if (isMC) W += static_cast<const LNewtonianMC&>(Sim->dynamics.getLiouvillean()).W(E);
      
      double Pc = static_cast<double>(p1.second)
	/ (intEnergyHist.data.binWidth * intEnergyHist.sampleCount 
	   * Sim->dynamics.units().unitEnergy());

      //We only try to optimize parts of the histogram with greater
      //than 1% probability
      if (Pc > 0.01)
	retval[lrint(E / intEnergyHist.data.binWidth)] = W + std::log(Pc);
    }
  
  //Now center the energy warps about 0 to not cause funny changes in the tails.
  typedef std::pair<const int, double> locpair;
  double avg = 0;
  BOOST_FOREACH(const locpair& p, retval)
    avg += p.second;

  avg /= retval.size();

  BOOST_FOREACH(locpair& p, retval)
    p.second -= avg;

  return retval;
}


void 
OPIntEnergyHist::output(magnet::xml::XmlStream& XML)
{
  XML << magnet::xml::tag("EnergyHist")
      << magnet::xml::attr("BinWidth") << binwidth;

  if (dynamic_cast<const dynamo::EnsembleNVT*>(Sim->ensemble.get()) != NULL)
    XML << magnet::xml::attr("T") 
	<< static_cast<const dynamo::EnsembleNVT&>(*(Sim->ensemble)).getReducedEnsembleVals()[2];
  
  intEnergyHist.outputClearHistogram(XML, Sim->dynamics.units().unitEnergy());
  
  if (Sim->dynamics.liouvilleanTypeTest<LNewtonianMC>())
    {
      dout << "Detected a Multi-canonical Liouvillean, outputting w parameters" << std::endl;
      const LNewtonianMC& liouvillean(static_cast<const LNewtonianMC&>(Sim->dynamics.getLiouvillean()));
  
#ifdef DYNAMO_DEBUG      
      if (!dynamic_cast<const dynamo::EnsembleNVT*>(Sim->ensemble.get()))
	M_throw() << "Multi-canonical simulations require an NVT ensemble";
#endif
      
      XML << magnet::xml::tag("PotentialDeformation")
	  << magnet::xml::attr("EnergyStep") << intEnergyHist.data.binWidth * Sim->dynamics.units().unitEnergy();
      
      typedef std::pair<const long, double> lv1pair;
      BOOST_FOREACH(const lv1pair &p1, intEnergyHist.data.data)
	{
	  double E = p1.first * intEnergyHist.data.binWidth;
	  
	  //Fetch the current W value
	  double W = liouvillean.W(E);
	  
	  double Pc = static_cast<double>(p1.second)
	    / (intEnergyHist.data.binWidth * intEnergyHist.sampleCount 
	       * Sim->dynamics.units().unitEnergy());
	  
	  XML << magnet::xml::tag("W")
	      << magnet::xml::attr("Energy") << E * Sim->dynamics.units().unitEnergy()
	      << magnet::xml::attr("Value") << W + std::log(Pc)
	      << magnet::xml::attr("OldValue") << W
	      << magnet::xml::endtag("W");
	}
      
      XML << magnet::xml::endtag("PotentialDeformation");
      
    }
  XML << magnet::xml::endtag("EnergyHist");

}
