/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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

#include <dynamo/outputplugins/intEnergyHist.hpp>
#include <dynamo/include.hpp>
#include <dynamo/dynamics/multicanonical.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/outputplugins/misc.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <fstream>

namespace dynamo {
  OPIntEnergyHist::OPIntEnergyHist(const dynamo::Simulation* tmp, const magnet::xml::Node& XML):
    OutputPlugin(tmp,"InternalEnergyHistogram", 10),//Before OPEnergy
    intEnergyHist(1.0),
    binwidth(1.0)
  {
    operator<<(XML);
  }

  void 
  OPIntEnergyHist::eventUpdate(const IntEvent &event, 
			    const PairEventData &) 
  {
    intEnergyHist.addVal(_ptrOPMisc->getConfigurationalU(), event.getdt());
  }

  void 
  OPIntEnergyHist::eventUpdate(const GlobalEvent &event, const NEventData&) 
  {
    intEnergyHist.addVal(_ptrOPMisc->getConfigurationalU(), event.getdt());
  }

  void 
  OPIntEnergyHist::eventUpdate(const LocalEvent &event, const NEventData&) 
  {
    intEnergyHist.addVal(_ptrOPMisc->getConfigurationalU(), event.getdt());
  }

  void 
  OPIntEnergyHist::eventUpdate(const System&, const NEventData&, const double& dt)
  {
    intEnergyHist.addVal(_ptrOPMisc->getConfigurationalU(), dt);
  }

  void 
  OPIntEnergyHist::operator<<(const magnet::xml::Node& XML)
  {  
    if (XML.hasAttribute("BinWidth"))
      binwidth = XML.getAttribute("BinWidth").as<double>();
  }

  void 
  OPIntEnergyHist::initialise() 
  {
    _ptrOPMisc = Sim->getOutputPlugin<OPMisc>();
    if (!_ptrOPMisc) M_throw() << "IntEnergyHist requires Misc plugin!";
    intEnergyHist = magnet::math::HistogramWeighted<>(binwidth * Sim->units.unitEnergy());
  }

  void 
  OPIntEnergyHist::changeSystem(OutputPlugin* EHist2)
  {
    std::swap(Sim, static_cast<OPIntEnergyHist*>(EHist2)->Sim);
  }

  std::unordered_map<int, double>
  OPIntEnergyHist::getImprovedW() const
  {
    if (!std::dynamic_pointer_cast<const DynNewtonianMC>(Sim->dynamics))
      M_throw() << "Cannot improve an non-Multicanonical Dynamics";

  const DynNewtonianMC& dynamics = static_cast<const DynNewtonianMC&>(*Sim->dynamics);

    if (dynamics.getEnergyStep() != intEnergyHist.getBinWidth())
      M_throw() << "Cannot improve the W potential when there is a mismatch between the"
		<< " internal energy histogram and MC potential bin widths.";

    std::unordered_map<int, double> retval;

    typedef std::pair<const int, double> lv1pair;
    for (const lv1pair &p1 : intEnergyHist)
      {
	double E = p1.first * intEnergyHist.getBinWidth();
	
	double Pc = static_cast<double>(p1.second)
	  / (intEnergyHist.getBinWidth() * intEnergyHist.getSampleCount()
	     * Sim->units.unitEnergy());

	//We only try to optimize parts of the histogram with greater
	//than 1% probability
	if (Pc > 0.01)
	  retval[lrint(E / intEnergyHist.getBinWidth())]
	    = dynamics.W(E) + std::log(Pc);
      }
  
    //Now center the energy warps about 0 to not cause funny changes in the tails.
    typedef std::pair<const int, double> locpair;
    double avg = 0;
    for (const locpair& p : retval)
      avg += p.second;

    avg /= retval.size();

    for (locpair& p : retval)
      p.second -= avg;

    return retval;
  }


  void 
  OPIntEnergyHist::output(magnet::xml::XmlStream& XML)
  {
    XML << magnet::xml::tag("EnergyHist")
	<< magnet::xml::attr("BinWidth") << binwidth;

    if (std::dynamic_pointer_cast<const dynamo::EnsembleNVT>(Sim->ensemble))
      XML << magnet::xml::attr("T") 
	  << static_cast<const dynamo::EnsembleNVT&>(*(Sim->ensemble)).getReducedEnsembleVals()[2];
  
    intEnergyHist.outputClearHistogram(XML, Sim->units.unitEnergy());
  
    if (std::dynamic_pointer_cast<const DynNewtonianMC>(Sim->dynamics))
      {
	dout << "Detected a Multi-canonical Dynamics, outputting w parameters" << std::endl;
	const DynNewtonianMC& dynamics(static_cast<const DynNewtonianMC&>(*Sim->dynamics));
	
	XML << magnet::xml::tag("PotentialDeformation")
	    << magnet::xml::attr("EnergyStep")
	    << dynamics.getEnergyStep() * Sim->units.unitEnergy();
	
	typedef std::pair<const int, double> locpair;
	for (const locpair& p1 : dynamics.getMap())
	  XML << magnet::xml::tag("W")
	      << magnet::xml::attr("Energy")
	      << p1.first * dynamics.getEnergyStep() * Sim->units.unitEnergy()
	      << magnet::xml::attr("Value") << p1.second
	      << magnet::xml::endtag("W");
	
	XML << magnet::xml::endtag("PotentialDeformation");
      
      }
    XML << magnet::xml::endtag("EnergyHist");

  }
}
