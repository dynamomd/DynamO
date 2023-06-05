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

#include <dynamo/outputplugins/tickerproperty/radialdist.hpp>
#include <dynamo/outputplugins/misc.hpp>
#include <dynamo/include.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <Eigen/Dense>

namespace dynamo {
  OPRadialDistribution::OPRadialDistribution(const dynamo::Simulation* tmp, 
					     const magnet::xml::Node& XML):
    OPTicker(tmp,"RadialDistribution"),
    binWidth(0.01),
    length(100),
    sampleCount(0),
    sample_energy(0),
    sample_energy_bin_width(0)
  { operator<<(XML); }

  void 
  OPRadialDistribution::operator<<(const magnet::xml::Node& XML)
  {
    try {
      if (XML.hasAttribute("BinWidth"))
	    binWidth = XML.getAttribute("BinWidth").as<double>();
        binWidth *= Sim->units.unitLength();
    
      if (XML.hasAttribute("Length"))
	    length = XML.getAttribute("Length").as<size_t>();
      else {
        size_t mindir = 0;
        for (size_t iDim = 1; iDim < NDIM; ++iDim)
          if (Sim->primaryCellSize[iDim] > Sim->primaryCellSize[mindir])
            mindir = iDim;
      
        //Times 2 as the max dist is half a box length
        //+2 for rounding truncation and for a zero bin       
        length = 2 + static_cast<size_t>(Sim->primaryCellSize[mindir] / (2 * binWidth));
      }

      if (XML.hasAttribute("SampleEnergy")) {
        sample_energy = XML.getAttribute("SampleEnergy").as<double>() 
          * Sim->units.unitEnergy();
        sample_energy_bin_width = 1 / Sim->units.unitEnergy();
      
        if (XML.hasAttribute("SampleEnergyWidth"))
          sample_energy_bin_width = XML.getAttribute("SampleEnergyWidth").as<double>() * Sim->units.unitEnergy();
      }
      
      dout << "BinWidth = " << binWidth / Sim->units.unitLength()
	         << "\nLength = " << length;
    }
    catch (std::exception& excep) {
	  M_throw() << "Error while parsing output plugin options\n" << excep.what();
    }
  }

  void 
  OPRadialDistribution::initialise()
  {
    gr_accumulator.resize(Sim->species.size(), 
		  std::vector<std::vector<unsigned long> >(Sim->species.size(), std::vector<unsigned long>(length, 0)));
    
    accumulator.resize(Sim->species.size() * Sim->species.size() * length, 0u);

    moments.resize(N_moments * Sim->species.size() * Sim->species.size() * length, 0u);

    if (!(Sim->getOutputPlugin<OPMisc>()))
      M_throw() << "Radial Distribution requires the Misc output plugin";

    ticker();
  }

  void 
  OPRadialDistribution::ticker()
  {
    //A test to ensure we only sample at a target energy (if
    //specified)
    if (sample_energy_bin_width)
      {
	if (std::abs(sample_energy - Sim->getOutputPlugin<OPMisc>()->getConfigurationalU())
	    > sample_energy_bin_width * 0.5)
	  return;
	else
	  dout << "Sampling radial distribution as configurational energy is" 
	       << Sim->getOutputPlugin<OPMisc>()->getConfigurationalU()
	    / Sim->units.unitEnergy()
	       << " sample_energy is "
	       << sample_energy / Sim->units.unitEnergy()
	       << " and sample_energy_bin_width is "
	       << sample_energy_bin_width / Sim->units.unitEnergy()
	       << std::endl;
      }
    
    ++sampleCount;
  
    std::fill(accumulator.begin(), accumulator.end(), 0u);

      for (const shared_ptr<Species>& sp1 : Sim->species)
        for (const shared_ptr<Species>& sp2 : Sim->species)
          for (const size_t& p1 : *sp1->getRange())
            for (const size_t& p2 : *sp2->getRange())
              {
                Vector rij = Sim->particles[p1].getPosition() - Sim->particles[p2].getPosition();
                Sim->BCs->applyBC(rij);
                const size_t i = static_cast<size_t>(rij.nrm() / binWidth + 0.5);
                if (i < length) {
                  ++gr_accumulator[sp1->getID()][sp2->getID()][i];
                  ++accumulator[(sp1->getID() * Sim->species.size() + sp2->getID()) * length + i];
                }
              }
  }

  std::vector<std::pair<double, double> > 
  OPRadialDistribution::getgrdata(size_t species1ID, size_t species2ID) const
  {
    std::vector<std::pair<double, double> > retval;
    const double density = (Sim->species[species2ID]->getCount() - (species1ID == species2ID)) / Sim->getSimVolume();
    const size_t originsTaken = sampleCount * Sim->species[species2ID]->getCount();
    
    //Skip the zero bin
    retval.reserve(length);
    for (size_t i = 0; i < length; ++i)
      {
	const double radius = binWidth * i;
	const double volshell =  M_PI * (4.0 * binWidth * radius * radius + binWidth * binWidth * binWidth / 3.0);
	const double GR = static_cast<double>(gr_accumulator[species1ID][species2ID][i]) / (density * originsTaken * volshell);
	retval.push_back(std::pair<double, double>(radius, GR));
      }

    return retval;
  }

  void
  OPRadialDistribution::output(magnet::xml::XmlStream& XML)
  {
    XML << magnet::xml::tag("RadialDistribution")
	      << magnet::xml::attr("SampleCount")
      	<< sampleCount;
  
    for (const shared_ptr<Species>& sp1 : Sim->species)
      for (const shared_ptr<Species>& sp2 : Sim->species)
        {
	        const double density = (sp2->getCount() - (sp1 == sp2)) / Sim->getSimVolume();
	        const size_t originsTaken = sampleCount * sp1->getCount();

	        XML << magnet::xml::tag("Species")
	            << magnet::xml::attr("Name1")
	            << sp1->getName()
	            << magnet::xml::attr("Name2")
	            << sp2->getName()
	            << magnet::xml::attr("Samples") << originsTaken
	            << magnet::xml::chardata();

	        //Skip the zero bin
	        for (size_t i = 1; i < length; ++i)
	        {
	          const double radius = binWidth * i;
	          const double volshell =  M_PI * (4.0 * binWidth * radius * radius + binWidth * binWidth * binWidth / 3.0);
	          const double GR = static_cast<double>(gr_accumulator[sp1->getID()][sp2->getID()][i]) / (density * originsTaken * volshell);
	          XML << radius / Sim->units.unitLength() << " " << GR << "\n";
	        }
	        XML << magnet::xml::endtag("Species");
        }

    XML << magnet::xml::endtag("RadialDistribution");
  }
}
