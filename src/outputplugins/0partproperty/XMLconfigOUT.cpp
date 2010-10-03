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

#include "XMLconfig.hpp"

# include <boost/iostreams/device/file.hpp>
# include <boost/iostreams/filtering_stream.hpp>
# include <boost/iostreams/filter/bzip2.hpp>
# include <boost/iostreams/chain.hpp>

#include <boost/foreach.hpp>
#include <boost/progress.hpp>
#include <iomanip>
#include "../../extcode/xmlwriter.hpp"
#include "../../simulation/particle.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_simdata.hpp"
#include "../../schedulers/scheduler.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"
#include "misc.hpp"
#include <iomanip>

OPConfig::OPConfig(const DYNAMO::SimData* tmp):
  OutputPlugin(tmp,"XMLConfig"),
  rounding(false),
  compressedOutput(true)
{}

OPConfig::~OPConfig()
{ I_cout() << "Unloaded"; }

void
OPConfig::output(xml::XmlStream &XML)
{  
  Sim->dynamics.getLiouvillean().updateAllParticles();

  XML << std::scientific
    //This has a minus one due to the digit in front of the decimal
    //An extra one is added if we're rounding
      << std::setprecision(std::numeric_limits<double>::digits10 - 1 - rounding)
      << xml::prolog() << xml::tag("DYNAMOconfig") 
      << xml::attr("version") << configFileVersion
      << xml::tag("Simulation")
      << xml::tag("Trajectory")
      << xml::attr("Coll") << Sim->endEventCount
      << xml::attr("nCollPrint") << Sim->eventPrintInterval;

  //Allow this block to fail if need be
  try {
    XML << xml::attr("lastMFT") 
	<< Sim->getOutputPlugin<OPMisc>()->getMFT();
  }
  catch (std::exception&)
    {}

  XML << xml::endtag("Trajectory")
      << *Sim->ensemble
      << xml::tag("Scheduler")
      << *Sim->ptrScheduler
      << xml::endtag("Scheduler")
      << xml::tag("History") 
      << xml::chardata()
      << Sim->ssHistory.str()
      << "\nRun for " << Sim->eventCount << " collisions"
      << xml::endtag("History") << xml::endtag("Simulation")
      << Sim->dynamics;

  Sim->dynamics.getLiouvillean().outputParticleXMLData(XML);

  XML << xml::endtag("DYNAMOconfig");

  I_cout() << "Configuration written out";
}

void 
OPConfig::fileOutput(const char *fileName)
{
  namespace io = boost::iostreams;

  if (compressedOutput) 
    {
      io::filtering_ostream coutputFile;
      coutputFile.push(io::bzip2_compressor());
      
      coutputFile.push(io::file_sink(fileName));
      
      xml::XmlStream XML(coutputFile);
      
      XML.setFormatXML(true);
      output(XML);
    }
  else
    {
      std::ofstream coutputFile(fileName, std::ios::out | std::ios::trunc);
      xml::XmlStream XML(coutputFile);
      XML.setFormatXML(true);
      output(XML);
    }
}

void 
OPConfig::setRounding()
{
  rounding = true;
}
