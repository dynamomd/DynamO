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

#include "XMLconfig.hpp"
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/chain.hpp>
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

COPConfig::COPConfig(const DYNAMO::SimData* tmp):
  COutputPlugin(tmp,"XMLConfig"),
  rounding(false),
  compressedOutput(true)
{}

COPConfig::~COPConfig()
{ I_cout() << "Unloaded"; }

void
COPConfig::output(xmlw::XmlStream &XML)
{  
  Sim->Dynamics.Liouvillean().updateAllParticles();
  
  XML << std::scientific
    //This has a minus one due to the digit in front of the decimal
    //An extra one is added if we're rounding
      << std::setprecision(std::numeric_limits<Iflt>::digits10 - 1 - rounding)
      << xmlw::prolog() << xmlw::tag("DYNAMOconfig") 
      << xmlw::attr("version") << configFileVersion
      << xmlw::tag("Simulation")
      << xmlw::tag("Trajectory")
      << xmlw::attr("Coll") << Sim->lMaxNColl
      << xmlw::attr("nCollPrint") << Sim->lNPrint;

  //Allow this block to fail if need be
  try {
    XML << xmlw::attr("lastMFT") 
	<< Sim->getOutputPlugin<COPMisc>()->getMFT();
  }
  catch (std::exception&)
    {}

  XML << xmlw::endtag("Trajectory")
      << *Sim->Ensemble
      << xmlw::tag("Scheduler")
      << *Sim->ptrScheduler
      << xmlw::endtag("Scheduler")
      << xmlw::tag("History") 
      << xmlw::chardata()
      << Sim->ssHistory.str()
      << "\nRun for " << Sim->lNColl << " collisions"
      << xmlw::endtag("History") << xmlw::endtag("Simulation")
      << Sim->Dynamics;

  Sim->Dynamics.Liouvillean().outputParticleXMLData(XML);

  XML << xmlw::tag("AppendedBinaryData")
      << xmlw::chardata();

  Sim->Dynamics.Liouvillean().outputParticleBin64Data(XML.getUnderlyingStream());
  XML << "\n" 
      << xmlw::endtag("AppendedBinaryData")
      << xmlw::endtag("DYNAMOconfig");

  I_cout() << "Configuration written out";
}

void 
COPConfig::fileOutput(const char *fileName)
{
  namespace io = boost::iostreams;

  io::filtering_ostream coutputFile;

  if (compressedOutput) coutputFile.push(io::bzip2_compressor());
  
  coutputFile.push(io::file_sink(fileName));

  xmlw::XmlStream XML(coutputFile);
  output(XML);
}

void 
COPConfig::setRounding()
{
  rounding = true;
}
