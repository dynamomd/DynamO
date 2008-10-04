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

COPConfig::COPConfig(const DYNAMO::SimData* tmp):
  COutputPlugin(tmp,"XMLConfig")
{}

COPConfig::~COPConfig()
{ I_cout() << "Unloaded"; }

void
COPConfig::output(xmlw::XmlStream &XML)
{  
  Sim->Dynamics.Liouvillean().updateAllParticles();
  
  XML << std::setprecision(std::numeric_limits<Iflt>::digits10)  
      << xmlw::prolog() << xmlw::tag("DYNAMOconfig") 
      << xmlw::attr("version") << "1.0"
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
      << Sim->Dynamics
      << xmlw::tag("ParticleData");

  I_cout() << "Writing Particles ";
  boost::progress_display prog(Sim->lN);

  for (unsigned long i = 0; i < Sim->lN; ++i)
    {
      CParticle tmp(Sim->vParticleList[i]);
      Sim->Dynamics.BCs().setPBC(tmp.getPosition(), tmp.getVelocity());

      tmp.scaleVelocity(1.0 / Sim->Dynamics.units().unitVelocity());
      tmp.scalePosition(1.0 / Sim->Dynamics.units().unitLength());

      XML << tmp;
      ++prog;
    }
    
  XML << xmlw::endtag("ParticleData")
      << xmlw::endtag("DYNAMOconfig");
  
  I_cout() << "Configuration written out";
}

void 
COPConfig::fileOutput(const char *fileName)
{
  namespace io = boost::iostreams;

  io::filtering_ostream coutputFile;
  coutputFile.push(io::bzip2_compressor());
  coutputFile.push(io::file_sink(fileName));

  xmlw::XmlStream XML(coutputFile);
  output(XML);
}
