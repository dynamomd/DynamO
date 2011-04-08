/*  DYNAMO:- Event driven molecular dynamics simulator 
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

#include "XMLconfig.hpp"

#include "../dynamics/dynamics.hpp"
#include "../schedulers/scheduler.hpp"
#include "../dynamics/BC/LEBC.hpp"
#include "../simulation/particle.hpp"
#include "../dynamics/units/units.hpp"
#include "../base/is_ensemble.hpp"
#include "../base/is_simdata.hpp"
#include "../dynamics/liouvillean/liouvillean.hpp"

#include <magnet/xmlreader.hpp>

CIPConfig::CIPConfig(std::string fn, DYNAMO::SimData* Sim):
  CInputPlugin(Sim,"initXMLFile"),
  fileName(fn)
{}

void
CIPConfig::initialise()
{
  using namespace magnet::xml;
  Document doc(fileName);  
  Node mainNode = doc.getNode("DYNAMOconfig");

  {
    std::string version(mainNode.getAttribute("version"));
    
    I_cout() << "Parsing XML file v" << version;
    
    if (version != configFileVersion)
      M_throw() << "This version of the config file is obsolete"
		<< "\nThe current version is " << configFileVersion
		<< "\nPlease look at the XMLFILE.VERSION file in the root directory of the dynamo source."
	;
  }

  Node subNode= mainNode.getNode("Simulation");
  Node browseNode = subNode.getNode("Trajectory");
  
  if (browseNode.getAttribute("lastMFT").valid())
    Sim->lastRunMFT = browseNode.getAttribute("lastMFT").as<double>();

  Sim->ssHistory << subNode.getNode("History");

  I_cout() << "Loading dynamics";

  Sim->dynamics << mainNode;

  I_cout() << "Loading Scheduler";

  Sim->ptrScheduler 
    = CScheduler::getClass(subNode.getNode("Scheduler"), Sim);

  I_cout() << "Loading Ensemble";
  if (subNode.getNode("Ensemble").valid())
    Sim->ensemble.reset
      (DYNAMO::CEnsemble::getClass(subNode.getNode("Ensemble"), Sim));
  else
    //Try and determine the Ensemble
    try {
      Sim->dynamics.getSystem("Thermostat");
      Sim->ensemble.reset(new DYNAMO::CENVT(Sim));
    }
    catch (std::exception&)
      {
	Sim->ensemble.reset(new DYNAMO::CENVE(Sim));
      }

  I_cout() << "Loading Particle data";

  Sim->dynamics.getLiouvillean().loadParticleXMLData(mainNode);
  
  //Fixes or conversions once system is loaded
  Sim->lastRunMFT *= Sim->dynamics.units().unitTime();

  I_cout() << "Configuration loaded";
}
