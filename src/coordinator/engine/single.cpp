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

#include "single.hpp"

CESingle::CESingle(const boost::program_options::variables_map& nVM, CThreadPool& tp):
  CEngine(nVM, "config.out.xml.bz2", "output.xml.bz2", tp),
  peekMode(false)
{}

void 
CESingle::peekData()
{
  peekMode = true;
  Simulation.simShutdown();
}

void
CESingle::runSimulation()
{
  do
    {
      if (peekMode)
	{
	  Simulation.setTrajectoryLength(vm["ncoll"].as<unsigned long long>());
	  Simulation.outputData("peek.data.xml.bz2");	    
	  peekMode = false;
	}
      
      Simulation.runSimulation();
    } 
  while (peekMode);
}

void
CESingle::initialisation()
{
  preSimInit();

  if (!(vm.count("config-file")) || 
      (vm["config-file"].as<std::vector<std::string> >().size() != 1))
    D_throw() << "You must only provide one input file in single mode";

  setupSim(Simulation, vm["config-file"].as<std::vector<std::string> >()[0]);

  Simulation.initialise();

  postSimInit(Simulation);

  std::cout << "\nMain: Loading Output Plugins";
  fflush(stdout);

  Simulation.initPlugins();

  if (vm.count("ticker-period"))
    Simulation.setTickerPeriod(vm["ticker-period"].as<Iflt>());

  if (vm.count("scale-ticker"))
    Simulation.scaleTickerPeriod(vm["scale-ticker"].as<Iflt>());
      
}

void
CESingle::outputData()
{
  Simulation.outputData(outputFormat.c_str());
}

void
CESingle::outputConfigs()
{
  Simulation.writeXMLfile(configFormat.c_str());
}

void 
CESingle::forceShutdown()
{
  Simulation.simShutdown();
}
