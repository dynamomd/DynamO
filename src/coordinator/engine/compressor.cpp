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

#include "compressor.hpp"

void 
ECompressingSimulation::getOptions(boost::program_options::options_description& opts)
{
  boost::program_options::options_description ropts("Compression Engine");

  ropts.add_options()
    ("growth-rate",boost::program_options::value<Iflt>()->default_value(1.0),
     "Compression rate for the simulation")
    ("check-system", "Check that the system has not violated any interaction"
     " information")
    ("target-pack-frac",boost::program_options::value<Iflt>(),
     "Target packing fraction that compression has to attain to exit")
    ("target-density",boost::program_options::value<Iflt>(),
     "Target number density that compression has to attain to exit")
    ;
  opts.add(ropts);
}

ECompressingSimulation::ECompressingSimulation(const boost::program_options::variables_map& nVM, 
					       magnet::thread::ThreadPool& tp):
  ESingleSimulation(nVM, tp)
{
  if (vm.count("target-pack-frac") && vm.count("target-density"))
    M_throw() << "Shouldn't specify both the packing fraction and density.";
    
}

void 
ECompressingSimulation::preSimInit()
{
  ESingleSimulation::preSimInit();

  compressPlug.set_ptr(new CIPCompression
		       (&simulation, vm["growth-rate"].as<Iflt>()));
}

void 
ECompressingSimulation::setupSim(Simulation& Sim, const std::string filename)
{
  ESingleSimulation::setupSim(Sim,filename);
  
  compressPlug->MakeGrowth();
            
  if (vm.count("target-pack-frac"))
    compressPlug->limitPackingFraction
      (vm["target-pack-frac"].as<Iflt>());
  else if (vm.count("target-density"))
    compressPlug->limitDensity
      (vm["target-density"].as<Iflt>());

  //This adds a system event to prevent the cellular scheduler 
  //failing during compression
  compressPlug->CellSchedulerHack();
}

void ECompressingSimulation::finaliseRun()
{
  compressPlug->RestoreSystem();

  if (vm.count("check-system"))
    compressPlug->checkOverlaps();
}
