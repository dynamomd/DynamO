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

#include "compressor.hpp"

void 
CECompressor::getOptions(boost::program_options::options_description& opts)
{
  boost::program_options::options_description ropts("Compression Engine");

  ropts.add_options()
    ("growth-rate",boost::program_options::value<Iflt>()->default_value(1.0),
     "Compression rate for the simulation")
    ("check-system", "Check that the system has not violated any interaction"
     " information")
    ("target-pack-frac",boost::program_options::value<Iflt>(),
     "Target packing fraction that compression has to attain")
    ;
  opts.add(ropts);
}

CECompressor::CECompressor(const boost::program_options::variables_map& nVM):
  CESingle(nVM)
{}

void 
CECompressor::preSimInit()
{
  CESingle::preSimInit();

  compressPlug.set_ptr(new CIPCompression
		       (&Simulation, vm["growth-rate"].as<Iflt>()));
}

void 
CECompressor::setupSim(CSimulation& Sim, const std::string filename)
{
  CESingle::setupSim(Sim,filename);
  
  //This adds a system event to prevent the cellular scheduler 
  //failing during compression
  compressPlug->MakeGrowth();
        
  compressPlug->CellSchedulerHack();
    
  if (vm.count("target-pack-frac"))
    compressPlug->limitPackingFraction
      (vm["target-pack-frac"].as<Iflt>());

}

void CECompressor::finaliseRun()
{
  compressPlug->RestoreSystem();

  if (vm.count("check-system"))
    compressPlug->checkOverlaps();
}
