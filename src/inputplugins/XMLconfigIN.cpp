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

#ifndef DYNAMO_CONDOR
# include <boost/iostreams/device/file.hpp>
# include <boost/iostreams/filtering_stream.hpp>
# include <boost/iostreams/filter/bzip2.hpp>
# include <boost/iostreams/chain.hpp>
namespace io = boost::iostreams;
#else
# include <fstream>
#endif

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

#include "../extcode/xmlParser.h"
#include "../base/is_exception.hpp"
#include "../dynamics/dynamics.hpp"
#include "../schedulers/scheduler.hpp"
#include "../dynamics/BC/LEBC.hpp"
#include "../simulation/particle.hpp"
#include "../dynamics/units/units.hpp"
#include "../base/is_ensemble.hpp"
#include "../base/is_simdata.hpp"
#include "../dynamics/liouvillean/liouvillean.hpp"


CIPConfig::CIPConfig(std::string fn, DYNAMO::SimData* Sim):
  CInputPlugin(Sim,"initXMLFile"),
  fileName(fn)
{}

void
CIPConfig::initialise()
{
  XMLNode xMainNode;

  if (!boost::filesystem::exists(fileName))
    D_throw() << "Could not open XML configuration file";

  //This scopes out the file objects
  {
#ifndef DYNAMO_CONDOR
    io::filtering_istream inputFile;
    
    if (std::string(fileName.end()-8, fileName.end()) == ".xml.bz2")
      {
	I_cout() << "Bzip compressed XML input file " << fileName << " loading";
	inputFile.push(io::bzip2_decompressor());
      }
    else if (std::string(fileName.end()-4, fileName.end()) == ".xml")
      I_cout() << "Uncompressed XML input file " << fileName << " loading";
    else
      D_throw() << "Unrecognized extension for input files";
    
    inputFile.push(io::file_source(fileName));
#else
    if (std::string(fileName.end()-8, fileName.end()) == ".xml.bz2")
      D_throw() << "Cannot load a compressed file when built for condor!";
    else if (std::string(fileName.end()-4, fileName.end()) == ".xml")
      I_cout() << "Uncompressed XML input file " << fileName << " loading";
    else
      D_throw() << "Unrecognized extension for input files";
    
    std::ifstream inputFile(fileName.c_str());
#endif
    std::cout.flush();
    
    //Copy file to a string
    std::string fileString; //The file is loaded into this
    std::string line;
    {
      while(getline(inputFile,line))
	fileString.append(line + "\n");
    }

    I_cout() << "File loaded, parsing XML";
    std::cout.flush();  

    XMLNode tmpNode = XMLNode::parseString(fileString.c_str());
    xMainNode = tmpNode.getChildNode("DYNAMOconfig");
  }


  {
    std::string version(xMainNode.getAttribute("version"));
    
    I_cout() << "Parsing XML file v" << version;
    
    if (version != configFileVersion)
      D_throw() << "This version of the config file is obsolete"
		<< "\nThe current version is " << configFileVersion;
  }

  XMLNode xSubNode= xMainNode.getChildNode("Simulation");
  XMLNode xBrowseNode = xSubNode.getChildNode("Trajectory");

  if (xBrowseNode.isAttributeSet("lastMFT"))
    Sim->lastRunMFT = atof(xBrowseNode.getAttribute("lastMFT"));

  xBrowseNode = xSubNode.getChildNode("History");
  Sim->ssHistory << xBrowseNode.getText();

  I_cout() << "Loading dynamics";

  Sim->dynamics << xMainNode;

  I_cout() << "Loading Scheduler";

  Sim->ptrScheduler = 
    CScheduler::getClass(xSubNode.getChildNode("Scheduler"),Sim);

  I_cout() << "Loading Ensemble";
  if (xSubNode.nChildNode("Ensemble"))
    Sim->Ensemble.reset
      (DYNAMO::CEnsemble::getClass(xSubNode.getChildNode("Ensemble"), Sim));
  else
    //Try and determine the Ensemble
    try {
      Sim->dynamics.getSystem("Thermostat");
      Sim->Ensemble.reset(new DYNAMO::CENVT(Sim));
    }
    catch (std::exception&)
      {
	Sim->Ensemble.reset(new DYNAMO::CENVE(Sim));
      }

  I_cout() << "Loading Particle data";

  xSubNode = xMainNode.getChildNode("ParticleData");

  Sim->dynamics.getLiouvillean().loadParticleXMLData(xMainNode);
  
  //Fixes or conversions once system is loaded
  Sim->lastRunMFT *= Sim->dynamics.units().unitTime();

  I_cout() << "Configuration loaded";
}
