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
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <boost/lexical_cast.hpp>
namespace io = boost::iostreams;

#include "../extcode/xmlParser.h"
#include "../base/is_exception.hpp"
#include "../dynamics/dynamics.hpp"
//#include "../schedulers/cellularLE.hpp"
#include "../schedulers/scheduler.hpp"
#include "../dynamics/BC/LEBC.hpp"
#include "../simulation/particle.hpp"
#include "../dynamics/units/units.hpp"
#include "../base/is_ensemble.hpp"
#include "../base/is_simdata.hpp"

CIPConfig::CIPConfig(std::string fn, DYNAMO::SimData* Sim):
  CInputPlugin(Sim,"initXMLFile"),
  fileName(fn)
{}

void
CIPConfig::initialise()
{
  XMLNode xMainNode;
  if (std::string(fileName.end()-4, fileName.end()) == ".xml")
    {
      if (!boost::filesystem::exists(fileName))
	D_throw() << "Could not open XML configuration file";

      I_cout() << "Uncompressed XML input file " << fileName << " loading";

      try {
	xMainNode=XMLNode::openFileHelper(fileName.c_str(), "DYNAMOconfig");
      }
      catch (std::exception&)
	{
	  xMainNode=XMLNode::openFileHelper(fileName.c_str(), "ISSSconfig");
	}
    }
  else if (std::string(fileName.end()-8, fileName.end()) == ".xml.bz2")
    {
      if (!boost::filesystem::exists(fileName))
	D_throw() << "Could not open XML configuration file";

      io::filtering_istream inputFile;
      inputFile.push(io::bzip2_decompressor());
      inputFile.push(io::file_source(fileName));
      //Copy file to a string
      std::string line, fileString;

      I_cout() << "Bzip compressed XML input file found\nDecompressing file "
	       << fileName;

      while(getline(inputFile,line)) 
	{
	  fileString.append(line);
	  fileString.append("\n");
	}

      I_cout() << "File Decompressed, parsing XML";
      fflush(stdout);
      
      XMLNode tmpNode = XMLNode::parseString(fileString.c_str());

      try {
	xMainNode = tmpNode.getChildNode("DYNAMOconfig");
      }
      catch (std::exception&)
	{
	  xMainNode = tmpNode.getChildNode("ISSSconfig");
	}
    }
  else
    D_throw() << "Unrecognised extension for input file";

  I_cout() << "Parsing XML file";
  XMLNode xSubNode= xMainNode.getChildNode("Simulation");
  XMLNode xBrowseNode = xSubNode.getChildNode("Trajectory");
  
  if (xBrowseNode.isAttributeSet("lastMFT"))
    Sim->lastRunMFT = atof(xBrowseNode.getAttribute("lastMFT"));

  xBrowseNode = xSubNode.getChildNode("History");
  Sim->ssHistory << xBrowseNode.getText();
  
  Sim->Dynamics << xMainNode;

  Sim->ptrScheduler = 
    CScheduler::getClass(xSubNode.getChildNode("Scheduler"),Sim);

  if (xSubNode.nChildNode("Ensemble"))
  Sim->Ensemble.reset
    (DYNAMO::CEnsemble::getClass(xSubNode.getChildNode("Ensemble"), Sim));
  else
    //Try and determine the Ensemble
    try {
      Sim->Dynamics.getSystem("Thermostat");
      Sim->Ensemble.reset(new DYNAMO::CENVT(Sim));
    }
    catch (std::exception&)
      {
	Sim->Ensemble.reset(new DYNAMO::CENVE(Sim));
      }

  xSubNode = xMainNode.getChildNode("ParticleData");
  unsigned long nPart = xSubNode.nChildNode("Particle");

  I_cout() << "Loading Particle Data ";
  fflush(stdout);

  int xml_iter = 0;
  bool outofsequence = false;
  
  if (nPart)
    {
      boost::progress_display prog(nPart);
      
      for (unsigned long i = 0; i < nPart; i++)
	{
	  xBrowseNode = xSubNode.getChildNode("Particle", &xml_iter);
	  
	  if (boost::lexical_cast<unsigned long>
	      (xBrowseNode.getAttribute("ID")) != i)
	    outofsequence = true;
	  
	  CParticle part(xBrowseNode, i);
	  part.scaleVelocity(Sim->Dynamics.units().unitVelocity());
	  part.scalePosition(Sim->Dynamics.units().unitLength());
	  Sim->vParticleList.push_back(part);
	  ++prog;
	}
    }
  else
    {  
      nPart = xSubNode.nChildNode("Pt");
      boost::progress_display prog(nPart);

      for (unsigned long i = 0; i < nPart; i++)
	{
	  xBrowseNode = xSubNode.getChildNode("Pt", &xml_iter);
	  
	  if (boost::lexical_cast<unsigned long>
	      (xBrowseNode.getAttribute("ID")) != i)
	    outofsequence = true;
	  
	  CParticle part(xBrowseNode, i);
	  part.scaleVelocity(Sim->Dynamics.units().unitVelocity());
	  part.scalePosition(Sim->Dynamics.units().unitLength());
	  Sim->vParticleList.push_back(part);
	  ++prog;
	}
      
    }

  I_cout() << Sim->vParticleList.size() 
	   << " Particles found";  
  
  if (outofsequence)
    I_cout() << IC_red << "Particle ID's out of sequence!\n"
	     << IC_red << "This can result in incorrect capture map loads etc.\n"
	     << IC_red << "Erase any capture maps in the configuration file so they are regenerated."
	     << IC_reset;

  //Fixes or conversions once system is loaded
  Sim->lastRunMFT *= Sim->Dynamics.units().unitTime();
}
