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

#include "outputplugin.hpp"
#include <boost/foreach.hpp>
#include "../dynamics/include.hpp"
#include "../datatypes/complex.hpp"
#include "../base/is_exception.hpp"
#include "../base/constants.hpp"
#include "../simulation/particle.hpp"
#include "../extcode/xmlParser.h"
#include "include.hpp"
#include "correlations/include.hpp"
#include "general/include.hpp"
#include <boost/tokenizer.hpp>

COutputPlugin::COutputPlugin(const DYNAMO::SimData* tmp, const char *aName, unsigned char order, const char *aColor):
  SimBase_const(tmp, aName, aColor),
  updateOrder(order)
{
  I_cout() << "Loaded";
}

void 
COutputPlugin::output(xmlw::XmlStream&) 
{}

void 
COutputPlugin::periodicOutput() 
{}

DYNAMO::Colorise_Text_Stream_Operator 
COutputPlugin::I_Pcout() const
{
  return DYNAMO::Colorise_Text_Stream_Operator(IC_blue);
}

COutputPlugin* 
COutputPlugin::getPlugin(std::string Details, const DYNAMO::SimData* Sim)
{
  XMLNode XML = XMLNode::createXMLTopNode("Plugin");

  typedef boost::tokenizer<boost::char_separator<char> > 
    tokenizer;
  
  boost::char_separator<char> DetailsSep(":");
  boost::char_separator<char> OptionsSep(",");
  boost::char_separator<char> ValueSep("=", "", boost::keep_empty_tokens);

  tokenizer tokens(Details, DetailsSep);
  tokenizer::iterator details_iter = tokens.begin();
  
  XML.addAttribute("Type", details_iter->c_str());
  ++details_iter;
  if (details_iter != tokens.end())
    {
      tokenizer option_tokens(*details_iter, OptionsSep);

      if (++details_iter != tokens.end())
	D_throw() << "Two colons in outputplugin options " << *details_iter;

      for (tokenizer::iterator options_iter = option_tokens.begin();
	   options_iter != option_tokens.end(); ++options_iter)
	{ 
	  tokenizer value_tokens(*options_iter, ValueSep);
	  
	  tokenizer::iterator value_iter = value_tokens.begin();
	  std::string opName(*value_iter);
	  if (++value_iter == value_tokens.end())
	    D_throw() << "Option " << opName << "with no value!";

	  XML.addAttribute(opName.c_str(), value_iter->c_str());	  
	}
    }

  return getPlugin(XML,Sim);
}

template<class T> COutputPlugin* 
COutputPlugin::testGeneratePlugin(const DYNAMO::SimData* Sim, const XMLNode& XML)
{
  try {
    Sim->getOutputPlugin<T>();
  } catch (std::exception&)
    {
      return new T(Sim, XML);
    }  

  //It's already in the simulation
  D_throw() << "Plugin is already loaded";
}

COutputPlugin* 
COutputPlugin::getPlugin(const XMLNode& XML, const DYNAMO::SimData* Sim)
{
  std::string Name = XML.getAttribute("Type");

  if (!Name.compare("MSD"))
    return testGeneratePlugin<COPMSD>(Sim, XML);
  else if (!Name.compare("PeriodicMSD"))
    return testGeneratePlugin<COPPeriodicMSD>(Sim, XML);
  else if (!Name.compare("EstTime"))
    return testGeneratePlugin<COPETA>(Sim, XML);
  else if (!Name.compare("ReplexTrace"))
    return testGeneratePlugin<COPReplexTrace>(Sim, XML);
  else if (!Name.compare("IntEnergyHist"))
    return testGeneratePlugin<COPIntEnergyHist>(Sim, XML);
  else if (!Name.compare("RadiusGyration"))
    return testGeneratePlugin<COPRGyration>(Sim, XML);
  else if (!Name.compare("Torsion"))
    return testGeneratePlugin<COPCTorsion>(Sim, XML);
  else if (!Name.compare("Geomview"))
    return testGeneratePlugin<COPGeomview>(Sim, XML);
  else if (!Name.compare("KEnergy"))
    return testGeneratePlugin<COPKEnergy>(Sim, XML);
  else if (!Name.compare("UEnergy"))
    return testGeneratePlugin<COPUEnergy>(Sim, XML);
  else if (!Name.compare("Misc"))
    return testGeneratePlugin<COPMisc>(Sim, XML);
  else if (!Name.compare("TinkerXYZ"))
    return testGeneratePlugin<COPTinkerXYZ>(Sim, XML);
  else if (!Name.compare("PackingFraction"))
    return testGeneratePlugin<COPPackingFraction>(Sim, XML);
  else if (!Name.compare("CollisionMatrix"))
    return testGeneratePlugin<COPCollMatrix>(Sim, XML);
  else if (!Name.compare("RdotV"))
    return testGeneratePlugin<COPRdotV>(Sim, XML);
  else if (!Name.compare("Momentum"))
    return testGeneratePlugin<COPMomentum>(Sim, XML);
  else if (!Name.compare("QMGA"))
    return testGeneratePlugin<COPQMGA>(Sim, XML);
#ifdef DYNAMO_VTK
  else if (!Name.compare("VTK"))
    return testGeneratePlugin<COPVTK>(Sim, XML);
#endif
  else if (!Name.compare("Povray"))
    return testGeneratePlugin<COPPovray>(Sim, XML);
  else if (!Name.compare("ContactMap"))
    return testGeneratePlugin<COPCContactMap>(Sim, XML);
  else if (!Name.compare("OverlapTester"))
    return testGeneratePlugin<COPOverlapTest>(Sim, XML);
  /*else if (!Name.compare("CollDistCheck"))
    return testGeneratePlugin<COPCollDistCheck>(Sim, XML);*/
  else if (!Name.compare("ChainBondAngles"))
    return testGeneratePlugin<COPChainBondAngles>(Sim, XML);
  else if (!Name.compare("ChainBondLength"))
    return testGeneratePlugin<COPChainBondLength>(Sim, XML);
  else if (!Name.compare("ReverseEventsCheck"))
    return testGeneratePlugin<COPReverseEventsCheck>(Sim, XML);
  else if (!Name.compare("VACF"))
    return testGeneratePlugin<COPVACF>(Sim,XML);
  else if (!Name.compare("ViscosityE"))
    return testGeneratePlugin<COPViscosityE>(Sim, XML);
  else if (!Name.compare("ThermalConductivityE"))
    return testGeneratePlugin<COPThermalConductivityE>(Sim, XML);
  else if (!Name.compare("MutualDiffusionGK"))
    return testGeneratePlugin<COPMutualDiffusionGK>(Sim, XML);
  else if (!Name.compare("ThermalDiffusionE"))
    return testGeneratePlugin<COPThermalDiffusionE>(Sim, XML);
  else if (!Name.compare("MFL"))
    return testGeneratePlugin<COPMFL>(Sim, XML);
  else if (!Name.compare("MFT"))
    return testGeneratePlugin<COPMFT>(Sim, XML);
  else if (!Name.compare("CollEnergyChange"))
    return testGeneratePlugin<COPCollEnergyChange>(Sim, XML);
  else if (!Name.compare("VelDist"))
    return testGeneratePlugin<COPVelDist>(Sim, XML);
  else if (!Name.compare("CollisionCorrelators"))
    return testGeneratePlugin<COPCollisionCorrelator>(Sim, XML);
#ifndef CBT
  else if (!Name.compare("BoundedPQStats"))
    return testGeneratePlugin<COPBoundedQStats>(Sim, XML);
#endif
  else 
    D_throw() << "Unknown type of OutputPlugin encountered\n"
	      << Name;
}
