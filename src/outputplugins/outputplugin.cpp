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

#include "include.hpp"
#include "outputplugin.hpp"
#include <boost/foreach.hpp>
#include "../dynamics/include.hpp"
#include "../datatypes/complex.hpp"
#include "../base/is_exception.hpp"
#include "../base/constants.hpp"
#include "../simulation/particle.hpp"
#include "../extcode/xmlParser.h"
#include "correlations/include.hpp"
#include "general/include.hpp"
#include <boost/tokenizer.hpp>

OutputPlugin::OutputPlugin(const DYNAMO::SimData* tmp, const char *aName, unsigned char order, const char *aColor):
  SimBase_const(tmp, aName, aColor),
  updateOrder(order)
{
  I_cout() << "Loaded";
}

void
OutputPlugin::output(xmlw::XmlStream&)
{}

void
OutputPlugin::periodicOutput()
{}

DYNAMO::Colorise_Text_Stream_Operator
OutputPlugin::I_Pcout() const
{
  return DYNAMO::Colorise_Text_Stream_Operator(IC_blue);
}

OutputPlugin*
OutputPlugin::getPlugin(std::string Details, const DYNAMO::SimData* Sim)
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
	  std::string val;
	  if (++value_iter == value_tokens.end())
	    //There is no value to save, must be a flag
	    val = "";
	  else
	    val = *value_iter;

	  XML.addAttribute(opName.c_str(), val.c_str());
	}
    }

  return getPlugin(XML,Sim);
}

template<class T> OutputPlugin*
OutputPlugin::testGeneratePlugin(const DYNAMO::SimData* Sim, const XMLNode& XML)
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

OutputPlugin*
OutputPlugin::getPlugin(const XMLNode& XML, const DYNAMO::SimData* Sim)
{
  std::string Name = XML.getAttribute("Type");

  {
    std::string XMLstring(XML.createXMLString());

    std::string::size_type pos = XMLstring.find_last_of("\n");

    if(pos != std::string::npos)
      {
	XMLstring.erase(pos);
      }

    std::cout << IC_purple << "\nOutputPluginParser:" << IC_reset
	      <<  " Parsing XML " << XMLstring;
  }

  if (!Name.compare("MSD"))
    return testGeneratePlugin<OPMSD>(Sim, XML);
  else if (!Name.compare("PeriodicMSD"))
    return testGeneratePlugin<OPPeriodicMSD>(Sim, XML);
  else if (!Name.compare("EstTime"))
    return testGeneratePlugin<OPETA>(Sim, XML);
  else if (!Name.compare("ReplexTrace"))
    return testGeneratePlugin<OPReplexTrace>(Sim, XML);
  else if (!Name.compare("IntEnergyHist"))
    return testGeneratePlugin<OPIntEnergyHist>(Sim, XML);
  else if (!Name.compare("RadiusGyration"))
    return testGeneratePlugin<OPRGyration>(Sim, XML);
  else if (!Name.compare("Torsion"))
    return testGeneratePlugin<OPCTorsion>(Sim, XML);
  else if (!Name.compare("Geomview"))
    return testGeneratePlugin<OPGeomview>(Sim, XML);
  else if (!Name.compare("Streamticker"))
    return testGeneratePlugin<OPStreamTicker>(Sim, XML);
  else if (!Name.compare("KEnergy"))
    return testGeneratePlugin<OPKEnergy>(Sim, XML);
  else if (!Name.compare("UEnergy"))
    return testGeneratePlugin<OPUEnergy>(Sim, XML);
  else if (!Name.compare("Misc"))
    return testGeneratePlugin<OPMisc>(Sim, XML);
  else if (!Name.compare("TinkerXYZ"))
    return testGeneratePlugin<OPTinkerXYZ>(Sim, XML);
  else if (!Name.compare("CollisionMatrix"))
    return testGeneratePlugin<OPCollMatrix>(Sim, XML);
  else if (!Name.compare("RdotV"))
    return testGeneratePlugin<OPRdotV>(Sim, XML);
  else if (!Name.compare("Momentum"))
    return testGeneratePlugin<OPMomentum>(Sim, XML);
  else if (!Name.compare("QMGA"))
    return testGeneratePlugin<OPQMGA>(Sim, XML);
  else if (!Name.compare("VTK"))
    return testGeneratePlugin<OPVTK>(Sim, XML);
  else if (!Name.compare("Povray"))
    return testGeneratePlugin<OPPovray>(Sim, XML);
  else if (!Name.compare("ContactMap"))
    return testGeneratePlugin<OPCContactMap>(Sim, XML);
  else if (!Name.compare("OverlapTester"))
    return testGeneratePlugin<OPOverlapTest>(Sim, XML);
  else if (!Name.compare("CollDistCheck"))
    return testGeneratePlugin<OPCollDistCheck>(Sim, XML);
  else if (!Name.compare("ChainBondAngles"))
    return testGeneratePlugin<OPChainBondAngles>(Sim, XML);
  else if (!Name.compare("Trajectory"))
    return testGeneratePlugin<OPTrajectory>(Sim, XML);
  else if (!Name.compare("ChainBondLength"))
    return testGeneratePlugin<OPChainBondLength>(Sim, XML);
  else if (!Name.compare("ReverseEventsCheck"))
    return testGeneratePlugin<OPReverseEventsCheck>(Sim, XML);
  else if (!Name.compare("VACF"))
    return testGeneratePlugin<OPVACF>(Sim,XML);
  else if (!Name.compare("ViscosityE"))
    return testGeneratePlugin<OPViscosityE>(Sim, XML);
  else if (!Name.compare("ViscosityCollisionalE"))
    return testGeneratePlugin<OPViscosityCollisionalE>(Sim, XML);
  else if (!Name.compare("ThermalConductivityE"))
    return testGeneratePlugin<OPThermalConductivityE>(Sim, XML);
  else if (!Name.compare("ThermalConductivitySpeciesSpeciesE"))
    return testGeneratePlugin<OPThermalConductivitySpeciesSpeciesE>(Sim, XML);
  else if (!Name.compare("MutualDiffusionGK"))
    return testGeneratePlugin<OPMutualDiffusionGK>(Sim, XML);
  else if (!Name.compare("MutualDiffusionE"))
    return testGeneratePlugin<OPMutualDiffusionE>(Sim, XML);
  else if (!Name.compare("ThermalDiffusionE"))
    return testGeneratePlugin<OPThermalDiffusionE>(Sim, XML);
  else if (!Name.compare("MFL"))
    return testGeneratePlugin<OPMFL>(Sim, XML);
  else if (!Name.compare("MFT"))
    return testGeneratePlugin<OPMFT>(Sim, XML);
  else if (!Name.compare("CollEnergyChange"))
    return testGeneratePlugin<OPCollEnergyChange>(Sim, XML);
  else if (!Name.compare("VelDist"))
    return testGeneratePlugin<OPVelDist>(Sim, XML);
  else if (!Name.compare("VelProfile"))
    return testGeneratePlugin<OPVelProfile>(Sim, XML);
  else if (!Name.compare("RadialDistribution"))
    return testGeneratePlugin<OPRadialDistribution>(Sim, XML);
  else if (!Name.compare("CollisionCorrelators"))
    return testGeneratePlugin<OPCollisionCorrelator>(Sim, XML);
  else if (!Name.compare("BoundedPQStats"))
    return testGeneratePlugin<OPBoundedQStats>(Sim, XML);
  else if (!Name.compare("MSDCorrelator"))
    return testGeneratePlugin<OPMSDCorrelator>(Sim, XML);
  else if (!Name.compare("RijVijComponents"))
    return testGeneratePlugin<OPRijVij>(Sim, XML);
  else if (!Name.compare("KEnergyTicker"))
    return testGeneratePlugin<OPKEnergyTicker>(Sim, XML);
  else if (!Name.compare("StructureImage"))
    return testGeneratePlugin<OPStructureImaging>(Sim, XML);
  else if (!Name.compare("EventEffects"))
    return testGeneratePlugin<OPEventEffects>(Sim, XML);
  else if (!Name.compare("SHCrystal"))
    return testGeneratePlugin<OPSHCrystal>(Sim, XML);
  else if (!Name.compare("SCParameter"))
    return testGeneratePlugin<OPSCParameter>(Sim, XML);
  else if (!Name.compare("CubeComponents"))
    return testGeneratePlugin<OPCubeComp>(Sim, XML);
  else if (!Name.compare("PlateMotion"))
    return testGeneratePlugin<OPPlateMotion>(Sim, XML);
  else if (!Name.compare("SelfDiffusionOrientationalGK"))
    return testGeneratePlugin<OPSelfDiffusionOrientationalGK>(Sim, XML);
  else if (!Name.compare("MSDOrientational"))
    return testGeneratePlugin<OPMSDOrientational>(Sim, XML);
  else if (!Name.compare("MSDOrientationalCorrelator"))
    return testGeneratePlugin<OPMSDOrientationalCorrelator>(Sim, XML);
  else if (!Name.compare("ChatteringCorrelator"))
    return testGeneratePlugin<OPChatteringCorrelator>(Sim, XML);
  else
    D_throw() << Name << ", Unknown type of OutputPlugin encountered";
}
