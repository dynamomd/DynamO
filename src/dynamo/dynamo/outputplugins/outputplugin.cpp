/*  dynamo:- Event driven molecular dynamics simulator
    http://www.dynamomd.org
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

#include <dynamo/outputplugins/outputplugin.hpp>
#include <dynamo/outputplugins/include.hpp>
#include <dynamo/outputplugins/include.hpp>
#include <dynamo/outputplugins/include.hpp>
#include <dynamo/outputplugins/tickerproperty/include.hpp>
#include <dynamo/include.hpp>
#include <dynamo/particle.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/tokenizer.hpp>

namespace dynamo {
  OutputPlugin::OutputPlugin(const dynamo::Simulation* tmp, const char *aName, unsigned char order):
    SimBase_const(tmp, aName),
    updateOrder(order)
  {
    dout << "Loaded" << std::endl;
  }

  void
  OutputPlugin::output(magnet::xml::XmlStream&)
  {}

  void
  OutputPlugin::periodicOutput()
  {}

  std::ostream&
  OutputPlugin::I_Pcout() const
  {
    return std::cout;
  }

  shared_ptr<OutputPlugin>
  OutputPlugin::getPlugin(std::string Details, const dynamo::Simulation* Sim)
  {
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

    rapidxml::xml_document<> doc;
    rapidxml::xml_node<> *node = doc.allocate_node(rapidxml::node_element, 
						   doc.allocate_string("OP"));

    boost::char_separator<char> DetailsSep(":");
    boost::char_separator<char> OptionsSep(",");
    boost::char_separator<char> ValueSep("=", "", boost::keep_empty_tokens);

    tokenizer tokens(Details, DetailsSep);
    tokenizer::iterator details_iter = tokens.begin();

    node->append_attribute(doc.allocate_attribute
			   ("Type", doc.allocate_string(details_iter->c_str())));

    ++details_iter;
    if (details_iter != tokens.end())
      {
	tokenizer option_tokens(*details_iter, OptionsSep);

	if (++details_iter != tokens.end())
	  M_throw() << "Two colons in outputplugin options " << *details_iter;

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

	    node->append_attribute(doc.allocate_attribute
				   (doc.allocate_string(opName.c_str()), 
				    doc.allocate_string(val.c_str())));
	  }
      }

    return getPlugin(magnet::xml::Node(node, NULL), Sim);
  }
  
  namespace {
    template<class T> shared_ptr<dynamo::OutputPlugin>
    testGeneratePlugin(const dynamo::Simulation* Sim, const magnet::xml::Node& XML)
    {
      if (!(Sim->getOutputPlugin<T>()))
	return shared_ptr<dynamo::OutputPlugin>(new T(Sim, XML));

      //It's already in the simulation
      M_throw() << "Plugin is already loaded";
    }
  }

  shared_ptr<OutputPlugin>
  OutputPlugin::getPlugin(const magnet::xml::Node& XML, const dynamo::Simulation* Sim)
  {
    std::string Name(XML.getAttribute("Type"));

    if (!Name.compare("MSD"))
      return testGeneratePlugin<OPMSD>(Sim, XML);
    else if (!Name.compare("PeriodicMSD"))
      return testGeneratePlugin<OPPeriodicMSD>(Sim, XML);
    else if (!Name.compare("ReplexTrace"))
      return testGeneratePlugin<OPReplexTrace>(Sim, XML);
    else if (!Name.compare("IntEnergyHist"))
      return testGeneratePlugin<OPIntEnergyHist>(Sim, XML);
    else if (!Name.compare("RadiusGyration"))
      return testGeneratePlugin<OPRGyration>(Sim, XML);
    else if (!Name.compare("Torsion"))
      return testGeneratePlugin<OPCTorsion>(Sim, XML);
    else if (!Name.compare("Misc"))
      return testGeneratePlugin<OPMisc>(Sim, XML);
    else if (!Name.compare("CollisionMatrix"))
      return testGeneratePlugin<OPCollMatrix>(Sim, XML);
    else if (!Name.compare("ContactMap"))
      return testGeneratePlugin<OPCContactMap>(Sim, XML);
    else if (!Name.compare("Contactmap"))
      return testGeneratePlugin<OPContactMap>(Sim, XML);
    else if (!Name.compare("OverlapTester"))
      return testGeneratePlugin<OPOverlapTest>(Sim, XML);
    else if (!Name.compare("ChainBondAngles"))
      return testGeneratePlugin<OPChainBondAngles>(Sim, XML);
    else if (!Name.compare("Trajectory"))
      return testGeneratePlugin<OPTrajectory>(Sim, XML);
    else if (!Name.compare("ChainBondLength"))
      return testGeneratePlugin<OPChainBondLength>(Sim, XML);
    else if (!Name.compare("VelDist"))
      return testGeneratePlugin<OPVelDist>(Sim, XML);
    else if (!Name.compare("VelProfile"))
      return testGeneratePlugin<OPVelProfile>(Sim, XML);
    else if (!Name.compare("RadialDistribution"))
      return testGeneratePlugin<OPRadialDistribution>(Sim, XML);
    else if (!Name.compare("MSDCorrelator"))
      return testGeneratePlugin<OPMSDCorrelator>(Sim, XML);
    else if (!Name.compare("VACF"))
      return testGeneratePlugin<OPVACF>(Sim, XML);
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
    else if (!Name.compare("MSDOrientational"))
      return testGeneratePlugin<OPMSDOrientational>(Sim, XML);
    else if (!Name.compare("MSDOrientationalCorrelator"))
      return testGeneratePlugin<OPMSDOrientationalCorrelator>(Sim, XML);
    else if (!Name.compare("OrientationalOrder"))
      return testGeneratePlugin<OPOrientationalOrder>(Sim, XML);
    else
      M_throw() << Name << ", Unknown type of OutputPlugin encountered";
  }
}
