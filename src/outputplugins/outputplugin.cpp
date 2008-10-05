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
COutputPlugin::getPlugin(std::string Name, const DYNAMO::SimData* Sim)
{
  XMLNode XML = XMLNode::createXMLTopNode("Plugin");
  XML.addAttribute("Type", Name.c_str());
  return getPlugin(XML,Sim);
}

COutputPlugin* 
COutputPlugin::getPlugin(const XMLNode& XML, const DYNAMO::SimData* Sim)
{
  std::string Name = XML.getAttribute("Type");

  if (!Name.compare("MSD"))
    return new COPMSD(Sim);
  else if (!Name.compare("PeriodicMSD"))
    return new COPPeriodicMSD(Sim);
  else if (!Name.compare("EstTime"))
    return new COPETA(Sim);
  else if (!Name.compare("ReplexTrace"))
    return new COPReplexTrace(Sim);
  else if (!Name.compare("IntEnergyHist"))
    return new COPIntEnergyHist(Sim);
  else if (!Name.compare("RadiusGyration"))
    return new COPRGyration(Sim);
  else if (!Name.compare("Torsion"))
    return new COPCTorsion(Sim);
  else if (!Name.compare("Geomview"))
    return new COPGeomview(Sim);
  else if (!Name.compare("KEnergy"))
    return new COPKEnergy(Sim);
  else if (!Name.compare("UEnergy"))
    return new COPUEnergy(Sim);
  else if (!Name.compare("Misc"))
    return new COPMisc(Sim);
  else if (!Name.compare("TinkerXYZ"))
    return new COPTinkerXYZ(Sim);
  else if (!Name.compare("ThermalConductivity"))
    return new COPThermalCon(Sim, XML);
  else if (!Name.compare("ThermalDiffusion"))
    return new   COPThermalDiffusion(Sim, XML);
  else if (!Name.compare("Viscosity"))
    return new COPViscosity(Sim, XML);
  else if (!Name.compare("MutualDiffusion"))
    return new COPMutualDiffusion(Sim, XML);
  else if (!Name.compare("PackingFraction"))
    return new COPPackingFraction(Sim);
  else if (!Name.compare("CollisionMatrix"))
    return new COPCollMatrix(Sim);
  else if (!Name.compare("RdotV"))
    return new COPRdotV(Sim);
  else if (!Name.compare("Momentum"))
    return new COPMomentum(Sim);
  else if (!Name.compare("QMGA"))
    return new COPQMGA(Sim);
#ifdef DYNAMO_VTK
  else if (!Name.compare("VTK"))
    return new COPVTK(Sim);
#endif
  else if (!Name.compare("VACF"))
    return new COPVACF(Sim,XML);
  else if (!Name.compare("Povray"))
    return new COPPovray(Sim);
  else if (!Name.compare("ContactMap"))
    return new COPCContactMap(Sim);
  else if (!Name.compare("OverlapTester"))
    return new COPOverlapTest(Sim);
  else if (!Name.compare("CollDistCheck"))
    return new COPCollDistCheck(Sim);
  else if (!Name.compare("ChainBondAngles"))
    return new COPChainBondAngles(Sim);
  else if (!Name.compare("ChainBondLength"))
    return new COPChainBondLength(Sim);
  else if (!Name.compare("ReverseEventsCheck"))
    return new COPReverseEventsCheck(Sim);
#ifndef CBT
  else if (!Name.compare("BoundedPQStats"))
    return new COPBoundedQStats(Sim);
#endif
  else 
    D_throw() << "Unknown type of OutputPlugin encountered\n"
	      << Name;
}
