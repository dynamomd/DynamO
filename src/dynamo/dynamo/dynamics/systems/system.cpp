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

#include <dynamo/dynamics/systems/system.hpp>
#include <dynamo/dynamics/systems/andersenThermostat.hpp>
#include <dynamo/dynamics/systems/rescale.hpp>
#include <dynamo/dynamics/systems/DSMCspheres.hpp>
#include <dynamo/dynamics/systems/RingDSMC.hpp>
#include <dynamo/dynamics/systems/umbrella.hpp>
#include <dynamo/dynamics/systems/visualizer.hpp>
#include <dynamo/dynamics/systems/sleep.hpp>
#include <dynamo/dynamics/../simulation/particle.hpp>
#include <dynamo/dynamics/interactions/intEvent.hpp>
#include <dynamo/dynamics/globals/globEvent.hpp>
#include <dynamo/dynamics/ranges/1RAll.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>

namespace dynamo {
  bool 
  System::operator<(const IntEvent& iEvent) const
  {
    return dt < iEvent.getdt();
  }

  bool 
  System::operator<(const GlobalEvent& gEvent) const
  {
    return dt < gEvent.getdt();
  }

  bool 
  System::operator<(const System& sEvent) const
  {
    return dt < sEvent.dt;
  }


  System::System(dynamo::SimData* tmp):
    SimBase(tmp, "SystemInteraction"),
    dt(HUGE_VAL)
  {
    type = NONE;
  }

  magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream& XML, 
				     const System& g)
  {
    g.outputXML(XML);
    return XML;
  }

  shared_ptr<System>
  System::getClass(const magnet::xml::Node& XML, dynamo::SimData* Sim)
  {
    if (!strcmp(XML.getAttribute("Type"),"Andersen"))
      return shared_ptr<System>(new SysAndersen(XML,Sim));
    else if (!strcmp(XML.getAttribute("Type"), "DSMCSpheres"))
      return shared_ptr<System>(new SysDSMCSpheres(XML, Sim));
    else if (!strcmp(XML.getAttribute("Type"), "Rescale"))
      return shared_ptr<System>(new SysRescale(XML, Sim));
    else if (!strcmp(XML.getAttribute("Type"), "RingDSMC"))
      return shared_ptr<System>(new SysRingDSMC(XML, Sim));
    else if (!strcmp(XML.getAttribute("Type"), "Umbrella"))
      return shared_ptr<System>(new SysUmbrella(XML, Sim));
    else if (!strcmp(XML.getAttribute("Type"), "Sleep"))
      return shared_ptr<System>(new SSleep(XML, Sim));
    else
      M_throw() << XML.getAttribute("Type")
		<< ", Unknown type of System event encountered";
  }
}
