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

#include <dynamo/outputplugins/1partproperty/MFL.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/dynamics/species/species.hpp>
#include <dynamo/dynamics/1particleEventData.hpp>
#include <dynamo/dynamics/units/units.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>

namespace dynamo {
  OPMFL::OPMFL(const dynamo::SimData* tmp, const magnet::xml::Node& XML):
    OP1PP(tmp,"MeanFreeLength", 250),
    binwidth(0.01)
  {
    operator<<(XML);
  }

  void 
  OPMFL::operator<<(const magnet::xml::Node& XML)
  {
    try 
      {
	if (XML.hasAttribute("binwidth"))
	  binwidth = XML.getAttribute("binwidth").as<double>();
      }
    catch (boost::bad_lexical_cast&)
      {
	M_throw() << "Failed a lexical cast in OPMFL";
      }
  }

  void
  OPMFL::initialise()
  {
    lastTime.resize(Sim->N, 0.0);
    data.resize(Sim->species.size(), 
		magnet::math::Histogram<>(Sim->units.unitLength() * binwidth));
  }

  void 
  OPMFL::A1ParticleChange(const ParticleEventData& PDat)
  {
    //We ignore stuff that hasn't had an event yet
    if (lastTime[PDat.getParticle().getID()] != 0.0)
      { 
	data[PDat.getSpecies().getID()]
	  .addVal(PDat.getParticle().getVelocity().nrm() 
		  * (Sim->dSysTime 
		     - lastTime[PDat.getParticle().getID()]));
      }
  
    lastTime[PDat.getParticle().getID()] = Sim->dSysTime;
  }

  void
  OPMFL::output(magnet::xml::XmlStream &XML)
  {
    XML << magnet::xml::tag("MFL");
  
    for (size_t id = 0; id < data.size(); ++id)
      {
	XML << magnet::xml::tag("Species")
	    << magnet::xml::attr("Name")
	    << Sim->species[id]->getName();

	data[id].outputHistogram(XML, 1.0 / Sim->units.unitLength());
      
	XML << magnet::xml::endtag("Species");
      }

    XML << magnet::xml::endtag("MFL");
  }
}
