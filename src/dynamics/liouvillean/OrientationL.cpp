/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>
    Copyright (C) 2011  Sebastian Gonzalez <tsuresuregusa@gmail.com>

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

#include "../species/inertia.hpp"
#include "OrientationL.hpp"
#include "../2particleEventData.hpp"
#include "../../datatypes/vector.xml.hpp"
#include "../units/units.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/iostreams/filter/base64.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/chain.hpp>
#include <boost/iostreams/device/stream_sink.hpp>
#include <boost/iostreams/device/stream_source.hpp>
#include <boost/iostreams/filter/base64cleaner.hpp>
#include <boost/iostreams/filter/linewrapout.hpp>

void 
LNOrientation::initialise() 
{
  Liouvillean::initialise();

  double sumEnergy(0.0);

  BOOST_FOREACH(const Particle& part, Sim->particleList)  
    sumEnergy += Sim->dynamics.getSpecies(part).getScalarMomentOfInertia(part.getID())
    * orientationData[part.getID()].angularVelocity.nrm2();
  
  //Check if any of the species are overridden
  bool hasInertia(false);
  BOOST_FOREACH(const magnet::ClonePtr<Species>& spec, Sim->dynamics.getSpecies())
    if (dynamic_cast<const SpInertia*>(spec.get_ptr()) != NULL)
      hasInertia = true;

  if (!hasInertia)
    M_throw() << "No species have inertia, using the orientational liouvillean is pointless";

  sumEnergy *= 0.5 / Sim->dynamics.units().unitEnergy();
  
  I_cout() << "System Rotational Energy " << sumEnergy
	   << "\nRotational kT " << sumEnergy / Sim->N;
}

void
LNOrientation::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Type")
      << "NOrientation";
}
  
void 
LNOrientation::loadParticleXMLData(const magnet::xml::Node& XML)
{
  Liouvillean::loadParticleXMLData(XML);

  orientationData.resize(Sim->N);
  
  size_t i(0);
  for (magnet::xml::Node node = XML.getNode("ParticleData").getNode("Pt"); 
       node.valid(); ++node, ++i)
    {
      orientationData[i].orientation << node.getNode("U");
      orientationData[i].angularVelocity << node.getNode("O");
      
      double oL = orientationData[i].orientation.nrm();
      
      if (!(oL > 0.0))
	M_throw() << "Particle ID " << i 
		  << " orientation vector is zero!";
      
      //Makes the vector a unit vector
      orientationData[i].orientation /= oL;
    }
}
 
