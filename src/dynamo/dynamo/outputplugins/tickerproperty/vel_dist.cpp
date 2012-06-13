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

#include <dynamo/outputplugins/tickerproperty/vel_dist.hpp>
#include <dynamo/include.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <dynamo/outputplugins/1partproperty/kenergy.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>

namespace dynamo {
  OPVelDist::OPVelDist(const dynamo::Simulation* tmp, 
		       const magnet::xml::Node& XML):
    OPTicker(tmp,"VelDist"),
    binWidth(0.01)
  { operator<<(XML); }

  void 
  OPVelDist::operator<<(const magnet::xml::Node& XML)
  {
    try {
      if (XML.hasAttribute("binWidth"))
	binWidth = XML.getAttribute("binWidth").as<double>();
    }
    catch (std::exception& excep)
      {
	M_throw() << "Error while parsing " << name << "options\n"
		  << excep.what();
      }
  }

  void 
  OPVelDist::initialise()
  {
    _ptrOPEnergy = Sim->getOutputPlugin<OPKEnergy>();

    if (!_ptrOPEnergy)
      M_throw() << "VelDist requires UEnergy plugin";

    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      data[iDim].resize(Sim->species.size(), 
			magnet::math::Histogram<>(Sim->units.unitVelocity() 
				     * binWidth));
  }

  void 
  OPVelDist::ticker()
  {
    double factor = std::sqrt(Sim->units.unitMass() 
			      / _ptrOPEnergy->getCurrentkT());

    BOOST_FOREACH(const shared_ptr<Species>& sp, Sim->species)
      BOOST_FOREACH(const size_t& ID, *sp->getRange())
      for (size_t iDim = 0; iDim < NDIM; ++iDim)
	data[iDim][sp->getID()]
	  .addVal(Sim->particleList[ID].getVelocity()[iDim] * factor);
  }

  void
  OPVelDist::output(magnet::xml::XmlStream& XML)
  {
    XML << magnet::xml::tag("VelDist");
  
    for (size_t id = 0; id < Sim->species.size(); ++id)
      {
	XML << magnet::xml::tag("Species")
	    << magnet::xml::attr("Name")
	    << Sim->species[id]->getName();
     
	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  {
	    XML << magnet::xml::tag("Dimension")
		<< magnet::xml::attr("val")
		<< iDim;
	  
	    data[iDim][id].outputHistogram(XML, 1.0);
	  
	    XML << magnet::xml::endtag("Dimension");
	  }
      
	XML << magnet::xml::endtag("Species");
      }
  
    XML << magnet::xml::endtag("VelDist");
  }
}
