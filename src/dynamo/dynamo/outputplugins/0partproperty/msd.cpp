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

#include <dynamo/outputplugins/0partproperty/msd.hpp>
#include <dynamo/include.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <boost/foreach.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
  OPMSD::OPMSD(const dynamo::Simulation* tmp, const magnet::xml::Node&):
    OutputPlugin(tmp,"MSD")
  {}

  OPMSD::~OPMSD()
  {}

  void
  OPMSD::initialise()
  {
    initPos.clear();
    initPos.resize(Sim->N);
  
    for (size_t ID = 0; ID < Sim->N; ++ID)
      initPos[ID] = Sim->particleList[ID].getPosition();
  }

  void
  OPMSD::output(magnet::xml::XmlStream &XML)
  {
    //Required to get the correct results
    Sim->dynamics->updateAllParticles();
  
    XML << magnet::xml::tag("MSD");
  
    BOOST_FOREACH(const shared_ptr<Species>& sp, Sim->species)
      {
	double MSD(calcMSD(*(sp->getRange())));
      
	XML << magnet::xml::tag("Species")
	    << magnet::xml::attr("Name") << sp->getName()
	    << magnet::xml::attr("val") << MSD
	    << magnet::xml::attr("diffusionCoeff") 
	    << MSD * Sim->units.unitTime() / (2 * NDIM * Sim->dSysTime)
	    << magnet::xml::endtag("Species");
      }

    if (!Sim->topology.empty())
      {
	XML << magnet::xml::tag("Structures");

	BOOST_FOREACH(const shared_ptr<Topology>& topo, Sim->topology)
	  {
	    double MSD(calcStructMSD(*topo));

	    XML << magnet::xml::tag("Structure")
		<< magnet::xml::attr("Name") << topo->getName()
		<< magnet::xml::attr("val") << MSD
		<< magnet::xml::attr("diffusionCoeff") 
		<< MSD * Sim->units.unitTime() / (2 * NDIM * Sim->dSysTime)
		<< magnet::xml::endtag("Structure");
	  }

	XML << magnet::xml::endtag("Structures");
      }

    XML << magnet::xml::endtag("MSD");
  }

  double
  OPMSD::calcMSD(const Range& range) const
  {
    double acc = 0.0;

    BOOST_FOREACH(const size_t ID, range)
      acc += (Sim->particleList[ID].getPosition() - initPos[ID]).nrm2();
  
    return acc / (range.size() * Sim->units.unitArea());
  }

  double
  OPMSD::calcStructMSD(const Topology& Itop) const
  {
    //Required to get the correct results
    Sim->dynamics->updateAllParticles();

    double acc = 0.0;
    BOOST_FOREACH(const shared_ptr<Range>& molRange, Itop.getMolecules())
      {
	Vector  origPos(0,0,0), currPos(0,0,0);
	double totmass = 0.0;
	BOOST_FOREACH(const unsigned long& ID, *molRange)
	  {
	    double pmass = Sim->species[Sim->particleList[ID]]->getMass(ID);

	    totmass += pmass;
	    currPos += Sim->particleList[ID].getPosition() * pmass;
	    origPos += initPos[ID] * pmass;
	  }
      
	currPos /= totmass;
	origPos /= totmass;

	acc += (currPos - origPos).nrm2();
      }

    acc /= Itop.getMoleculeCount() * Sim->units.unitArea();
  
    return acc;
  }
}
