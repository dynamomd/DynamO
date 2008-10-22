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

#include "collisionCorrelator.hpp"
#include "../../dynamics/2particleEventData.hpp"
#include "../../dynamics/include.hpp"

COPCollisionCorrelator::COPCollisionCorrelator(const DYNAMO::SimData* t1,
					       const XMLNode& XML):
  COP2PP(t1,"CollisionCorrelator"),
  collisionHistoryLength(20)
{ operator<<(XML); }

void 
COPCollisionCorrelator::operator<<(const XMLNode& XML)
{
  try {
    if (XML.isAttributeSet("length"))
      collisionHistoryLength 
	= boost::lexical_cast<size_t>
	(XML.getAttribute("length"));
      }
  catch (std::exception& excep)
    {
      D_throw() << "Error while parsing " << name << "options\n"
		<< excep.what();
    }
}

void 
COPCollisionCorrelator::initialise()
{
  //Set the history size and clear it
  partnerHist.resize
    (Sim->lN, boost::circular_buffer<size_t>(collisionHistoryLength));

  //Set the initial value to its own id to indicate no history
  for (size_t id = 0; id < Sim->lN; ++id)
    BOOST_FOREACH(size_t& val, partnerHist[id])
      val = id;
  
  counter.resize(Sim->Dynamics.getSpecies().size(),
		 std::vector<std::pair<size_t,size_t> >
		 (collisionHistoryLength, std::pair<size_t, size_t>(0, 0)));
}

void 
COPCollisionCorrelator::performSweep(const size_t& id, const size_t& specid,
				     const size_t& nid)
{
  for (size_t n = 0; n < collisionHistoryLength; ++n)
    //Check that the history is valid
    if (partnerHist[id][n] != id)
      {
	++(counter[specid][n].second);

	if (partnerHist[id][n] == nid)
	  ++(counter[specid][n].first);	
      }

  //Update the history
  partnerHist[id]
    .push_front(nid);
}

void 
COPCollisionCorrelator::A2ParticleChange(const C2ParticleData& PDat)
{
  performSweep(PDat.particle1_.getParticle().getID(),
	       PDat.particle1_.getSpecies().getID(),
	       PDat.particle2_.getParticle().getID());
  
  performSweep(PDat.particle2_.getParticle().getID(),
	       PDat.particle2_.getSpecies().getID(),
	       PDat.particle1_.getParticle().getID());
}

void
COPCollisionCorrelator::output(xmlw::XmlStream &XML)
{
  XML << xmlw::tag("CollisionCorrelator");

  for (size_t id = 0; id < counter.size(); ++id)
    {
      XML << xmlw::tag("Species")
	  << xmlw::attr("Name")
	  << Sim->Dynamics.getSpecies()[id].getName()
	  << xmlw::chardata();

      for (size_t n = 0; n < collisionHistoryLength; ++n)
	XML << n + 1  << " " << static_cast<Iflt>(counter[id][n].first)
	  / static_cast<Iflt>(counter[id][n].second)  
	    << "\n";

      XML << xmlw::endtag("Species");
    }

  XML << xmlw::endtag("CollisionCorrelator");
}
