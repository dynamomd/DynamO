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

#include <dynamo/outputplugins/tickerproperty/OrientationalOrder.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/simulation.hpp>
#include <boost/foreach.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <fstream>
#include <cmath>
#include <limits>

namespace dynamo {
  OPOrientationalOrder::OPOrientationalOrder(const dynamo::Simulation* tmp, const magnet::xml::Node& XML):
    OPTicker(tmp,"OrientationalOrder"), 
    _axis(1,0,0),
    _rg(1)
  {
    operator<<(XML);
  }

  void 
  OPOrientationalOrder::operator<<(const magnet::xml::Node& XML)
  {
    if (XML.hasAttribute("CutOffR"))
      _rg = XML.getAttribute("CutOffR").as<double>();

    _rg *= Sim->units.unitLength();

    dout << "Cut off radius set to " 
	 << _rg / Sim->units.unitLength() << std::endl;
  }


  void 
  OPOrientationalOrder::initialise() 
  { 
    ticker();
  }

  namespace {
    struct Neighbours {
      void addNeighbour(const Particle& p1, const size_t& p2)
      {
	if (p1.getID() != p2)
	  _neighbours.push_back(p2);
      }
      
      std::vector<size_t> _neighbours;
    };

    struct MagVec : public Vector
    {
      MagVec(const Vector& vec):
	Vector(vec)
      {}

      bool operator<(const MagVec& ovec) const
      {
	return nrm() < ovec.nrm();
      }

      bool operator>(const MagVec& ovec) const
      {
	return nrm() > ovec.nrm();
      }
    };
  }

  void 
  OPOrientationalOrder::ticker()
  {
    size_t count(0);
    ComplexNum sum(0,0);
    BOOST_FOREACH(const Particle& part, Sim->particles)
      {
	Neighbours nbs;
	
	std::auto_ptr<IDRange> ids(Sim->ptrScheduler->getParticleNeighbours(part));
	BOOST_FOREACH(const size_t& id1, *ids)
	  nbs.addNeighbour(part, id1);
	
	if (nbs._neighbours.size() >= 6)
	  {
	    std::vector<MagVec> bonds;
	    BOOST_FOREACH(const size_t& id2, nbs._neighbours)
	      {
		Vector bond = Sim->particles[id2].getPosition() - part.getPosition();
		Sim->BCs->applyBC(bond);
		bonds.push_back(bond);
	      }
	    std::sort(bonds.begin(), bonds.end());

	    if (bonds.size() > 6)
	      bonds.erase(bonds.begin()+6, bonds.end());

	    //Normalise the vectors.
	    BOOST_FOREACH(MagVec& vec, bonds)
	      {
		vec /= vec.nrm();
		
		double angle;
		if (vec[0] != 0)
		  angle = atan(vec[1]/vec[0]);
		else
		  if (vec[1] > 0)
		    angle = M_PI / 2;
		  else
		    angle = - M_PI / 2;
		
		if (vec[0] < 0) angle += M_PI;

		sum += std::exp(ComplexNum(0, 6 * angle));
	      }
	    
	    ++count;
	  }
      }
    _history.push_back(sum / ComplexNum(0, 6.0 * count));
  }

  void 
  OPOrientationalOrder::output(magnet::xml::XmlStream& XML)
  {
    XML << magnet::xml::tag("OrientationalOrder")
	<< magnet::xml::chardata();
    
    BOOST_FOREACH(const ComplexNum& val, _history)
      XML << "\n" << val.real() << " " << val.imag();
    
    XML << magnet::xml::endtag("OrientationalOrder");
  }
}
