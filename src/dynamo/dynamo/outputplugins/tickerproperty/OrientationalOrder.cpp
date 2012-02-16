/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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
#include <dynamo/dynamics/globals/neighbourList.hpp>
#include <dynamo/dynamics/units/units.hpp>
#include <dynamo/dynamics/BC/BC.hpp>
#include <boost/foreach.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <fstream>
#include <cmath>
#include <limits>

namespace dynamo {
  OPOrientationalOrder::OPOrientationalOrder(const dynamo::SimData* tmp, const magnet::xml::Node& XML):
    OPTicker(tmp,"OrientationalOrder"), 
    _axis(1,0,0),
    _rg(1),
    _nblistID(std::numeric_limits<size_t>::max())
  {
    operator<<(XML);
  }

  void 
  OPOrientationalOrder::operator<<(const magnet::xml::Node& XML)
  {
    try
      {
	if (XML.hasAttribute("CutOffR"))
	  _rg = XML.getAttribute("CutOffR").as<double>();
      }
    catch (boost::bad_lexical_cast &)
      {
	M_throw() << "Failed a lexical cast in OPOrientationalOrder";
      }

    _rg *= Sim->dynamics.units().unitLength();

    dout << "Cut off radius set to " 
	 << _rg / Sim->dynamics.units().unitLength() << std::endl;
  }


  void 
  OPOrientationalOrder::initialise() 
  { 
    double smallestlength = HUGE_VAL;

    BOOST_FOREACH(const shared_ptr<Global>& pGlob, Sim->dynamics.getGlobals())
      if (std::tr1::dynamic_pointer_cast<GNeighbourList>(pGlob))
	{
	  const double l(static_cast<const GNeighbourList*>(pGlob.get())
			 ->getMaxSupportedInteractionLength());
	  if ((l >= _rg) && (l < smallestlength))
	    {
	      //this neighbourlist is better suited
	      smallestlength = l;
	      _nblistID = pGlob->getID();
	    }
	}

    if (_nblistID == std::numeric_limits<size_t>::max())
      M_throw() << "There is not a suitable neighbourlist for the cut-off radius selected."
	"\nR_g = " << _rg / Sim->dynamics.units().unitLength();

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
    BOOST_FOREACH(const Particle& part, Sim->particleList)
      {
	Neighbours nbs;
	
	static_cast<const GNeighbourList*>
	  (Sim->dynamics.getGlobals()[_nblistID].get())
	  ->getParticleNeighbourhood
	  (part, magnet::function::MakeDelegate
	   (&nbs, &Neighbours::addNeighbour));
	
	if (nbs._neighbours.size() >= 6)
	  {
	    std::vector<MagVec> bonds;
	    BOOST_FOREACH(const size_t& id2, nbs._neighbours)
	      {
		Vector bond = Sim->particleList[id2].getPosition() - part.getPosition();
		Sim->dynamics.BCs().applyBC(bond);
		bonds.push_back(bond);
	      }
	    std::sort(bonds.begin(), bonds.end());

	    if (bonds.size() > 6)
	      bonds.erase(bonds.begin()+6, bonds.end());

	    //Normalise the vectors.
	    BOOST_FOREACH(MagVec& vec, bonds)
	      {
		vec /= vec.nrm();
		ComplexNum ang1(vec[0], vec[1]);
		ComplexNum ang2 = ang1 * ang1;
		sum += ang2 * ang2 * ang2;
	      }
	    
	    ++count;
	  }
      }
    sum /= 6.0 * count;

    _history.push_back(sum);
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
