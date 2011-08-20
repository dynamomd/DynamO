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


#include <dynamo/outputplugins/1partproperty/rijvijdirection.hpp>
#include <dynamo/dynamics/include.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/outputplugins/0partproperty/collMatrix.hpp>
#include <boost/foreach.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
  OPRijVij::OPRijVij(const dynamo::SimData* tmp, const magnet::xml::Node&):
    OutputPlugin(tmp, "RdotV")
  {}

  void 
  OPRijVij::initialise()
  {}

  void 
  OPRijVij::process2PED(mapdata& ref, const PairEventData& PDat)
  {
    Vector  rijnorm(PDat.rij / PDat.rij.nrm());
    Vector  vijnorm(PDat.vijold / PDat.vijold.nrm());

    double rvdot(rijnorm | vijnorm);

    for (size_t iDim(0); iDim < NDIM; ++iDim)
      {
	ref.rij[iDim].addVal(rijnorm[iDim]);
	ref.vij[iDim].addVal(vijnorm[iDim]);

	size_t id1(static_cast<size_t>((rijnorm[iDim] + 1.0) * 1000));
	size_t id2(static_cast<size_t>(-rvdot * 1000.0));

	++ref.rijcostheta[iDim].at(id1).first;
	ref.rijcostheta[iDim].at(id1).second += rvdot;

	++ref.costhetarij[iDim].at(id2).first;
	ref.costhetarij[iDim].at(id2).second += std::fabs(rijnorm[iDim]);


	id1 = static_cast<size_t>((rijnorm[iDim] + 1.0) * 100);
	id2 = static_cast<size_t>(-rvdot * 100.0);

	++ref.anglemapcount;
	++ref.anglemap[iDim][id1][id2];
      }
  }

  void 
  OPRijVij::eventUpdate(const IntEvent& iEvent, const PairEventData& pDat)
  {
  
    process2PED(rvdotacc[mapKey(iEvent.getType(), 
				getClassKey(iEvent))], pDat);
  }

  void 
  OPRijVij::eventUpdate(const GlobalEvent& globEvent, const NEventData& SDat)
  {
    BOOST_FOREACH(const PairEventData& pDat, SDat.L2partChanges)
      process2PED(rvdotacc[mapKey(globEvent.getType(), getClassKey(globEvent))],
		  pDat);
  }

  void 
  OPRijVij::eventUpdate(const LocalEvent& localEvent, const NEventData& SDat)
  {
    BOOST_FOREACH(const PairEventData& pDat, SDat.L2partChanges)
      process2PED(rvdotacc[mapKey(localEvent.getType(), getClassKey(localEvent))],
		  pDat);
  }

  void
  OPRijVij::eventUpdate(const System& sysEvent, const NEventData& SDat, const double&)
  {
    BOOST_FOREACH(const PairEventData& pDat, SDat.L2partChanges)
      process2PED(rvdotacc[mapKey(sysEvent.getType(), getClassKey(sysEvent))],
		  pDat);
  }

  void
  OPRijVij::output(magnet::xml::XmlStream &XML)
  {
    XML << magnet::xml::tag("RijVijComponents");
  
    typedef std::pair<const mapKey, mapdata> mappair;

    BOOST_FOREACH(const mappair& pair1, rvdotacc)
      {
	XML << magnet::xml::tag("Element")
	    << magnet::xml::attr("Type") 
	    << pair1.first.first
	    << magnet::xml::attr("EventName") 
	    << getName(pair1.first.second, Sim);

      
	for (size_t iDim(0); iDim < NDIM; ++iDim)
	  {
	    XML << magnet::xml::tag("Rij")
		<< magnet::xml::attr("dimension")
		<< iDim
		<< magnet::xml::chardata();

	    pair1.second.rij[iDim].outputHistogram(XML, 1.0);

	    XML << magnet::xml::endtag("Rij");
	  }

	for (size_t iDim(0); iDim < NDIM; ++iDim)
	  {
	    XML << magnet::xml::tag("Vij")
		<< magnet::xml::attr("dimension")
		<< iDim
		<< magnet::xml::chardata();

	    pair1.second.vij[iDim].outputHistogram(XML, 1.0);

	    XML << magnet::xml::endtag("Vij");
	  }

	for (size_t iDim(0); iDim < NDIM; ++iDim)
	  {
	    XML << magnet::xml::tag("RijVijvsRij")
		<< magnet::xml::attr("dimension")
		<< iDim
		<< magnet::xml::chardata();

	    for (size_t i(0); i < 2000; ++i)
	      XML << ((i - 1000.0) / 1000.0) << " "
		  << pair1.second.rijcostheta[iDim][i].second 
		/ pair1.second.rijcostheta[iDim][i].first
		  << "\n";

	    XML << magnet::xml::endtag("RijVijvsRij");
	  }

	for (size_t iDim(0); iDim < NDIM; ++iDim)
	  {
	    XML << magnet::xml::tag("RijvsRijVij")
		<< magnet::xml::attr("dimension")
		<< iDim
		<< magnet::xml::chardata();

	    for (size_t i(0); i < 1000; ++i)
	      XML << ( static_cast<double>(i) / -1000.0) << " "
		  << pair1.second.costhetarij[iDim][i].second 
		/ pair1.second.costhetarij[iDim][i].first
		  << "\n";

	    XML << magnet::xml::endtag("RijvsRijVij");
	  }

	for (size_t iDim(0); iDim < NDIM; ++iDim)
	  {
	    XML << magnet::xml::tag("XijRvdot")
		<< magnet::xml::attr("dimension")
		<< iDim
		<< magnet::xml::chardata();

	    for (size_t i1(0); i1 < 200; ++i1)
	      {	      
		for (size_t i2(0); i2 < 100; ++i2)
		  XML << ( (static_cast<double>(i1) - 100.0) / 100.0) << " "
		      << ( static_cast<double>(i2) / -100.0) << " "
		      << static_cast<double>(pair1.second.anglemap[iDim][i1][i2])
		    / static_cast<double>(pair1.second.anglemapcount)
		      << "\n";

		XML << "\n";
	      }
	  

	    XML << magnet::xml::endtag("XijRvdot");
	  }
      
	XML << magnet::xml::endtag("Element");
      }
  
  
    XML << magnet::xml::endtag("RijVijComponents");
  }
}
