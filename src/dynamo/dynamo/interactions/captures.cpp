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

#include <dynamo/interactions/captures.hpp>
#include <dynamo/particle.hpp>
#include <dynamo/schedulers/scheduler.hpp>
#include <dynamo/simulation.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>

namespace dynamo {
  void 
  ISingleCapture::testAddToCaptureMap(const Particle& p1, const size_t& p2) const
  {
    if (captureTest(p1, Sim->particleList[p2]))
      addToCaptureMap(p1, Sim->particleList[p2]);
  }   


  void 
  ICapture::initCaptureMap()
  {
    //If not loaded or invalidated
    if (noXmlLoad)
      {      
	clear();

	for (std::vector<Particle>::const_iterator iPtr1 = Sim->particleList.begin();
	     iPtr1 != Sim->particleList.end(); iPtr1++)
	  for (std::vector<Particle>::const_iterator iPtr2 = iPtr1+1;
	       iPtr2 != Sim->particleList.end(); iPtr2++)
	    testAddToCaptureMap(*iPtr1, iPtr2->getID());
      }
  }

  void 
  ISingleCapture::loadCaptureMap(const magnet::xml::Node& XML)
  {
    if (XML.hasNode("CaptureMap"))
      {
	noXmlLoad = false;
	clear();

	for (magnet::xml::Node node = XML.getNode("CaptureMap").fastGetNode("Pair");
	     node.valid(); ++node)
	  captureMap.insert(cMapKey(node.getAttribute("ID1").as<size_t>(),
				    node.getAttribute("ID2").as<size_t>()));
      }
  }

  void 
  ISingleCapture::outputCaptureMap(magnet::xml::XmlStream& XML) const 
  {
    XML << magnet::xml::tag("CaptureMap");

    BOOST_FOREACH(const cMapKey& IDs, captureMap)
      XML << magnet::xml::tag("Pair")
	  << magnet::xml::attr("ID1") << IDs.first
	  << magnet::xml::attr("ID2") << IDs.second
	  << magnet::xml::endtag("Pair");
  
    XML << magnet::xml::endtag("CaptureMap");
  }

  //////////////////////////////////////////////////////
  //////////////////////////////////////////////////////

  void 
  IMultiCapture::testAddToCaptureMap(const Particle& p1, const size_t& p2) const
  {
    int capval = captureTest(p1, Sim->particleList[p2]);
    if (capval) captureMap[cMapKey(p1.getID(), p2)] = capval; 
  }

  void 
  IMultiCapture::loadCaptureMap(const magnet::xml::Node& XML)
  {
    if (XML.hasNode("CaptureMap"))
      {
	noXmlLoad = false;
	clear();

	for (magnet::xml::Node node = XML.getNode("CaptureMap").fastGetNode("Pair");
	     node.valid(); ++node)
	  captureMap[cMapKey(node.getAttribute("ID1").as<size_t>(),
			     node.getAttribute("ID2").as<size_t>())]
	    = node.getAttribute("val").as<size_t>();
      }
  }

  void 
  IMultiCapture::outputCaptureMap(magnet::xml::XmlStream& XML) const 
  {
    XML << magnet::xml::tag("CaptureMap");

    typedef std::pair<const cMapKey, int> locpair;

    BOOST_FOREACH(const locpair& IDs, captureMap)
      XML << magnet::xml::tag("Pair")
	  << magnet::xml::attr("ID1") << IDs.first.first
	  << magnet::xml::attr("ID2") << IDs.first.second
	  << magnet::xml::attr("val") << IDs.second
	  << magnet::xml::endtag("Pair");
  
    XML << magnet::xml::endtag("CaptureMap");
  }
}
