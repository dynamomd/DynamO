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

#include "captures.hpp"
#include "../../simulation/particle.hpp"
#include "../../base/is_simdata.hpp"
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////

void 
ISingleCapture::initCaptureMap(const std::vector<Particle>& particleList)
{
  //If not loaded or invalidated
  if (noXmlLoad)
    {
      captureMap.clear();
      
      for (std::vector<Particle>::const_iterator iPtr
	     = particleList.begin();
	   iPtr != particleList.end(); iPtr++)
	for (std::vector<Particle>::const_iterator iPtr2 = iPtr + 1;
	     iPtr2 != particleList.end(); iPtr2++)
	  if (captureTest(*iPtr,*iPtr2))
	    addToCaptureMap(*iPtr, *iPtr2);
    }
}

void 
ISingleCapture::loadCaptureMap(const magnet::xml::Node& XML)
{
  if (XML.hasNode("CaptureMap"))
    {
      noXmlLoad = false;
      captureMap.clear();

      for (magnet::xml::Node node = XML.getNode("CaptureMap").fastGetNode("Pair");
	   node.valid(); ++node)
	captureMap.insert(std::pair<size_t, size_t>(node.getAttribute("ID1").as<size_t>(),
						    node.getAttribute("ID2").as<size_t>()));
    }
}

void 
ISingleCapture::outputCaptureMap(magnet::xml::XmlStream& XML) const 
{
  XML << magnet::xml::tag("CaptureMap");

  typedef std::pair<size_t, size_t> locpair;

  BOOST_FOREACH(const locpair& IDs, captureMap)
    XML << magnet::xml::tag("Pair")
	<< magnet::xml::attr("ID1") << IDs.first
	<< magnet::xml::attr("ID2") << IDs.second
	<< magnet::xml::endtag("Pair");
  
  XML << magnet::xml::endtag("CaptureMap");
}

void
ISingleCapture::addToCaptureMap(const Particle& p1, const Particle& p2) const
{
#ifdef DYNAMO_DEBUG
  if (p1.getID() == p2.getID())
    M_throw() << "Particle captured itself";

  if (p1.getID() < p2.getID())
    {
      if  (captureMap.count(std::pair<size_t, size_t>(p1.getID(), p2.getID())))
	M_throw() << "Insert found " << p1.getID()
		  << " and " << p2.getID() << " in the capture map";
    }      
  else
    if  (captureMap.count(std::pair<size_t, size_t>(p2.getID(), p1.getID())))
      M_throw() << "Insert found " << p2.getID() 
		<< " and " << p1.getID() << " in the capture map";
#endif 
  
  (p1.getID() < p2.getID())
    ? captureMap.insert(std::pair<size_t, size_t>(p1.getID(), p2.getID()))
    : captureMap.insert(std::pair<size_t, size_t>(p2.getID(), p1.getID()));
}

void 
ISingleCapture::removeFromCaptureMap(const Particle& p1, const Particle& p2) const
{
#ifdef DYNAMO_DEBUG
  if (p1.getID() == p2.getID())
    M_throw() << "Particle disassociated itself";

  if  (!((p1.getID() < p2.getID())
	 ? captureMap.erase(std::pair<size_t, size_t>(p1.getID(), p2.getID()))
	 : captureMap.erase(std::pair<size_t, size_t>(p2.getID(), p1.getID()))
	 ))
    M_throw() << "Erase did not find " << p2.getID() 
	      << " and " << p1.getID() << " in the capture map";    

#else

  (p1.getID() < p2.getID())
    ? captureMap.erase(std::pair<size_t, size_t>(p1.getID(), p2.getID()))
    : captureMap.erase(std::pair<size_t, size_t>(p2.getID(), p1.getID()));
  
#endif
}

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////

void 
IMultiCapture::initCaptureMap(const std::vector<Particle>& particleList)
{
  //If not loaded or invalidated
  if (noXmlLoad)
    {      
      captureMap.clear();
      
      for (std::vector<Particle>::const_iterator iPtr 
	     = particleList.begin();
	   iPtr != particleList.end(); iPtr++) 
	for (std::vector<Particle>::const_iterator iPtr2 = iPtr + 1;
	     iPtr2 != particleList.end(); iPtr2++)
	  {
	    int capval = captureTest(*iPtr,*iPtr2);
	    if (captureTest(*iPtr,*iPtr2))
	      captureMap[cMapKey(iPtr->getID(), iPtr2->getID())] 
		= capval; 
	  }
    }
}

void 
IMultiCapture::loadCaptureMap(const magnet::xml::Node& XML)
{
  if (XML.hasNode("CaptureMap"))
    {
      noXmlLoad = false;
      captureMap.clear();

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


bool 
IMultiCapture::isCaptured(const Particle& p1, const Particle& p2) const
{
#ifdef DYNAMO_DEBUG
  if (p1.getID() == p2.getID())
    M_throw() << "Particle is testing if it captured itself";
#endif 

  return captureMap.count(cMapKey(p1.getID(), p2.getID()));
}
