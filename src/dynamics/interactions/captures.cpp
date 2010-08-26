/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
#include <boost/foreach.hpp>
#include "../../base/is_exception.hpp"
#include "../../simulation/particle.hpp"
#include "../../base/is_simdata.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"

ICapture::ICapture(DYNAMO::SimData* tmp,C2Range* nR): 
  Interaction(tmp,nR)
{}

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////

void 
ISingleCapture::initCaptureMap()
{
  //If not loaded or invalidated
  if (noXmlLoad)
    {      
      I_cout() << "Capture map reinitialising";
      
      captureMap.clear();
      
      for (std::vector<Particle>::const_iterator iPtr 
	     = Sim->particleList.begin();
	   iPtr != Sim->particleList.end(); iPtr++) 
	for (std::vector<Particle>::const_iterator iPtr2 = iPtr + 1;
	     iPtr2 != Sim->particleList.end(); iPtr2++)
	  if (range->isInRange(*iPtr, *iPtr2))
	    if (captureTest(*iPtr,*iPtr2))	      
	      addToCaptureMap(*iPtr, *iPtr2); 
    }
}

void 
ISingleCapture::loadCaptureMap(const XMLNode& XML)
{
  if (XML.nChildNode("CaptureMap"))
    {
      XMLNode browseNode, subNode;
      subNode = XML.getChildNode("CaptureMap");

      if (!subNode.isAttributeSet("Size"))
	{
	  I_cout() << "Could not find size in capture map";
	  noXmlLoad = true;
	  return;
	}

      noXmlLoad = false;

      captureMap.clear();

      int xml_iter = 0;
      long counter = subNode.nChildNode("Pair");
      for (long i = 0; i < counter; i++)
	{
	  browseNode = subNode.getChildNode("Pair",&xml_iter);
	  
	  captureMap.insert
	    (std::pair<size_t, size_t>
	     (boost::lexical_cast<size_t>(browseNode.getAttribute("ID1")),
	      boost::lexical_cast<size_t>(browseNode.getAttribute("ID2"))));
	}
    }
}

void 
ISingleCapture::outputCaptureMap(xml::XmlStream& XML) const 
{
  XML << xml::tag("CaptureMap") << xml::attr("Size") << Sim->N;

  typedef std::pair<size_t, size_t> locpair;

  BOOST_FOREACH(const locpair& IDs, captureMap)
    XML << xml::tag("Pair")
	<< xml::attr("ID1") << IDs.first
	<< xml::attr("ID2") << IDs.second
	<< xml::endtag("Pair");
  
  XML << xml::endtag("CaptureMap");
}

void
ISingleCapture::addToCaptureMap(const Particle& p1, const Particle& p2) const
{
#ifdef DYNAMO_DEBUG
  if (p1.getID() == p2.getID())
    D_throw() << "Particle captured itself";

  if (p1.getID() < p2.getID())
    {
      if  (captureMap.count(std::pair<size_t, size_t>(p1.getID(), p2.getID())))
	D_throw() << "Insert found " << p1.getID()
		  << " and " << p2.getID() << " in the capture map";
    }      
  else
    if  (captureMap.count(std::pair<size_t, size_t>(p2.getID(), p1.getID())))
      D_throw() << "Insert found " << p2.getID() 
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
    D_throw() << "Particle disassociated itself";

  if  (!((p1.getID() < p2.getID())
	 ? captureMap.erase(std::pair<size_t, size_t>(p1.getID(), p2.getID()))
	 : captureMap.erase(std::pair<size_t, size_t>(p2.getID(), p1.getID()))
	 ))
    D_throw() << "Erase did not find " << p2.getID() 
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
IMultiCapture::initCaptureMap()
{
  //If not loaded or invalidated
  if (noXmlLoad)
    {      
      I_cout() << "Capture map reinitialising";
      
      captureMap.clear();
      
      for (std::vector<Particle>::const_iterator iPtr 
	     = Sim->particleList.begin();
	   iPtr != Sim->particleList.end(); iPtr++) 
	for (std::vector<Particle>::const_iterator iPtr2 = iPtr + 1;
	     iPtr2 != Sim->particleList.end(); iPtr2++)
	  if (range->isInRange(*iPtr, *iPtr2))
	    {
	      int capval = captureTest(*iPtr,*iPtr2);
	      if (captureTest(*iPtr,*iPtr2))
		captureMap[cMapKey(iPtr->getID(), iPtr2->getID())] 
		  = capval; 
	    }
    }
}

void 
IMultiCapture::loadCaptureMap(const XMLNode& XML)
{
  if (XML.nChildNode("CaptureMap"))
    {
      XMLNode browseNode, subNode;
      subNode = XML.getChildNode("CaptureMap");

      if (!subNode.isAttributeSet("Size"))
	{
	  I_cout() << "Could not find size in capture map";
	  noXmlLoad = true;
	  return;
	}

      noXmlLoad = false;

      captureMap.clear();

      int xml_iter = 0;
      long counter = subNode.nChildNode("Pair");
      for (long i = 0; i < counter; i++)
	{
	  browseNode = subNode.getChildNode("Pair",&xml_iter);
	  
	  captureMap
	    [cMapKey
	     (boost::lexical_cast<size_t>(browseNode.getAttribute("ID1")),
	      boost::lexical_cast<size_t>(browseNode.getAttribute("ID2")))]
	    = boost::lexical_cast<size_t>(browseNode.getAttribute("val"));
	}
    }
}

void 
IMultiCapture::outputCaptureMap(xml::XmlStream& XML) const 
{
  XML << xml::tag("CaptureMap") << xml::attr("Size") << Sim->N;

  typedef std::pair<const cMapKey, int> locpair;

  BOOST_FOREACH(const locpair& IDs, captureMap)
    XML << xml::tag("Pair")
	<< xml::attr("ID1") << IDs.first.first
	<< xml::attr("ID2") << IDs.first.second
	<< xml::attr("val") << IDs.second
	<< xml::endtag("Pair");
  
  XML << xml::endtag("CaptureMap");
}


bool 
IMultiCapture::isCaptured(const Particle& p1, const Particle& p2) const
{
#ifdef DYNAMO_DEBUG
  if (p1.getID() == p2.getID())
    D_throw() << "Particle is testing if it captured itself";
#endif 

  return captureMap.count(cMapKey(p1.getID(), p2.getID()));
}
