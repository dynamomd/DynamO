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

#include "captures.hpp"
#include <boost/foreach.hpp>
#include "../../base/is_exception.hpp"
#include "../../simulation/particle.hpp"
#include "../../base/is_simdata.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"

CICapture::CICapture(const DYNAMO::SimData* tmp,C2Range* nR): 
CInteraction(tmp,nR),
noXmlLoad(true),
captures(0)
{}

size_t
CICapture::getTotalCaptureCount() const
{
  return captures;
}

void 
CICapture::initCaptureMap()
{
  //If not loaded or invalidated
  if (noXmlLoad || (captureMap.size() != Sim->lN))
    {      
      if (captureMap.size() != Sim->lN)
	I_cout() << "Change in the number of particles invalidates the capture map!";
      
      captureMap.clear();
      captureMap.resize(Sim->lN);
      
      //Presume that the capture map was not loaded as its empty
      for (std::vector<CParticle>::const_iterator iPtr 
	     = Sim->vParticleList.begin();
	   iPtr != Sim->vParticleList.end(); iPtr++) 
	for (std::vector<CParticle>::const_iterator iPtr2 = iPtr + 1;
	     iPtr2 != Sim->vParticleList.end(); iPtr2++)
	  if (range->isInRange(*iPtr, *iPtr2))
	    if (captureTest(*iPtr,*iPtr2))	      
	      addToCaptureMap(*iPtr, *iPtr2); 
    }
}

void 
CICapture::loadCaptureMap(const XMLNode& XML)
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

      captures = 0;
      captureMap.clear();
      captureMap.resize(boost::lexical_cast<unsigned long>(subNode.getAttribute("Size")));

      int xml_iter = 0;
      long counter = subNode.nChildNode("Pair");
      for (long i = 0; i < counter; i++)
	{
	  browseNode = subNode.getChildNode("Pair",&xml_iter);
	  
	  captureMap.at(boost::lexical_cast<unsigned long>
			(browseNode.getAttribute("ID1")))
	    .insert(boost::lexical_cast<unsigned long>(browseNode.getAttribute("ID2")));
	  
	  ++captures;
	}
    }
}

void 
CICapture::outputCaptureMap(xmlw::XmlStream& XML) const 
{
  XML << xmlw::tag("CaptureMap") << xmlw::attr("Size") << Sim->lN;

  for (size_t i = 0; i < Sim->lN; ++i)
    if (!captureMap[i].empty())
      BOOST_FOREACH(const unsigned long& ID2, captureMap[i])
	XML << xmlw::tag("Pair")
	    << xmlw::attr("ID1") << i 
	    << xmlw::attr("ID2") << ID2
	    << xmlw::endtag("Pair");
  
  XML << xmlw::endtag("CaptureMap");
}


bool 
CICapture::isCaptured(const CParticle& p1, const CParticle& p2) const
{
#ifdef DYNAMO_DEBUG
  if (p1.getID() == p2.getID())
    D_throw() << "Particle is testing if it captured itself";
#endif 

  return (p1.getID() < p2.getID())
    ? captureMap[p1.getID()].count(p2.getID())
    : captureMap[p2.getID()].count(p1.getID());
}

void
CICapture::addToCaptureMap(const CParticle& p1, const CParticle& p2) const
{
#ifdef DYNAMO_DEBUG
  if (p1.getID() == p2.getID())
    D_throw() << "Particle captured itself";

  if (p1.getID() < p2.getID())
    {
      if  (captureMap[p1.getID()].count(p2.getID()))
	D_throw() << "Insert found " << p1.getID()
		  << " and " << p2.getID() << " in the capture map";
    }      
  else
    if  (captureMap[p2.getID()].count(p1.getID()))
      D_throw() << "Insert found " << p2.getID() 
		<< " and " << p1.getID() << " in the capture map";
#endif 
  
  (p1.getID() < p2.getID())
    ? captureMap[p1.getID()].insert(p2.getID())
    : captureMap[p2.getID()].insert(p1.getID());
  
    ++captures;
}

void 
CICapture::removeFromCaptureMap(const CParticle& p1, const CParticle& p2) const
{
#ifdef DYNAMO_DEBUG
  if (p1.getID() == p2.getID())
    D_throw() << "Particle disassociated itself";

  if  (!((p1.getID() < p2.getID())
	 ? captureMap[p1.getID()].erase(p2.getID())
	 : captureMap[p2.getID()].erase(p1.getID())))
    D_throw() << "Erase did not find " << p2.getID() 
	      << " and " << p1.getID() << " in the capture map";    

#else

  (p1.getID() < p2.getID())
    ? captureMap[p1.getID()].erase(p2.getID())
    : captureMap[p2.getID()].erase(p1.getID());
  
#endif
  
  --captures;
}
