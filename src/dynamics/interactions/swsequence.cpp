/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "swsequence.hpp"
#include <boost/lexical_cast.hpp>
#include <cmath>
#include "../../base/is_exception.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../BC/BC.hpp"
#include "../dynamics.hpp"
#include "../units/units.hpp"
#include "../globals/global.hpp"
#include "../../simulation/particle.hpp"
#include "../interactions/intEvent.hpp"
#include "../species/species.hpp"
#include "../2particleEventData.hpp"
#include "../liouvillean/liouvillean.hpp"
#include "../../base/is_simdata.hpp"
#include "../../schedulers/scheduler.hpp"
#include "../NparticleEventData.hpp"
#include <iomanip>

CISWSequence::CISWSequence(DYNAMO::SimData* tmp, Iflt nd, Iflt nl,
			   Iflt ne, std::vector<size_t> seq, C2Range* nR):
  CISingleCapture(tmp,nR),
  diameter(nd),d2(nd*nd),lambda(nl),ld2(nd*nd*nl*nl),e(ne), sequence(seq) 
{
  std::set<size_t> letters;

  BOOST_FOREACH(const size_t& id, seq)
    if (letters.find(id) == letters.end())
      //Count the letter
      letters.insert(id);
  
  alphabet.resize(letters.size());
  
  BOOST_FOREACH(std::vector<Iflt>& vec, alphabet)
    vec.resize(letters.size(), 0.0);
}

CISWSequence::CISWSequence(const XMLNode& XML, DYNAMO::SimData* tmp):
  CISingleCapture(tmp, NULL) //A temporary value!
{
  operator<<(XML);
}

Iflt 
CISWSequence::getColourFraction(const CParticle& part) const
{
  return (sequence[part.getID() % sequence.size()] + 0.5) / static_cast<Iflt>(alphabet.size());   
}

void 
CISWSequence::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") << "SquareWellSeq"
      << xmlw::attr("Diameter") 
      << diameter / Sim->dynamics.units().unitLength() 
      << xmlw::attr("Elasticity") << e
      << xmlw::attr("Lambda") << lambda
      << xmlw::attr("Name") << intName
      << range;

  XML << xmlw::tag("Sequence");
  
  for (size_t i = 0; i < sequence.size(); ++i)
    XML << xmlw::tag("Element")
	<< xmlw::attr("seqID") << i
	<< xmlw::attr("Letter") << sequence[i]
	<< xmlw::endtag("Element");

  XML << xmlw::endtag("Sequence")
      << xmlw::tag("Alphabet");

  for (size_t i = 0; i < alphabet.size(); ++i)
    for (size_t j = i; j < alphabet[i].size(); ++j)
      XML << xmlw::tag("Word")
	  << xmlw::attr("Letter1") << i
	  << xmlw::attr("Letter2") << j
	  << xmlw::attr("Depth") << alphabet[i][j]
	  << xmlw::endtag("Word");

  XML << xmlw::endtag("Alphabet");

  
  CISingleCapture::outputCaptureMap(XML);  
}

void 
CISWSequence::operator<<(const XMLNode& XML)
{
  if (strcmp(XML.getAttribute("Type"),"SquareWellSeq"))
    D_throw() << "Attempting to load SquareWell from non SquareWell entry";
  
  range.set_ptr(C2Range::loadClass(XML,Sim));
  
  try { 
    diameter = Sim->dynamics.units().unitLength() 
      * boost::lexical_cast<Iflt>(XML.getAttribute("Diameter"));
    
    e = boost::lexical_cast<Iflt>(XML.getAttribute("Elasticity"));
        
    lambda = boost::lexical_cast<Iflt>(XML.getAttribute("Lambda"));
    
    d2 = diameter * diameter;
    
    ld2 = d2 * lambda * lambda;
    
    intName = XML.getAttribute("Name");

    CISingleCapture::loadCaptureMap(XML);
  
    //Load the sequence
    XMLNode subNode = XML.getChildNode("Sequence");

    sequence.clear();
    sequence.resize(subNode.nChildNode("Element"), 0);

    int xml_iter = 0;
    std::set<size_t> letters;

    for (size_t i = 0; i < sequence.size(); i++)
      {
	XMLNode browseNode = subNode.getChildNode("Element", &xml_iter);

	//Check if the letter is there before
	size_t letter = boost::lexical_cast<size_t>(browseNode.getAttribute("Letter"));
	if (letters.find(letter) == letters.end())
	  //Count the letter
	  letters.insert(letter);

	sequence[boost::lexical_cast<size_t>(browseNode.getAttribute("seqID"))] 
	  = letter;
      }

    //Initialise all the well depths to 1.0
    alphabet.resize(letters.size());

    BOOST_FOREACH(std::vector<Iflt>& vec, alphabet)
      vec.resize(letters.size(), 0.0);

    //Load the sequence
    subNode = XML.getChildNode("Alphabet");
    size_t count = subNode.nChildNode("Word");

    xml_iter = 0;
    for (size_t i = 0; i < count; ++i)
      {
	XMLNode browseNode = subNode.getChildNode("Word", &xml_iter);
	
	alphabet
	  .at(boost::lexical_cast<size_t>(browseNode.getAttribute("Letter1")))
	  .at(boost::lexical_cast<size_t>(browseNode.getAttribute("Letter2")))
	  = boost::lexical_cast<Iflt>(browseNode.getAttribute("Depth"));

	alphabet
	  .at(boost::lexical_cast<size_t>(browseNode.getAttribute("Letter2")))
	  .at(boost::lexical_cast<size_t>(browseNode.getAttribute("Letter1")))
	  = boost::lexical_cast<Iflt>(browseNode.getAttribute("Depth"));
      }

  }
  catch (boost::bad_lexical_cast &)
    {
      D_throw() << "Failed a lexical cast in CISWSequence";
    }
}

CInteraction* 
CISWSequence::Clone() const 
{ return new CISWSequence(*this); }

Iflt 
CISWSequence::getInternalEnergy() const 
{ 
  //Once the capture maps are loaded just iterate through that determining energies
  Iflt Energy = 0.0;
  typedef std::pair<size_t, size_t> locpair;

  BOOST_FOREACH(const locpair& IDs, captureMap)
    Energy += alphabet
    [sequence[IDs.first % sequence.size()]]
    [sequence[IDs.second % sequence.size()]];
  
  return -Energy; 
}

Iflt 
CISWSequence::hardCoreDiam() const 
{ return diameter; }

Iflt 
CISWSequence::maxIntDist() const 
{ return diameter * lambda; }

void 
CISWSequence::rescaleLengths(Iflt scale) 
{ 
  diameter += scale * diameter; 
  d2 = diameter * diameter;
  ld2 = diameter * diameter * lambda * lambda;
}

void 
CISWSequence::initialise(size_t nID)
{
  ID = nID;
  CISingleCapture::initCaptureMap();
}

bool 
CISWSequence::captureTest(const CParticle& p1, const CParticle& p2) const
{
  Vector  rij = p1.getPosition() - p2.getPosition();
  Sim->dynamics.BCs().applyBC(rij);
  
  if ((rij.nrm2() <= ld2) && (rij.nrm2() >= d2))
    return true;
  
  return false;
}

CIntEvent 
CISWSequence::getEvent(const CParticle &p1, 
		       const CParticle &p2) const 
{    
#ifdef DYNAMO_DEBUG
  if (!Sim->dynamics.getLiouvillean().isUpToDate(p1))
    D_throw() << "Particle 1 is not up to date";
  
  if (!Sim->dynamics.getLiouvillean().isUpToDate(p2))
    D_throw() << "Particle 2 is not up to date";

  if (p1 == p2)
    D_throw() << "You shouldn't pass p1==p2 events to the interactions!";
#endif 

  CPDData colldat(*Sim, p1, p2);
  
  if (isCaptured(p1, p2)) 
    {
      if (Sim->dynamics.getLiouvillean().SphereSphereInRoot(colldat, d2)) 
	{
#ifdef DYNAMO_OverlapTesting
	  //Check that there is no overlap 
	  if (Sim->dynamics.getLiouvillean().sphereOverlap(colldat, d2))
	    D_throw() << "Overlapping particles found" 
		      << ", particle1 " << p1.getID() 
		      << ", particle2 " 
		      << p2.getID() << "\nOverlap = " << (sqrt(colldat.r2) - sqrt(d2))/Sim->dynamics.units().unitLength();
#endif	  
	  return CIntEvent(p1, p2, colldat.dt, CORE, *this);
	}
      else
	if (Sim->dynamics.getLiouvillean().SphereSphereOutRoot(colldat, ld2))
	  return CIntEvent(p1, p2, colldat.dt, WELL_OUT, *this);
    }
  else if (Sim->dynamics.getLiouvillean().SphereSphereInRoot(colldat, ld2)) 
    {
#ifdef DYNAMO_OverlapTesting
      if (Sim->dynamics.getLiouvillean().sphereOverlap(colldat,ld2))
	{
	  if (Sim->dynamics.getLiouvillean().sphereOverlap(colldat,d2))
	    D_throw() << "Overlapping cores (but not registerd as captured) particles found in square well" 
		      << "\nparticle1 " << p1.getID() << ", particle2 " 
		      << p2.getID() << "\nOverlap = " << (sqrt(colldat.r2) - sqrt(d2))/Sim->dynamics.units().unitLength();
	  else
	    D_throw() << "Overlapping wells (but not registerd as captured) particles found" 
		      << "\nparticle1 " << p1.getID() << ", particle2 " 
		      << p2.getID() << "\nOverlap = " << (sqrt(colldat.r2) - sqrt(ld2))/Sim->dynamics.units().unitLength();
	  
	}
#endif
      return CIntEvent(p1, p2, colldat.dt, WELL_IN, *this);
    }

  return CIntEvent(p1, p2, HUGE_VAL, NONE, *this);
}

void
CISWSequence::runEvent(const CParticle& p1,
		       const CParticle& p2,
		       const CIntEvent& iEvent) const
{  
  ++Sim->lNColl;

  switch (iEvent.getType())
    {
    case CORE:
      {
	C2ParticleData retVal(Sim->dynamics.getLiouvillean().SmoothSpheresColl(iEvent, e, d2, CORE));
	Sim->signalParticleUpdate(retVal);
	
	Sim->ptrScheduler->fullUpdate(p1, p2);
	
	BOOST_FOREACH(smrtPlugPtr<OutputPlugin> & Ptr, Sim->outputPlugins)
	  Ptr->eventUpdate(iEvent, retVal);

	break;
      }
    case WELL_IN:
      {
	C2ParticleData retVal(Sim->dynamics.getLiouvillean()
			      .SphereWellEvent
			      (iEvent, alphabet
			       [sequence[p1.getID() % sequence.size()]]
			       [sequence[p2.getID() % sequence.size()]], 
			       ld2));
	
	if (retVal.getType() != BOUNCE)
	  addToCaptureMap(p1, p2);      

	Sim->signalParticleUpdate(retVal);

	Sim->ptrScheduler->fullUpdate(p1, p2);
	
	BOOST_FOREACH(smrtPlugPtr<OutputPlugin> & Ptr, Sim->outputPlugins)
	  Ptr->eventUpdate(iEvent, retVal);

	break;
      }
    case WELL_OUT:
      {
	C2ParticleData retVal(Sim->dynamics.getLiouvillean()
			      .SphereWellEvent
			      (iEvent, -alphabet
			       [sequence[p1.getID() % sequence.size()]]
			       [sequence[p2.getID() % sequence.size()]], 
			       ld2));
	
	if (retVal.getType() != BOUNCE)
	  removeFromCaptureMap(p1, p2);
	
	Sim->signalParticleUpdate(retVal);

	Sim->ptrScheduler->fullUpdate(p1, p2);
	
	BOOST_FOREACH(smrtPlugPtr<OutputPlugin> & Ptr, Sim->outputPlugins)
	  Ptr->eventUpdate(iEvent, retVal);

	break;
      }
    default:
      D_throw() << "Unknown collision type";
    }
}

void
CISWSequence::checkOverlaps(const CParticle& part1, const CParticle& part2) const
{
  Vector  rij = part1.getPosition() - part2.getPosition();
  Sim->dynamics.BCs().applyBC(rij);
  Iflt r2 = rij.nrm2();

  if (isCaptured(part1, part2))
    {
      if (r2 < d2)
	I_cerr() << std::setprecision(std::numeric_limits<float>::digits10)
		 << "Possible captured overlap occured in diagnostics\n ID1=" << part1.getID() 
		 << ", ID2=" << part2.getID() << "\nR_ij^2=" 
		 << r2 / pow(Sim->dynamics.units().unitLength(),2)
		 << "\nd^2=" 
		 << d2 / pow(Sim->dynamics.units().unitLength(),2);

      if (r2 > lambda * lambda * d2)
	I_cerr() << std::setprecision(std::numeric_limits<float>::digits10)
		 << "Possible escaped captured pair in diagnostics\n ID1=" << part1.getID() 
		 << ", ID2=" << part2.getID() << "\nR_ij^2=" 
		 << r2 / pow(Sim->dynamics.units().unitLength(),2)
		 << "\n(lambda * d)^2=" 
		 << ld2 / pow(Sim->dynamics.units().unitLength(),2);
    }
  else 
    if (r2 < ld2)
      I_cerr() << std::setprecision(std::numeric_limits<float>::digits10)
	       << "Possible missed captured pair in diagnostics\n ID1=" << part1.getID() 
	       << ", ID2=" << part2.getID() << "\nR_ij^2=" 
	       << r2 / pow(Sim->dynamics.units().unitLength(),2)
	       << "\n(lambda * d)^2=" 
	       << ld2 / pow(Sim->dynamics.units().unitLength(),2);
}

void 
CISWSequence::write_povray_desc(const DYNAMO::RGB& rgb, 
				const size_t& specID, 
				std::ostream& os) const
{
  DYNAMO::ColorMap<size_t> seqmap(0, alphabet.size() * Sim->dynamics.getSpecies().size() - 1);
  
  for (size_t i(0); i < alphabet.size(); ++i)
    {
      DYNAMO::RGB col(seqmap.getColor(i * Sim->dynamics.getSpecies().size() + specID));

      os << "#declare intrep" << ID << "center"<< i << " = " 
	 << "sphere {\n <0,0,0> " 
	 << diameter / 2.0
	 << "\n texture { pigment { color rgb<" << col.R << "," << col.G
	 << "," << col.B << "> }}\nfinish { phong 0.9 phong_size 60 }\n}\n"
	 << "#declare intrep" << ID << "seqwell" << i
	 << " = sphere {\n <0,0,0> " << diameter * lambda * 0.5
	 << "\n texture { pigment { color rgbt <1,1,1,0.9> }}\n}\n";
    }

  BOOST_FOREACH(const size_t& part, *(Sim->dynamics.getSpecies()[specID]->getRange()))
    {
      Vector  pos(Sim->vParticleList[part].getPosition());
      Sim->dynamics.BCs().applyBC(pos);
      
      os << "object {\n intrep" << ID << "center"<< sequence[part % sequence.size()] << "\n translate < "
	 << pos[0] << ", " << pos[1] << ", " << pos[2] << ">\n}\n";
    }

  os << "merge {\n";
  BOOST_FOREACH(const size_t& part, *(Sim->dynamics.getSpecies()[specID]->getRange()))
    {
      Vector  pos(Sim->vParticleList[part].getPosition());
      Sim->dynamics.BCs().applyBC(pos);
      
      os << "object {\n intrep" << ID << "seqwell" << sequence[part % sequence.size()] << "\n translate < "
	 << pos[0] << ", " << pos[1] << ", " << pos[2] << ">\n}\n";
    }
  
  os << "}\n";

}
