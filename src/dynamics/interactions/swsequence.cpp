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

#include "swsequence.hpp"
#include <boost/lexical_cast.hpp>
#include <cmath>
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

ISWSequence::ISWSequence(DYNAMO::SimData* tmp, double nd, double nl,
			   double ne, std::vector<size_t> seq, C2Range* nR):
  ISingleCapture(tmp,nR),
  diameter(nd),d2(nd*nd),lambda(nl),ld2(nd*nd*nl*nl),e(ne), sequence(seq) 
{
  std::set<size_t> letters;

  BOOST_FOREACH(const size_t& id, seq)
    if (letters.find(id) == letters.end())
      //Count the letter
      letters.insert(id);
  
  alphabet.resize(letters.size());
  
  BOOST_FOREACH(std::vector<double>& vec, alphabet)
    vec.resize(letters.size(), 0.0);
}

ISWSequence::ISWSequence(const XMLNode& XML, DYNAMO::SimData* tmp):
  ISingleCapture(tmp, NULL) //A temporary value!
{
  operator<<(XML);
}

double 
ISWSequence::getColourFraction(const Particle& part) const
{
  return (sequence[part.getID() % sequence.size()] + 0.5) / static_cast<double>(alphabet.size());   
}

void 
ISWSequence::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Type") << "SquareWellSeq"
      << xml::attr("Diameter") 
      << diameter / Sim->dynamics.units().unitLength() 
      << xml::attr("Elasticity") << e
      << xml::attr("Lambda") << lambda
      << xml::attr("Name") << intName
      << range;

  XML << xml::tag("Sequence");
  
  for (size_t i = 0; i < sequence.size(); ++i)
    XML << xml::tag("Element")
	<< xml::attr("seqID") << i
	<< xml::attr("Letter") << sequence[i]
	<< xml::endtag("Element");

  XML << xml::endtag("Sequence")
      << xml::tag("Alphabet");

  for (size_t i = 0; i < alphabet.size(); ++i)
    for (size_t j = i; j < alphabet[i].size(); ++j)
      XML << xml::tag("Word")
	  << xml::attr("Letter1") << i
	  << xml::attr("Letter2") << j
	  << xml::attr("Depth") << alphabet[i][j]
	  << xml::endtag("Word");

  XML << xml::endtag("Alphabet");

  
  ISingleCapture::outputCaptureMap(XML);  
}

void 
ISWSequence::operator<<(const XMLNode& XML)
{
  if (strcmp(XML.getAttribute("Type"),"SquareWellSeq"))
    M_throw() << "Attempting to load SquareWell from non SquareWell entry";
  
  range.set_ptr(C2Range::loadClass(XML,Sim));
  
  try { 
    diameter = Sim->dynamics.units().unitLength() 
      * boost::lexical_cast<double>(XML.getAttribute("Diameter"));
    
    e = boost::lexical_cast<double>(XML.getAttribute("Elasticity"));
        
    lambda = boost::lexical_cast<double>(XML.getAttribute("Lambda"));
    
    d2 = diameter * diameter;
    
    ld2 = d2 * lambda * lambda;
    
    intName = XML.getAttribute("Name");

    ISingleCapture::loadCaptureMap(XML);
  
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

    BOOST_FOREACH(std::vector<double>& vec, alphabet)
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
	  = boost::lexical_cast<double>(browseNode.getAttribute("Depth"));

	alphabet
	  .at(boost::lexical_cast<size_t>(browseNode.getAttribute("Letter2")))
	  .at(boost::lexical_cast<size_t>(browseNode.getAttribute("Letter1")))
	  = boost::lexical_cast<double>(browseNode.getAttribute("Depth"));
      }

  }
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in CISWSequence";
    }
}

Interaction* 
ISWSequence::Clone() const 
{ return new ISWSequence(*this); }

double 
ISWSequence::getInternalEnergy() const 
{ 
  //Once the capture maps are loaded just iterate through that determining energies
  double Energy = 0.0;
  typedef std::pair<size_t, size_t> locpair;

  BOOST_FOREACH(const locpair& IDs, captureMap)
    Energy += alphabet
    [sequence[IDs.first % sequence.size()]]
    [sequence[IDs.second % sequence.size()]];
  
  return -Energy; 
}

double 
ISWSequence::hardCoreDiam() const 
{ return diameter; }

double 
ISWSequence::maxIntDist() const 
{ return diameter * lambda; }

void 
ISWSequence::rescaleLengths(double scale) 
{ 
  diameter += scale * diameter; 
  d2 = diameter * diameter;
  ld2 = diameter * diameter * lambda * lambda;
}

void 
ISWSequence::initialise(size_t nID)
{
  ID = nID;
  ISingleCapture::initCaptureMap();
}

bool 
ISWSequence::captureTest(const Particle& p1, const Particle& p2) const
{
  Vector  rij = p1.getPosition() - p2.getPosition();
  Sim->dynamics.BCs().applyBC(rij);

#ifdef DYNAMO_DEBUG
  if (rij.nrm2() >= d2)
    I_cerr() << "Warning! Two particles might be overlapping"
	     << "\nrij^2 = " << (rij | rij)
	     << "\nd^2 = " << d2;
#endif

  return (rij.nrm2() <= ld2);
}

IntEvent 
ISWSequence::getEvent(const Particle &p1, 
		       const Particle &p2) const 
{    
#ifdef DYNAMO_DEBUG
  if (!Sim->dynamics.getLiouvillean().isUpToDate(p1))
    M_throw() << "Particle 1 is not up to date";
  
  if (!Sim->dynamics.getLiouvillean().isUpToDate(p2))
    M_throw() << "Particle 2 is not up to date";

  if (p1 == p2)
    M_throw() << "You shouldn't pass p1==p2 events to the interactions!";
#endif 

  CPDData colldat(*Sim, p1, p2);
  
  if (isCaptured(p1, p2)) 
    {
      if (Sim->dynamics.getLiouvillean().SphereSphereInRoot(colldat, d2)) 
	{
#ifdef DYNAMO_OverlapTesting
	  //Check that there is no overlap 
	  if (Sim->dynamics.getLiouvillean().sphereOverlap(colldat, d2))
	    M_throw() << "Overlapping particles found" 
		      << ", particle1 " << p1.getID() 
		      << ", particle2 " 
		      << p2.getID() << "\nOverlap = " << (sqrt(colldat.r2) - sqrt(d2))/Sim->dynamics.units().unitLength();
#endif	  
	  return IntEvent(p1, p2, colldat.dt, CORE, *this);
	}
      else
	if (Sim->dynamics.getLiouvillean().SphereSphereOutRoot(colldat, ld2))
	  return IntEvent(p1, p2, colldat.dt, WELL_OUT, *this);
    }
  else if (Sim->dynamics.getLiouvillean().SphereSphereInRoot(colldat, ld2)) 
    {
#ifdef DYNAMO_OverlapTesting
      if (Sim->dynamics.getLiouvillean().sphereOverlap(colldat,ld2))
	{
	  if (Sim->dynamics.getLiouvillean().sphereOverlap(colldat,d2))
	    M_throw() << "Overlapping cores (but not registerd as captured) particles found in square well" 
		      << "\nparticle1 " << p1.getID() << ", particle2 " 
		      << p2.getID() << "\nOverlap = " << (sqrt(colldat.r2) - sqrt(d2))/Sim->dynamics.units().unitLength();
	  else
	    M_throw() << "Overlapping wells (but not registerd as captured) particles found" 
		      << "\nparticle1 " << p1.getID() << ", particle2 " 
		      << p2.getID() << "\nOverlap = " << (sqrt(colldat.r2) - sqrt(ld2))/Sim->dynamics.units().unitLength();
	  
	}
#endif
      return IntEvent(p1, p2, colldat.dt, WELL_IN, *this);
    }

  return IntEvent(p1, p2, HUGE_VAL, NONE, *this);
}

void
ISWSequence::runEvent(const Particle& p1,
		       const Particle& p2,
		       const IntEvent& iEvent) const
{  
  ++Sim->eventCount;

  switch (iEvent.getType())
    {
    case CORE:
      {
	PairEventData retVal(Sim->dynamics.getLiouvillean().SmoothSpheresColl(iEvent, e, d2, CORE));
	Sim->signalParticleUpdate(retVal);
	
	Sim->ptrScheduler->fullUpdate(p1, p2);
	
	BOOST_FOREACH(magnet::ClonePtr<OutputPlugin> & Ptr, Sim->outputPlugins)
	  Ptr->eventUpdate(iEvent, retVal);

	break;
      }
    case WELL_IN:
      {
	PairEventData retVal(Sim->dynamics.getLiouvillean()
			      .SphereWellEvent
			      (iEvent, alphabet
			       [sequence[p1.getID() % sequence.size()]]
			       [sequence[p2.getID() % sequence.size()]], 
			       ld2));
	
	if (retVal.getType() != BOUNCE)
	  addToCaptureMap(p1, p2);      

	Sim->signalParticleUpdate(retVal);

	Sim->ptrScheduler->fullUpdate(p1, p2);
	
	BOOST_FOREACH(magnet::ClonePtr<OutputPlugin> & Ptr, Sim->outputPlugins)
	  Ptr->eventUpdate(iEvent, retVal);

	break;
      }
    case WELL_OUT:
      {
	PairEventData retVal(Sim->dynamics.getLiouvillean()
			      .SphereWellEvent
			      (iEvent, -alphabet
			       [sequence[p1.getID() % sequence.size()]]
			       [sequence[p2.getID() % sequence.size()]], 
			       ld2));
	
	if (retVal.getType() != BOUNCE)
	  removeFromCaptureMap(p1, p2);
	
	Sim->signalParticleUpdate(retVal);

	Sim->ptrScheduler->fullUpdate(p1, p2);
	
	BOOST_FOREACH(magnet::ClonePtr<OutputPlugin> & Ptr, Sim->outputPlugins)
	  Ptr->eventUpdate(iEvent, retVal);

	break;
      }
    default:
      M_throw() << "Unknown collision type";
    }
}

void
ISWSequence::checkOverlaps(const Particle& part1, const Particle& part2) const
{
  Vector  rij = part1.getPosition() - part2.getPosition();
  Sim->dynamics.BCs().applyBC(rij);
  double r2 = rij.nrm2();

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
    {
      if (r2 < d2)
	I_cerr() << "Particles overlapping cores without even being captured."
		 << "\nProbably a bad initial configuration."
		 << std::setprecision(std::numeric_limits<float>::digits10)
		 << "\n ID1=" 
		 << part1.getID() 
		 << ", ID2=" << part2.getID() << "\nR_ij^2=" 
		 << r2 / pow(Sim->dynamics.units().unitLength(),2)
		 << "\nd^2=" 
		 << d2 / pow(Sim->dynamics.units().unitLength(),2);
      if (r2 < ld2)
	I_cerr() << std::setprecision(std::numeric_limits<float>::digits10)
		 << "Possible missed captured pair in diagnostics\n ID1=" 
		 << part1.getID() 
		 << ", ID2=" << part2.getID() << "\nR_ij^2=" 
		 << r2 / pow(Sim->dynamics.units().unitLength(),2)
		 << "\n(lambda * d)^2=" 
		 << ld2 / pow(Sim->dynamics.units().unitLength(),2);
    }
}

void 
ISWSequence::write_povray_desc(const DYNAMO::RGB& rgb, 
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
      Vector  pos(Sim->particleList[part].getPosition());
      Sim->dynamics.BCs().applyBC(pos);
      
      os << "object {\n intrep" << ID << "center"<< sequence[part % sequence.size()] << "\n translate < "
	 << pos[0] << ", " << pos[1] << ", " << pos[2] << ">\n}\n";
    }

  os << "merge {\n";
  BOOST_FOREACH(const size_t& part, *(Sim->dynamics.getSpecies()[specID]->getRange()))
    {
      Vector  pos(Sim->particleList[part].getPosition());
      Sim->dynamics.BCs().applyBC(pos);
      
      os << "object {\n intrep" << ID << "seqwell" << sequence[part % sequence.size()] << "\n translate < "
	 << pos[0] << ", " << pos[1] << ", " << pos[2] << ">\n}\n";
    }
  
  os << "}\n";

}
