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

#include "include.hpp"
#include "../../extcode/xmlwriter.hpp"
#include "../../extcode/xmlParser.h"
#include "../../base/is_exception.hpp"
#include "../../base/is_simdata.hpp"
#include "../2particleEventData.hpp"
#include "../units/units.hpp"
#include <boost/foreach.hpp>
#include <boost/progress.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/iostreams/filter/base64.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/chain.hpp>
#include <boost/iostreams/device/stream_sink.hpp>
#include <boost/iostreams/device/stream_source.hpp>
#include <boost/iostreams/filter/base64cleaner.hpp>
#include <boost/iostreams/filter/linewrapout.hpp>
#include "../../extcode/binaryHelper.hpp"

xml::XmlStream& operator<<(xml::XmlStream& XML, const Liouvillean& g)
{
  g.outputXML(XML);
  return XML;
}

Liouvillean* 
Liouvillean::loadClass(const XMLNode& XML, DYNAMO::SimData* tmp)
{
  if (!strcmp(XML.getAttribute("Type"),"Newtonian"))
    return new LNewtonian(tmp);
  if (!strcmp(XML.getAttribute("Type"),"NewtonianGravity"))
    return new LNewtonianGravity(tmp, XML);
  else if (!strcmp(XML.getAttribute("Type"),"NOrientation"))
    return new LNOrientation(tmp, XML);
  else if (!strcmp(XML.getAttribute("Type"),"SLLOD"))
    return new LSLLOD(tmp);
  else if (!strcmp(XML.getAttribute("Type"),"NewtonianMC"))
    return new LNewtonianMC(tmp, XML);
  else 
    D_throw() << XML.getAttribute("Type")
	      << ", Unknown type of Liouvillean encountered";
}

PairEventData 
Liouvillean::runLineLineCollision(const IntEvent&,
				   const Iflt&, const Iflt&) const
{ D_throw() << "Not implemented for this Liouvillean."; }

bool
Liouvillean::getLineLineCollision(CPDData&, const Iflt&, 
				   const Particle&, const Particle&
				   ) const
{ D_throw() << "Not implemented for this Liouvillean."; }

Iflt 
Liouvillean::getPBCSentinelTime(const Particle&, const Iflt&) const
{ D_throw() << "Not implemented for this Liouvillean."; }

void 
Liouvillean::loadParticleXMLData(const XMLNode& XML)
{
  I_cout() << "Loading Particle Data";
  std::cout.flush();
  
  XMLNode xSubNode = XML.getChildNode("ParticleData");

  if (xSubNode.isAttributeSet("AttachedBinary")
      && (std::toupper(xSubNode.getAttribute("AttachedBinary")[0]) == 'Y'))
    {
#ifndef DYNAMO_CONDOR
      Sim->binaryXML = true;

      unsigned long nPart 
	= boost::lexical_cast<unsigned long>(xSubNode.getAttribute("N"));

      boost::progress_display prog(nPart);

      boost::iostreams::filtering_istream base64Convertor;
      
      base64Convertor.push(boost::iostreams::base64_decoder());
      base64Convertor.push(boost::iostreams::base64cleaner_input_filter());
      
      {
	const char* start = XML.getChildNode("AppendedBinaryVelPos").getText();
	base64Convertor.push(boost::make_iterator_range(std::make_pair(start, start + strlen(start))));
      }

      for (unsigned long i = 0; i < nPart; ++i)
	{
	  unsigned long ID;
	  Vector  vel;
	  Vector  pos;
	  
	  binaryread(base64Convertor, ID);

	  if (i != ID) 
	    D_throw() << "Binary data corruption detected, id's don't match"
		      << "\nMight be because this file was generated on another architecture (32/64bit)"
		      << "\nTry using the --text option to avoid using the binary format";
	  
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    binaryread(base64Convertor, vel[iDim]);
	  
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    binaryread(base64Convertor, pos[iDim]);
	  
	  vel *= Sim->dynamics.units().unitVelocity();
	  pos *= Sim->dynamics.units().unitLength();
	  
	  Sim->particleList.push_back(Particle(pos, vel, ID));

	  ++prog;
	}
#else
      D_throw() << "Cannot use appended binary data in Condor executables!";
#endif
    }
  else
    {
      int xml_iter = 0;
      
      unsigned long nPart = xSubNode.nChildNode("Pt");
      boost::progress_display prog(nPart);
      bool outofsequence = false;  
      
      for (unsigned long i = 0; i < nPart; i++)
	{
	  XMLNode xBrowseNode = xSubNode.getChildNode("Pt", &xml_iter);
	  
	  if (boost::lexical_cast<unsigned long>
	      (xBrowseNode.getAttribute("ID")) != i)
	    outofsequence = true;
	  
	  Particle part(xBrowseNode, i);
	  part.scaleVelocity(Sim->dynamics.units().unitVelocity());
	  part.scalePosition(Sim->dynamics.units().unitLength());
	  Sim->particleList.push_back(part);
	  ++prog;
	}
      if (outofsequence)
	I_cout() << IC_red << "Particle ID's out of sequence!\n"
		 << IC_red << "This can result in incorrect capture map loads etc.\n"
		 << IC_red << "Erase any capture maps in the configuration file so they are regenerated."
		 << IC_reset;            
    }
  Sim->N = Sim->particleList.size();

  I_cout() << "Particle count " << Sim->N;
}

void 
Liouvillean::outputParticleXMLData(xml::XmlStream& XML) const
{
  if (Sim->binaryXML)
    {
      XML << xml::tag("ParticleData")
	  << xml::attr("N") << Sim->N
	  << xml::attr("AttachedBinary") << "Y"
	  << xml::endtag("ParticleData")
	  << xml::tag("AppendedBinaryVelPos")
	  << xml::chardata();

      {//have to scope out the iostream writes before closing the XML
	boost::iostreams::filtering_ostream base64Convertor;
	base64Convertor.push(boost::iostreams::base64_encoder());
	base64Convertor.push(boost::iostreams::line_wrapping_output_filter(80));
	base64Convertor.push(boost::iostreams::stream_sink<std::ostream>(XML.getUnderlyingStream()));
	
	boost::progress_display prog(Sim->N);
	
	BOOST_FOREACH(const Particle& part, Sim->particleList)
	  {
	    Particle tmp(part);
	    Sim->dynamics.BCs().applyBC(tmp.getPosition(), tmp.getVelocity());
	    
	    tmp.scaleVelocity(1.0 / Sim->dynamics.units().unitVelocity());
	    tmp.scalePosition(1.0 / Sim->dynamics.units().unitLength());	  
	    
	    binarywrite(base64Convertor, tmp.getID());
	    
	    for (size_t iDim(0); iDim < NDIM; ++iDim)
	      binarywrite(base64Convertor, tmp.getVelocity()[iDim]);
	    
	    for (size_t iDim(0); iDim < NDIM; ++iDim)
	      binarywrite(base64Convertor, tmp.getPosition()[iDim]);
	    
	    ++prog;
	  }
      }

      XML << "\n" << xml::endtag("AppendedBinaryVelPos");
    }
  else
    {
      XML << xml::tag("ParticleData")
	  << xml::attr("N") << Sim->N
	  << xml::attr("AttachedBinary") << ("N");

      I_cout() << "Writing Particles ";
      
      boost::progress_display prog(Sim->N);
      
      for (unsigned long i = 0; i < Sim->N; ++i)
	{
	  Particle tmp(Sim->particleList[i]);
	  Sim->dynamics.BCs().applyBC(tmp.getPosition(), tmp.getVelocity());
	  
	  tmp.scaleVelocity(1.0 / Sim->dynamics.units().unitVelocity());
	  tmp.scalePosition(1.0 / Sim->dynamics.units().unitLength());
	  
	  XML << xml::tag("Pt") << tmp;

	  extraXMLParticleData(XML, i);

	  XML << xml::endtag("Pt");

	  ++prog;
	}

      XML << xml::endtag("ParticleData");
    }

  extraXMLData(XML);
}

Iflt 
Liouvillean::getParticleKineticEnergy(const Particle& part) const
{
  return 0.5 * (part.getVelocity().nrm2()) * Sim->dynamics.getSpecies(part).getMass();
}

Vector  
Liouvillean::getVectorParticleKineticEnergy(const Particle& part) const
{
  Vector  tmp(0.5 * part.getVelocity() * Sim->dynamics.getSpecies(part).getMass());

  for (size_t i = 0; i < NDIM; ++i)
    tmp[i] *= part.getVelocity()[i];

  return tmp;
}

Iflt 
Liouvillean::getSystemKineticEnergy() const
{
  Iflt sumEnergy(0);

  BOOST_FOREACH(const Particle& part, Sim->particleList)
    sumEnergy += getParticleKineticEnergy(part);

  return sumEnergy;
}

Vector 
Liouvillean::getVectorSystemKineticEnergy() const
{
  Vector  sumEnergy(0,0,0);

  BOOST_FOREACH(const Particle& part, Sim->particleList)
    sumEnergy += getVectorParticleKineticEnergy(part);

  return sumEnergy;
}

Iflt 
Liouvillean::getkT() const
{
  return 2.0 * getSystemKineticEnergy() / (Sim->particleList.size() * static_cast<Iflt>(this->getParticleDOF()));
}

void
Liouvillean::rescaleSystemKineticEnergy(const Iflt& scale)
{
  Iflt scalefactor(sqrt(scale));

  BOOST_FOREACH(Particle& part, Sim->particleList)
    part.getVelocity() *= scalefactor;
}

void
Liouvillean::rescaleSystemKineticEnergy(const Vector& scalefactors)
{
  BOOST_FOREACH(Particle& part, Sim->particleList)
    for (size_t iDim(0); iDim < NDIM; ++iDim)
      part.getVelocity()[iDim] *= scalefactors[iDim];
}

PairEventData 
Liouvillean::parallelCubeColl(const IntEvent& event, 
			       const Iflt& e, 
			       const Iflt& d, 
			       const EEventType& eType) const
{ D_throw() << "Not Implemented"; }

PairEventData 
Liouvillean::parallelCubeColl(const IntEvent& event, 
			       const Iflt& e, 
			       const Iflt& d,
			       const Matrix& Rot,
			       const EEventType& eType) const
{ D_throw() << "Not Implemented"; }

Iflt 
Liouvillean::getPointPlateCollision(const Particle& np1, const Vector& nrw0,
				     const Vector& nhat, const Iflt& Delta,
				     const Iflt& Omega, const Iflt& Sigma,
				     const Iflt& t, bool) const
{
  D_throw() << "Not Implemented";
}

ParticleEventData 
Liouvillean::runOscilatingPlate
(const Particle& part, const Vector& rw0, const Vector& nhat, Iflt& delta, 
 const Iflt& omega0, const Iflt& sigma, const Iflt& mass, const Iflt& e,
 Iflt& t, bool strongPlate) const
{
  D_throw() << "Not Implemented";
}


Iflt 
Liouvillean::getCylinderWallCollision(const Particle& part, 
				       const Vector& origin, 
				       const Vector& norm,
				       const Iflt&
				       ) const
{
  D_throw() << "Not Implemented";
}

ParticleEventData 
Liouvillean::runCylinderWallCollision(const Particle&, 
				       const Vector &,
				       const Vector &,
				       const Iflt&
				       ) const
{
  D_throw() << "Not Implemented";
}

ParticleEventData 
Liouvillean::runSphereWallCollision(const Particle&, 
				    const Vector &,
				    const Iflt&
				    ) const
{
  D_throw() << "Not Implemented";
}

PairEventData 
Liouvillean::SmoothSpheresCollInfMassSafe(const IntEvent&, const Iflt&, 
					   const Iflt&, const EEventType&) const
{
  D_throw() << "Not Implemented";
}

PairEventData 
Liouvillean::RoughSpheresColl(const IntEvent& event, 
			      const Iflt& e, 
			      const Iflt& et, 
			      const Iflt& d2, 
			      const EEventType& eType
			      ) const
{
  D_throw() << "Not Implemented, you need rotational dynamics";
}

ParticleEventData 
Liouvillean::runRoughWallCollision(const Particle& part, 
				   const Vector & vNorm,
				   const Iflt& e,
				   const Iflt& et,
				   const Iflt& r
				   ) const
{
  D_throw() << "Not Implemented, you need rotational dynamics";
}
