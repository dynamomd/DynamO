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

xmlw::XmlStream& operator<<(xmlw::XmlStream& XML, const Liouvillean& g)
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

C2ParticleData 
Liouvillean::runLineLineCollision(const CIntEvent&,
				   const Iflt&, const Iflt&) const
{ D_throw() << "Not implemented for this Liouvillean."; }

bool
Liouvillean::getLineLineCollision(CPDData&, const Iflt&, 
				   const CParticle&, const CParticle&
				   ) const
{ D_throw() << "Not implemented for this Liouvillean."; }

Iflt 
Liouvillean::getPBCSentinelTime(const CParticle&, const Iflt&) const
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
	  
	  Sim->vParticleList.push_back(CParticle(pos, vel, ID));

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
	  
	  CParticle part(xBrowseNode, i);
	  part.scaleVelocity(Sim->dynamics.units().unitVelocity());
	  part.scalePosition(Sim->dynamics.units().unitLength());
	  Sim->vParticleList.push_back(part);
	  ++prog;
	}
      if (outofsequence)
	I_cout() << IC_red << "Particle ID's out of sequence!\n"
		 << IC_red << "This can result in incorrect capture map loads etc.\n"
		 << IC_red << "Erase any capture maps in the configuration file so they are regenerated."
		 << IC_reset;            
    }
  Sim->lN = Sim->vParticleList.size();

  I_cout() << "Particle count " << Sim->lN;
}

void 
Liouvillean::outputParticleXMLData(xmlw::XmlStream& XML) const
{
  if (Sim->binaryXML)
    {
      XML << xmlw::tag("ParticleData")
	  << xmlw::attr("N") << Sim->lN
	  << xmlw::attr("AttachedBinary") << "Y"
	  << xmlw::endtag("ParticleData")
	  << xmlw::tag("AppendedBinaryVelPos")
	  << xmlw::chardata();

      {//have to scope out the iostream writes before closing the XML
	boost::iostreams::filtering_ostream base64Convertor;
	base64Convertor.push(boost::iostreams::base64_encoder());
	base64Convertor.push(boost::iostreams::line_wrapping_output_filter(80));
	base64Convertor.push(boost::iostreams::stream_sink<std::ostream>(XML.getUnderlyingStream()));
	
	boost::progress_display prog(Sim->lN);
	
	BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
	  {
	    CParticle tmp(part);
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

      XML << "\n" << xmlw::endtag("AppendedBinaryVelPos");
    }
  else
    {
      XML << xmlw::tag("ParticleData")
	  << xmlw::attr("N") << Sim->lN
	  << xmlw::attr("AttachedBinary") << ("N");

      I_cout() << "Writing Particles ";
      
      boost::progress_display prog(Sim->lN);
      
      for (unsigned long i = 0; i < Sim->lN; ++i)
	{
	  CParticle tmp(Sim->vParticleList[i]);
	  Sim->dynamics.BCs().applyBC(tmp.getPosition(), tmp.getVelocity());
	  
	  tmp.scaleVelocity(1.0 / Sim->dynamics.units().unitVelocity());
	  tmp.scalePosition(1.0 / Sim->dynamics.units().unitLength());
	  
	  XML << xmlw::tag("Pt") << tmp;

	  extraXMLParticleData(XML, i);

	  XML << xmlw::endtag("Pt");

	  ++prog;
	}

      XML << xmlw::endtag("ParticleData");
    }

  extraXMLData(XML);
}

Iflt 
Liouvillean::getParticleKineticEnergy(const CParticle& part) const
{
  return 0.5 * (part.getVelocity().nrm2()) * Sim->dynamics.getSpecies(part).getMass();
}

Vector  
Liouvillean::getVectorParticleKineticEnergy(const CParticle& part) const
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

  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    sumEnergy += getParticleKineticEnergy(part);

  return sumEnergy;
}

Vector 
Liouvillean::getVectorSystemKineticEnergy() const
{
  Vector  sumEnergy(0,0,0);

  BOOST_FOREACH(const CParticle& part, Sim->vParticleList)
    sumEnergy += getVectorParticleKineticEnergy(part);

  return sumEnergy;
}

Iflt 
Liouvillean::getkT() const
{
  return 2.0 * getSystemKineticEnergy() / (Sim->vParticleList.size() * static_cast<Iflt>(this->getParticleDOF()));
}

void
Liouvillean::rescaleSystemKineticEnergy(const Iflt& scale)
{
  Iflt scalefactor(sqrt(scale));

  BOOST_FOREACH(CParticle& part, Sim->vParticleList)
    part.getVelocity() *= scalefactor;
}

void
Liouvillean::rescaleSystemKineticEnergy(const Vector& scalefactors)
{
  BOOST_FOREACH(CParticle& part, Sim->vParticleList)
    for (size_t iDim(0); iDim < NDIM; ++iDim)
      part.getVelocity()[iDim] *= scalefactors[iDim];
}

C2ParticleData 
Liouvillean::parallelCubeColl(const CIntEvent& event, 
			       const Iflt& e, 
			       const Iflt& d, 
			       const EEventType& eType) const
{ D_throw() << "Not Implemented"; }

C2ParticleData 
Liouvillean::parallelCubeColl(const CIntEvent& event, 
			       const Iflt& e, 
			       const Iflt& d,
			       const Matrix& Rot,
			       const EEventType& eType) const
{ D_throw() << "Not Implemented"; }

Iflt 
Liouvillean::getPointPlateCollision(const CParticle& np1, const Vector& nrw0,
				     const Vector& nhat, const Iflt& Delta,
				     const Iflt& Omega, const Iflt& Sigma,
				     const Iflt& t, bool) const
{
  D_throw() << "Not Implemented";
}

C1ParticleData 
Liouvillean::runOscilatingPlate
(const CParticle& part, const Vector& rw0, const Vector& nhat, Iflt& delta, 
 const Iflt& omega0, const Iflt& sigma, const Iflt& mass, const Iflt& e,
 Iflt& t, bool strongPlate) const
{
  D_throw() << "Not Implemented";
}


Iflt 
Liouvillean::getCylinderWallCollision(const CParticle& part, 
				       const Vector& origin, 
				       const Vector& norm,
				       const Iflt&
				       ) const
{
  D_throw() << "Not Implemented";
}

C1ParticleData 
Liouvillean::runCylinderWallCollision(const CParticle&, 
				       const Vector &,
				       const Vector &,
				       const Iflt&
				       ) const
{
  D_throw() << "Not Implemented";
}

C1ParticleData 
Liouvillean::runSphereWallCollision(const CParticle&, 
				    const Vector &,
				    const Iflt&
				    ) const
{
  D_throw() << "Not Implemented";
}

C2ParticleData 
Liouvillean::SmoothSpheresCollInfMassSafe(const CIntEvent&, const Iflt&, 
					   const Iflt&, const EEventType&) const
{
  D_throw() << "Not Implemented";
}

C2ParticleData 
Liouvillean::RoughSpheresColl(const CIntEvent& event, 
			      const Iflt& e, 
			      const Iflt& et, 
			      const Iflt& d2, 
			      const EEventType& eType
			      ) const
{
  D_throw() << "Not Implemented, you need rotational dynamics";
}
