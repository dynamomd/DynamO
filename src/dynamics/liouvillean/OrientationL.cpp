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

#include "OrientationL.hpp"
#include "../2particleEventData.hpp"
#include <boost/progress.hpp>
#include "../../datatypes/vector.xml.hpp"

void
CLNOrientation::operator<<(const XMLNode& XML)
{
  XMLNode xSubNode(XML.getParentNode().getParentNode()
		   .getChildNode("ParticleData"));
  
  size_t nPart = xSubNode.nChildNode("Pt");
  
  orientationData.resize(nPart);
  
  I_cout() << "Loading orientation data....";

  boost::progress_display prog(nPart);
  
  int xml_iter = 0;

  for (size_t i = 0; i < nPart; ++i)
    {
      XMLNode xBrowseNode(xSubNode.getChildNode("Pt", &xml_iter));
      
      orientationData[i].orientation << xBrowseNode.getChildNode("U");
      orientationData[i].angularVelocity << xBrowseNode.getChildNode("O");
      ++prog;
    }
}

void
CLNOrientation::outputXML(xmlw::XmlStream& XML) const
{
  XML << xmlw::attr("Type") 
      << "NOrientation";
}

void 
CLNOrientation::outputExtraPDatXML(xmlw::XmlStream& XML,
				   const CParticle& part) const
{
  XML << xmlw::tag("O")
      << orientationData[part.getID()].angularVelocity
      << xmlw::endtag("O")
      << xmlw::tag("U")
      << orientationData[part.getID()].orientation
      << xmlw::endtag("U");
}

bool 
CLNOrientation::getLineLineCollision(const CPDData&, const Iflt&, 
				     const CParticle&, const CParticle&,
				     const Iflt&
				     ) const
{ D_throw() << "Not implemented"; }

C2ParticleData 
CLNOrientation::runLineLineCollision(const CIntEvent&) const
{ D_throw() << "Not implemented"; }

void 
CLNOrientation::streamParticle(CParticle& part, const Iflt& dt) const
{
  if(NDIM != 3) { D_throw() << "Implemented only for 3D rotations"; }
  else
  {
  
    //First hand over to the newtonian dynamics
    CLNewton::streamParticle(part, dt);
  
    //Now stream the orientation dynamics
    
    Iflt angle = orientationData[part.getID()].angularVelocity.length() * dt;
	
	  Iflt v1 = (orientationData[part.getID()].angularVelocity.unitVector())[0];
  	Iflt v2 = (orientationData[part.getID()].angularVelocity.unitVector())[1];
	  Iflt v3 = (orientationData[part.getID()].angularVelocity.unitVector())[2];
	
	  // axis is not undefined and angle is not zero
  	if(!(v1 == 0 && v2 == 0 && v3 == 0) && angle != 0)
	  {	
		
		  Iflt matrix[3][3];
		
  	  Iflt cos_term = cos(angle);
	  	Iflt sin_term = sin(angle);
	
	    matrix[0][0] = pow(v1, 2) + (pow(v2 ,2) + pow(v3, 2))*(cos_term);
	    matrix[0][1] = (v1 * v2 * (1 - cos_term)) - (v3 * sin_term);
	    matrix[0][2] = (v1 * v3 * (1 - cos_term)) + (v2 * sin_term);
	  	matrix[1][0] = (v1 * v2 * (1 - cos_term)) + (v3 * sin_term);
	  	matrix[1][1] = pow(v2, 2) + (pow(v3 ,2) + pow(v1, 2))*(cos_term);
	  	matrix[1][2] = (v2 * v3 * (1 - cos_term)) - (v1 * sin_term);
	  	matrix[2][0] = (v3 * v1 * (1 - cos_term)) - (v2 * sin_term);
	  	matrix[2][1] = (v2 * v3 * (1 - cos_term)) + (v1 * sin_term);
	  	matrix[2][2] = pow(v3, 2) + (pow(v1 ,2) + pow(v2, 2))*(cos_term);
		
	  	CVector<> tempvec;
	  	tempvec[0] = (matrix[0][0] * orientationData[part.getID()].orientation[0]) + (matrix[0][1] * orientationData[part.getID()].orientation[1]) + (matrix[0][2] * orientationData[part.getID()].orientation[2]);
	  	tempvec[1] = (matrix[1][0] * orientationData[part.getID()].orientation[0]) + (matrix[1][1] * orientationData[part.getID()].orientation[1]) + (matrix[1][2] * orientationData[part.getID()].orientation[2]);
	  	tempvec[2] = (matrix[2][0] * orientationData[part.getID()].orientation[0]) + (matrix[2][1] * orientationData[part.getID()].orientation[1]) + (matrix[2][2] * orientationData[part.getID()].orientation[2]);
		
  		orientationData[part.getID()].orientation = tempvec;
		}
	}
}

C1ParticleData 
CLNOrientation::runAndersenWallCollision(const CParticle& part, 
			 const CVector<>& vNorm,
			 const Iflt& sqrtT
			 ) const
{
  D_throw() << "Need to implement thermostating of the rotational degrees"
    " of freedom";
}
  
C1ParticleData 
CLNOrientation::randomGaussianEvent(const CParticle& part, 
				    const Iflt& sqrtT) const
{
  D_throw() << "Need to implement thermostating of the rotational degrees"
    " of freedom";  
}
