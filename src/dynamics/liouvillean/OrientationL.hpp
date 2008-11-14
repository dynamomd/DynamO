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

#ifndef CLNOrientation_H
#define CLNOrientation_H

#include "NewtonL.hpp"
#include <vector>
#include "../../datatypes/vector.hpp"

class CLNOrientation: public CLNewton
{
public:  
  CLNOrientation(DYNAMO::SimData* Sim, const XMLNode& XML):
    CLNewton(Sim)
  {
    operator<<(XML);
  }

  virtual CLiouvillean* Clone() const { return new CLNOrientation(*this); }

  virtual bool getLineLineCollision(CPDData& PD, const Iflt& length, 
				    const CParticle& p1, const CParticle& p2
				    ) const;
  
  virtual C2ParticleData runLineLineCollision(const CIntEvent& eevent, 
					      const Iflt& length) const;
  
  virtual C1ParticleData runAndersenWallCollision(const CParticle& part, 
						  const CVector<>& vNorm,
						  const Iflt& sqrtT
						  ) const;
  
  virtual C1ParticleData randomGaussianEvent(const CParticle& part, 
					     const Iflt& sqrtT) const;

  virtual void operator<<(const XMLNode&);

  virtual void outputExtraPDatXML(xmlw::XmlStream&,
				  const CParticle&) const;

  struct rotData
  {
    CVector<> orientation;
    CVector<> angularVelocity;
  };

  const rotData& getRotData(const CParticle& part) const
  { return orientationData[part.getID()]; }

protected:
  virtual void outputXML(xmlw::XmlStream&) const;

  virtual void streamParticle(CParticle&, const Iflt&) const;

  struct orientationStreamType
  {
    rotData rot;
    CVector<> velocity;
    CVector<> position;
  };
  
  struct collisionPoints
  {
    Iflt alpha;
    Iflt beta;
  };
  
  virtual void performRotation(orientationStreamType&, const Iflt&) const;
  
  virtual bool recursiveRootFinder(orientationStreamType& A, orientationStreamType& B, const Iflt& length, 
                                   const Iflt& interpolationSize, const Iflt& windowOpen, const Iflt& windowClosed, Iflt& collisionTime) const;
  
  virtual collisionPoints getCollisionPoints(orientationStreamType& A, orientationStreamType& B) const;
  
  mutable std::vector<rotData> orientationData;
};
#endif
