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

  CLNOrientation(DYNAMO::SimData* Sim):
    CLNewton(Sim)
  {}

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

  struct rotData
  {
    CVector<> orientation;
    CVector<> angularVelocity;
  };

  const rotData& getRotData(const CParticle& part) const
  { return orientationData[part.getID()]; }
  
  void initLineOrientations(const Iflt&);

  enum
  {
    ROOT_SMALLEST_EITHER   =   1,
    ROOT_SMALLEST_POSITIVE =   2,
    ROOT_SMALLEST_NEGATIVE =   4,
    ROOT_LARGEST_EITHER    =   8,
    ROOT_LARGEST_POSITIVE  =  16,
    ROOT_LARGEST_NEGATIVE  =  32
  };

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

  virtual bool quadraticSolution(Iflt& returnVal, const int returnType, Iflt A, Iflt B, Iflt C) const;

  virtual Iflt frenkelRecursiveSearch(orientationStreamType A, orientationStreamType B, Iflt length, Iflt t_low, Iflt t_up, Iflt currentClock) const;

  virtual Iflt F_zeroDeriv(orientationStreamType A, orientationStreamType B) const;
  virtual Iflt F_firstDeriv(orientationStreamType A, orientationStreamType B) const;
  virtual Iflt F_secondDeriv(orientationStreamType A, orientationStreamType B) const;

  virtual Iflt F_firstDeriv_max(orientationStreamType A, orientationStreamType B, Iflt length) const;
  virtual Iflt F_secondDeriv_max(orientationStreamType A, orientationStreamType B, Iflt length) const;
  
  virtual void performRotation(orientationStreamType& osret, const Iflt& dt) const;

  virtual Iflt quadraticRootFinder(orientationStreamType A, orientationStreamType B, Iflt initialJump) const;
  
  virtual collisionPoints getCollisionPoints(orientationStreamType& A, orientationStreamType& B) const;
  
  mutable std::vector<rotData> orientationData;
};
#endif
