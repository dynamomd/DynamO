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

  virtual bool getLineLineCollision(const CPDData& PD, const Iflt& length, 
				    const CParticle& p1, const CParticle& p2,
				    const Iflt& twindow
				    ) const;
  
  virtual C2ParticleData runLineLineCollision(const CIntEvent& eevent) const;
  
  virtual C1ParticleData runAndersenWallCollision(const CParticle& part, 
						  const CVector<>& vNorm,
						  const Iflt& sqrtT
						  ) const;
  
  virtual C1ParticleData randomGaussianEvent(const CParticle& part, 
					     const Iflt& sqrtT) const;

  virtual void operator<<(const XMLNode&);

  virtual void outputExtraPDatXML(xmlw::XmlStream&,
				  const CParticle&) const;

protected:
  virtual void outputXML(xmlw::XmlStream&) const;

  virtual void streamParticle(CParticle&, const Iflt&) const;

  struct rotData
  {
    CVector<> orientation;
    CVector<> angularVelocity;
  };
  
  struct orientationStreamType
  {
    rotData rot;
    CVector<> velocity;
    CVector<> position;
  };
  
  virtual void performRotation(orientationStreamType&, const Iflt&) const;
  
  virtual bool recursiveRootFinder(const Iflt& interpolationSize, const Iflt& window_open, const Iflt& window_closed) const;
  
  mutable std::vector<rotData> orientationData;
};
#endif
