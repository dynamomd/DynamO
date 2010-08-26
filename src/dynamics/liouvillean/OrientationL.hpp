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

#ifndef LNOrientation_H
#define LNOrientation_H

#include "NewtonL.hpp"
#include <vector>
#include "../../datatypes/vector.hpp"

class CLinesFunc;
class CShape;

class LNOrientation: public LNewtonian
{
public:  
  LNOrientation(DYNAMO::SimData* Sim, const XMLNode& XML):
    LNewtonian(Sim),
    lastAbsoluteClock(-1),
    lastCollParticle1(0),
    lastCollParticle2(0)
  {}

  LNOrientation(DYNAMO::SimData* Sim):
    LNewtonian(Sim),
    lastAbsoluteClock(-1),
    lastCollParticle1(0),
    lastCollParticle2(0)
  {}

  virtual void initialise();

  virtual Liouvillean* Clone() const { return new LNOrientation(*this); }

  virtual void loadParticleXMLData(const XMLNode&);

  virtual bool getLineLineCollision(CPDData& PD, const Iflt& length, 
				    const Particle& p1, const Particle& p2
				    ) const;
  
  virtual PairEventData runLineLineCollision(const IntEvent& eevent, 
					      const Iflt& elasticity, const Iflt& length) const;
  
  virtual ParticleEventData runAndersenWallCollision(const Particle& part, 
						  const Vector & vNorm,
						  const Iflt& sqrtT
						  ) const;
  
  virtual ParticleEventData randomGaussianEvent(const Particle& part, 
					     const Iflt& sqrtT) const;

  struct rotData
  {
    Vector  orientation;
    Vector  angularVelocity;
  };

  const rotData& getRotData(const Particle& part) const
  { return orientationData[part.getID()]; }

  const std::vector<rotData>& getCompleteRotData() const
  { return orientationData; }
  
  void initLineOrientations(const Iflt&);

  virtual PairEventData RoughSpheresColl(const IntEvent& event, 
					  const Iflt& e, 
					  const Iflt& et, 
					  const Iflt& d2, 
					  const EEventType& eType = CORE
					  ) const;

  virtual ParticleEventData runRoughWallCollision(const Particle& part, 
					       const Vector & vNorm,
					       const Iflt& e,
					       const Iflt& et,
					       const Iflt& r
					       ) const;

protected:

  virtual void extraXMLParticleData(xmlw::XmlStream&, const size_t) const;

  virtual void extraXMLData(xmlw::XmlStream&) const;

  virtual void outputXML(xmlw::XmlStream&) const;

  virtual void streamParticle(Particle&, const Iflt&) const;
  
  virtual size_t getParticleDOF() const;
  virtual Iflt getParticleKineticEnergy(const Particle& part) const;
  virtual void rescaleSystemKineticEnergy(const Iflt&);
  
  mutable std::vector<rotData> orientationData;
  mutable lIflt lastAbsoluteClock;
  mutable unsigned int lastCollParticle1;
  mutable unsigned int lastCollParticle2;  
};
#endif
