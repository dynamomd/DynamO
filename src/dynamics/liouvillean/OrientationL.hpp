/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>
    Copyright (C) 2011  Sebastian Gonzalez <tsuresuregusa@gmail.com>

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

#pragma once
#include "NewtonL.hpp"
#include "../../datatypes/vector.hpp"
#include <vector>

class CLinesFunc;
class CDumbbellsFunc;
class CShape;

class LNOrientation: public LNewtonian
{
public:  
  LNOrientation(DYNAMO::SimData* Sim, const magnet::xml::Node& XML):
    LNewtonian(Sim) {}

  LNOrientation(DYNAMO::SimData* Sim): LNewtonian(Sim) {}

  virtual void initialise();

  virtual Liouvillean* Clone() const { return new LNOrientation(*this); }

  virtual void loadParticleXMLData(const magnet::xml::Node&);

  virtual bool getLineLineCollision(CPDData& PD, const double& length, 
				    const Particle& p1, const Particle& p2
				    ) const;
  
  virtual PairEventData runLineLineCollision(const IntEvent& eevent, 
					     const double& elasticity, const double& length) const;
  
  virtual bool getOffCenterSphereOffCenterSphereCollision(CPDData& PD, const double& length, const double& diameter,

							  const Particle& p1, const Particle& p2
							  ) const;
  
  virtual PairEventData runOffCenterSphereOffCenterSphereCollision(const IntEvent& eevent, 
								   const double& elasticity, const double& length, const double& diameter) const;

  
  virtual ParticleEventData runAndersenWallCollision(const Particle& part, 
						  const Vector & vNorm,
						  const double& sqrtT
						  ) const;
  
  virtual ParticleEventData randomGaussianEvent(const Particle& part, 
					     const double& sqrtT) const;

  struct rotData
  {
    Vector  orientation;
    Vector  angularVelocity;
  };

  const rotData& getRotData(const Particle& part) const
  { return orientationData[part.getID()]; }

  const std::vector<rotData>& getCompleteRotData() const
  { return orientationData; }
  
  void initLineOrientations(const double&);

  virtual PairEventData RoughSpheresColl(const IntEvent& event, 
					  const double& e, 
					  const double& et, 
					  const double& d2, 
					  const EEventType& eType = CORE
					  ) const;

  virtual ParticleEventData runRoughWallCollision(const Particle& part, 
					       const Vector & vNorm,
					       const double& e,
					       const double& et,
					       const double& r
					       ) const;

protected:

  virtual void extraXMLParticleData(xml::XmlStream&, const size_t) const;

  virtual void extraXMLData(xml::XmlStream&) const;

  virtual void outputXML(xml::XmlStream&) const;

  virtual void streamParticle(Particle&, const double&) const;
  
  virtual size_t getParticleDOF() const;
  virtual double getParticleKineticEnergy(const Particle& part) const;
  virtual void rescaleSystemKineticEnergy(const double&);
  
  mutable std::vector<rotData> orientationData;
};
