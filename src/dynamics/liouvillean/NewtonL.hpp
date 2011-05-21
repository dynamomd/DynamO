/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>
    Copyright (C) 2011  Sebastian Gonzalez

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
#include "liouvillean.hpp"

class LNewtonian: public Liouvillean
{
public:
  LNewtonian(dynamo::SimData*);

  //Pair particle dynamics
  virtual bool SphereSphereInRoot(CPDData&, const double&, bool p1Dynamic, bool p2Dynamic) const;
  virtual bool SphereSphereOutRoot(CPDData&, const double&, bool p1Dynamic, bool p2Dynamic) const;  
  virtual bool sphereOverlap(const CPDData&, const double&) const;

  virtual bool CubeCubeInRoot(CPDData&, const double&) const;

  virtual bool cubeOverlap(const CPDData&, const double&) const;

  virtual void streamParticle(Particle&, const double&) const;

  virtual double getSquareCellCollision2(const Particle&, 
				       const Vector &, 
				       const Vector &
				       ) const;

  virtual int getSquareCellCollision3(const Particle&, 
				      const Vector &, 
				      const Vector &
				      ) const;
  
  virtual std::pair<bool,double>
  getPointPlateCollision(const Particle& np1, const Vector& nrw0,
			 const Vector& nhat, const double& Delta,
			 const double& Omega, const double& Sigma,
			 const double& t, bool) const;

  virtual ParticleEventData runOscilatingPlate
  (const Particle& part, const Vector& rw0, const Vector& nhat, double& delta, 
   const double& omega0, const double& sigma, const double& mass, const double& e, 
   double& t, bool strongPlate) const;

  virtual double getPBCSentinelTime(const Particle&, const double&) const;

  virtual PairEventData SmoothSpheresColl(const IntEvent&, const double&, 
					  const double&, 
					  const EEventType& eType) const;

  virtual bool DSMCSpheresTest(const Particle&, const Particle&, 
			       double&, const double&, CPDData&) const;

  virtual PairEventData DSMCSpheresRun(const Particle&, const Particle&, 
					const double&, CPDData&) const;
  
  virtual PairEventData SphereWellEvent(const IntEvent&, const double&, 
					 const double&) const;

  virtual double getWallCollision(const Particle&, 
				const Vector &, 
				const Vector &) const;

  virtual double getParticleTriangleEvent(const Particle& part, 
					  const Vector & A, 
					  const Vector & B, 
					  const Vector & C,
					  const double dist
					  ) const;

  virtual double getCylinderWallCollision(const Particle&, 
					const Vector &, 
					const Vector &,
					const double&
					) const;

  virtual ParticleEventData runCylinderWallCollision(const Particle&, 
						  const Vector &,
						  const Vector &,
						  const double&
						  ) const;

  virtual ParticleEventData runSphereWallCollision(const Particle&, 
						const Vector &,
						const double&
						) const;

  virtual ParticleEventData runWallCollision(const Particle&, 
					  const Vector &,
					  const double&
					  ) const;

  virtual ParticleEventData runAndersenWallCollision(const Particle&, 
						  const Vector &,
						  const double& T
						  ) const;

  virtual ParticleEventData randomGaussianEvent(const Particle&, const double&) const;

  //Structure dynamics
  virtual NEventData multibdyCollision(const CRange&, const CRange&,
					   const double&,
					   const EEventType&) const;

  virtual NEventData multibdyWellEvent(const CRange&, const CRange&, 
					   const double&, const double&, 
					   EEventType&) const;

  //Cloning
  virtual Liouvillean* Clone() const { return new LNewtonian(*this); }

  virtual PairEventData parallelCubeColl(const IntEvent& event,
					  const double& e, const double& d,
					  const Matrix& rot,
					  const EEventType& eType = CORE) const;

  virtual bool getLineLineCollision(CPDData& PD, const double& length, 
				    const Particle& p1, const Particle& p2
				    ) const;
  
  virtual PairEventData runLineLineCollision(const IntEvent& eevent, 
					     const double& elasticity, const double& length) const;

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

  virtual bool getOffCenterSphereOffCenterSphereCollision(CPDData& PD, const double& length, 
							  const double& diameter,
							  const Particle& p1, const Particle& p2
							  ) const;
  
  virtual PairEventData runOffCenterSphereOffCenterSphereCollision(const IntEvent& eevent, 
								   const double& elasticity, 
								   const double& length, 
								   const double& diameter) const;

protected:
  virtual void outputXML(xml::XmlStream& ) const;

  mutable long double lastAbsoluteClock;
  mutable unsigned int lastCollParticle1;
  mutable unsigned int lastCollParticle2;  

};
