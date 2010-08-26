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

#ifndef LNewtonian_H
#define LNewtonian_H

#include "liouvillean.hpp"


class LNewtonian: public Liouvillean
{
public:
  LNewtonian(DYNAMO::SimData*);

  //Pair particle dynamics
  virtual bool SphereSphereInRoot(CPDData&, const Iflt&) const;
  virtual bool SphereSphereOutRoot(CPDData&, const Iflt&) const;  
  virtual bool sphereOverlap(const CPDData&, const Iflt&) const;

  virtual bool CubeCubeInRoot(CPDData&, const Iflt&) const;
  virtual bool CubeCubeInRoot(CPDData&, const Iflt&, const Matrix&) const;

  virtual bool cubeOverlap(const CPDData&, const Iflt&) const;
  virtual bool cubeOverlap(const CPDData&, const Iflt&, const Matrix&) const;

  virtual void streamParticle(Particle&, const Iflt&) const;

  virtual Iflt getSquareCellCollision2(const Particle&, 
				       const Vector &, 
				       const Vector &
				       ) const;

  virtual int getSquareCellCollision3(const Particle&, 
				      const Vector &, 
				      const Vector &
				      ) const;
  
  virtual Iflt getPointPlateCollision(const Particle& np1, const Vector& nrw0,
				      const Vector& nhat, const Iflt& Delta,
				      const Iflt& Omega, const Iflt& Sigma,
				      const Iflt& t, bool) const;

  virtual ParticleEventData runOscilatingPlate
  (const Particle& part, const Vector& rw0, const Vector& nhat, Iflt& delta, 
   const Iflt& omega0, const Iflt& sigma, const Iflt& mass, const Iflt& e, 
   Iflt& t, bool strongPlate) const;

  virtual Iflt getPBCSentinelTime(const Particle&, const Iflt&) const;

  virtual PairEventData SmoothSpheresColl(const IntEvent&, const Iflt&, 
					   const Iflt&, 
					   const EEventType& eType) const;

  virtual PairEventData SmoothSpheresCollInfMassSafe(const IntEvent&, const Iflt&, 
						      const Iflt&,
						      const EEventType&) const;

  virtual bool DSMCSpheresTest(const Particle&, const Particle&, 
			       Iflt&, const Iflt&, CPDData&) const;

  virtual PairEventData DSMCSpheresRun(const Particle&, const Particle&, 
					const Iflt&, CPDData&) const;
  
  virtual PairEventData SphereWellEvent(const IntEvent&, const Iflt&, 
					 const Iflt&) const;

  virtual Iflt getWallCollision(const Particle&, 
				const Vector &, 
				const Vector &) const;

  virtual Iflt getCylinderWallCollision(const Particle&, 
					const Vector &, 
					const Vector &,
					const Iflt&
					) const;

  virtual ParticleEventData runCylinderWallCollision(const Particle&, 
						  const Vector &,
						  const Vector &,
						  const Iflt&
						  ) const;

  virtual ParticleEventData runSphereWallCollision(const Particle&, 
						const Vector &,
						const Iflt&
						) const;

  virtual ParticleEventData runWallCollision(const Particle&, 
					  const Vector &,
					  const Iflt&
					  ) const;

  virtual ParticleEventData runAndersenWallCollision(const Particle&, 
						  const Vector &,
						  const Iflt& T
						  ) const;

  virtual ParticleEventData randomGaussianEvent(const Particle&, const Iflt&) const;

  //Structure dynamics
  virtual NEventData multibdyCollision(const CRange&, const CRange&,
					   const Iflt&,
					   const EEventType&) const;

  virtual NEventData multibdyWellEvent(const CRange&, const CRange&, 
					   const Iflt&, const Iflt&, 
					   EEventType&) const;

  //Cloning
  virtual Liouvillean* Clone() const { return new LNewtonian(*this); }

  virtual PairEventData parallelCubeColl(const IntEvent& event, 
					  const Iflt& e, 
					  const Iflt& d, 
					  const EEventType& eType) const;

  virtual PairEventData parallelCubeColl(const IntEvent& event,
					  const Iflt& e, const Iflt& d,
					  const Matrix& rot,
					  const EEventType& eType = CORE) const;

protected:
  virtual void outputXML(xml::XmlStream& ) const;
};
#endif
