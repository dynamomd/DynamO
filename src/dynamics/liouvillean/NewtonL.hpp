/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

  virtual void streamParticle(CParticle&, const Iflt&) const;

  virtual Iflt getSquareCellCollision2(const CParticle&, 
				       const Vector &, 
				       const Vector &
				       ) const;

  virtual int getSquareCellCollision3(const CParticle&, 
				      const Vector &, 
				      const Vector &
				      ) const;
  
  virtual Iflt getPointPlateCollision(const CParticle& np1, const Vector& nrw0,
				      const Vector& nhat, const Iflt& Delta,
				      const Iflt& Omega, const Iflt& Sigma,
				      const Iflt& t, bool) const;

  virtual C1ParticleData runOscilatingPlate
  (const CParticle& part, const Vector& rw0, const Vector& nhat, Iflt& delta, 
   const Iflt& omega0, const Iflt& sigma, const Iflt& mass, const Iflt& e, 
   Iflt& t, bool strongPlate) const;

  virtual Iflt getPBCSentinelTime(const CParticle&, const Iflt&) const;

  virtual C2ParticleData SmoothSpheresColl(const CIntEvent&, const Iflt&, 
					   const Iflt&, 
					   const EEventType& eType) const;

  virtual C2ParticleData SmoothSpheresCollInfMassSafe(const CIntEvent&, const Iflt&, 
						      const Iflt&,
						      const EEventType&) const;

  virtual bool DSMCSpheresTest(const CParticle&, const CParticle&, 
			       Iflt&, const Iflt&, CPDData&) const;

  virtual C2ParticleData DSMCSpheresRun(const CParticle&, const CParticle&, 
					const Iflt&, CPDData&) const;
  
  virtual C2ParticleData SphereWellEvent(const CIntEvent&, const Iflt&, 
					 const Iflt&) const;

  virtual Iflt getWallCollision(const CParticle&, 
				const Vector &, 
				const Vector &) const;

  virtual Iflt getCylinderWallCollision(const CParticle&, 
					const Vector &, 
					const Vector &,
					const Iflt&
					) const;

  virtual C1ParticleData runCylinderWallCollision(const CParticle&, 
						  const Vector &,
						  const Vector &,
						  const Iflt&
						  ) const;

  virtual C1ParticleData runWallCollision(const CParticle&, 
					  const Vector &,
					  const Iflt&
					  ) const;

  virtual C1ParticleData runAndersenWallCollision(const CParticle&, 
						  const Vector &,
						  const Iflt& T
						  ) const;

  virtual C1ParticleData randomGaussianEvent(const CParticle&, const Iflt&) const;

  //Structure dynamics
  virtual CNParticleData multibdyCollision(const CRange&, const CRange&,
					   const Iflt&,
					   const EEventType&) const;

  virtual CNParticleData multibdyWellEvent(const CRange&, const CRange&, 
					   const Iflt&, const Iflt&, 
					   EEventType&) const;

  //Cloning
  virtual Liouvillean* Clone() const { return new LNewtonian(*this); }

  virtual C2ParticleData parallelCubeColl(const CIntEvent& event, 
					  const Iflt& e, 
					  const Iflt& d, 
					  const EEventType& eType) const;

  virtual C2ParticleData parallelCubeColl(const CIntEvent& event,
					  const Iflt& e, const Iflt& d,
					  const Matrix& rot,
					  const EEventType& eType = CORE) const;

protected:
  virtual void outputXML(xmlw::XmlStream& ) const;
};
#endif
