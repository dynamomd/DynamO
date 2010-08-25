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

#ifndef CLSLLOD_H
#define CLSLLOD_H

#include "liouvillean.hpp"


class LSLLOD: public Liouvillean
{
public:
  LSLLOD(DYNAMO::SimData*);

  //Structure dynamics
  virtual CNParticleData multibdyCollision(const CRange&, const CRange&, 
					   const Iflt&, 
					   const EEventType&) const;
  
  virtual CNParticleData multibdyWellEvent(const CRange&, const CRange&, 
					   const Iflt&, const Iflt&, 
					   EEventType&) const;

  //Pair particle dynamics
  virtual bool SphereSphereInRoot(CPDData&, const Iflt&) const;
  virtual bool SphereSphereOutRoot(CPDData&, const Iflt&) const;  
  virtual bool sphereOverlap(const CPDData&, const Iflt&) const;

  virtual void streamParticle(Particle&, const Iflt&) const;

  virtual Iflt getSquareCellCollision2(const Particle&, 
				       const Vector &, 
				       const Vector &
				       ) const;

  virtual int getSquareCellCollision3(const Particle&, 
				      const Vector &, 
				      const Vector &
				      ) const;
  
  virtual C2ParticleData SmoothSpheresColl(const IntEvent&, const Iflt&, const Iflt&, const EEventType& eType) const;

  virtual bool DSMCSpheresTest(const Particle&, const Particle&, 
			       Iflt&, const Iflt&, CPDData&) const;

  virtual C2ParticleData DSMCSpheresRun(const Particle&, const Particle&, 
					const Iflt&, CPDData&) const;
  
  virtual C2ParticleData SphereWellEvent(const IntEvent&, const Iflt&, const Iflt&) const;

  virtual Iflt getWallCollision(const Particle&, 
				const Vector &, 
				const Vector &
				  ) const;

  virtual C1ParticleData runWallCollision(const Particle&, 
					  const Vector &,
					  const Iflt&
					  ) const;

  virtual C1ParticleData runAndersenWallCollision(const Particle&, 
						  const Vector &,
						  const Iflt& T
						  ) const;

  virtual C1ParticleData randomGaussianEvent(const Particle&, const Iflt&) const;

  virtual Liouvillean* Clone() const { return new LSLLOD(*this); }

protected:
  virtual void outputXML(xmlw::XmlStream& ) const;
};
#endif
