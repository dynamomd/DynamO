/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
  virtual NEventData multibdyCollision(const CRange&, const CRange&, 
					   const double&, 
					   const EEventType&) const;
  
  virtual NEventData multibdyWellEvent(const CRange&, const CRange&, 
					   const double&, const double&, 
					   EEventType&) const;

  //Pair particle dynamics
  virtual bool SphereSphereInRoot(CPDData&, const double&, bool p1Dynamic, bool p2Dynamic) const;
  virtual bool SphereSphereOutRoot(CPDData&, const double&, bool p1Dynamic, bool p2Dynamic) const;  
  virtual bool sphereOverlap(const CPDData&, const double&) const;

  virtual void streamParticle(Particle&, const double&) const;

  virtual double getSquareCellCollision2(const Particle&, 
				       const Vector &, 
				       const Vector &
				       ) const;

  virtual int getSquareCellCollision3(const Particle&, 
				      const Vector &, 
				      const Vector &
				      ) const;
  
  virtual PairEventData SmoothSpheresColl(const IntEvent&, const double&, const double&, const EEventType& eType) const;

  virtual bool DSMCSpheresTest(const Particle&, const Particle&, 
			       double&, const double&, CPDData&) const;

  virtual PairEventData DSMCSpheresRun(const Particle&, const Particle&, 
					const double&, CPDData&) const;
  
  virtual PairEventData SphereWellEvent(const IntEvent&, const double&, const double&) const;

  virtual double getWallCollision(const Particle&, 
				const Vector &, 
				const Vector &
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

  virtual Liouvillean* Clone() const { return new LSLLOD(*this); }

protected:
  virtual void outputXML(xml::XmlStream& ) const;
};
#endif
