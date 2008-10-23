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

#ifndef CLNewton_H
#define CLNewton_H

#include "liouvillean.hpp"
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_real.hpp>
#include "../../extcode/include/boost/random/01_normal_distribution.hpp"

class CLNewton: public CLiouvillean
{
public:
  CLNewton(DYNAMO::SimData*);

  //Pair particle dynamics
  virtual bool SphereSphereInRoot(CPDData&, const Iflt&) const;
  virtual bool SphereSphereOutRoot(CPDData&, const Iflt&) const;  
  virtual bool sphereOverlap(const CPDData&, const Iflt&) const;

  virtual void streamParticle(CParticle&, const Iflt&) const;

  virtual intPart getSquareCellCollision(const CParticle&, 
					 const CVector<>&, 
					 const CVector<>&
					 ) const;

  virtual Iflt getSquareCellCollision2(const CParticle&, 
				       const CVector<>&, 
				       const CVector<>&
				       ) const;

  virtual size_t getSquareCellCollision3(const CParticle&, 
				       const CVector<>&, 
				       const CVector<>&
				       ) const;
  
  virtual C2ParticleData SmoothSpheresColl(const CIntEvent&, const Iflt&, const Iflt&, const EEventType& eType) const;
  
  virtual C2ParticleData SphereWellEvent(const CIntEvent&, const Iflt&, const Iflt&) const;

  virtual Iflt getHalfBoxTraversalTime(const CParticle&) const;

  virtual Iflt getWallCollision(const CParticle&, 
				const CVector<>&, 
				const CVector<>&
				  ) const;

  virtual C1ParticleData runWallCollision(const CParticle&, 
					  const CVector<>&,
					  const Iflt&
					  ) const;

  virtual C1ParticleData runAndersenWallCollision(const CParticle&, 
						  const CVector<>&,
						  const Iflt& T
						  ) const;

  virtual C1ParticleData randomGaussianEvent(const CParticle&, const Iflt&) const;

  virtual CLiouvillean* Clone() const { return new CLNewton(*this); }

protected:
  virtual void outputXML(xmlw::XmlStream& ) const;
  
  mutable boost::variate_generator<DYNAMO::baseRNG&, boost::normal_distribution_01<Iflt> > normal_sampler;
  mutable boost::variate_generator<DYNAMO::baseRNG&, boost::uniform_real<Iflt> > uniform_sampler;  
};
#endif
