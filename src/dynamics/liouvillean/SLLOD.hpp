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

#ifndef CLSLLOD_H
#define CLSLLOD_H

#include "liouvillean.hpp"


class CLSLLOD: public CLiouvillean
{
public:
  CLSLLOD(DYNAMO::SimData*);

  //Structure Dynamics
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

  virtual void streamParticle(CParticle&, const Iflt&) const;

  virtual Iflt getSquareCellCollision2(const CParticle&, 
				       const CVector<>&, 
				       const CVector<>&
				       ) const;

  virtual size_t getSquareCellCollision3(const CParticle&, 
				       const CVector<>&, 
				       const CVector<>&
				       ) const;
  
  virtual C2ParticleData SmoothSpheresColl(const CIntEvent&, const Iflt&, const Iflt&, const EEventType& eType) const;

  virtual bool DSMCSpheresTest(const CParticle&, const CParticle&, 
			       Iflt&, const Iflt&, CPDData&) const;

  virtual C2ParticleData DSMCSpheresRun(const CParticle&, const CParticle&, 
					const Iflt&, CPDData&) const;
  
  virtual C2ParticleData SphereWellEvent(const CIntEvent&, const Iflt&, const Iflt&) const;

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

  virtual CLiouvillean* Clone() const { return new CLSLLOD(*this); }

protected:
  virtual void outputXML(xmlw::XmlStream& ) const;
};
#endif
