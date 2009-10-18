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

#ifndef CLCompression_H
#define CLCompression_H

#include "NewtonL.hpp"

class CLCompression: public CLNewton
{
public:
  CLCompression(DYNAMO::SimData*, Iflt);

  virtual bool SphereSphereInRoot(CPDData&, const Iflt&) const;
  virtual bool SphereSphereOutRoot(CPDData&, const Iflt&) const;  
  virtual bool sphereOverlap(const CPDData&, const Iflt&) const;

  virtual void streamParticle(CParticle&, const Iflt&) const;

  virtual C2ParticleData SmoothSpheresColl(const CIntEvent&, const Iflt&, const Iflt&, const EEventType&) const;

  virtual C2ParticleData SphereWellEvent(const CIntEvent&, const Iflt&, const Iflt&) const;
  
  virtual Liouvillean* Clone() const { return new CLCompression(*this); };

  Iflt getGrowthRate() const { return growthRate; }

  virtual Iflt getPBCSentinelTime(const CParticle&, const Iflt&) const;

  virtual bool CubeCubeInRoot(CPDData& pd, const Iflt& d) const { D_throw() << "Not Implemented"; }

  virtual bool CubeCubeOutRoot(CPDData&, const Iflt& d) const { D_throw() << "Not Implemented"; }

  virtual bool cubeOverlap(const CPDData& PD, const Iflt& d) const { D_throw() << "Not Implemented"; }

  virtual C2ParticleData parallelCubeColl(const CIntEvent& event, 
					  const Iflt& e, 
					  const Iflt& d, 
					  const EEventType& eType = CORE
					  ) const;
  
protected:
  virtual void outputXML(xmlw::XmlStream& ) const;
  
  Iflt growthRate;
};
#endif
