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

#ifndef LCompression_H
#define LCompression_H

#include "NewtonL.hpp"

class LCompression: public LNewtonian
{
public:
  LCompression(DYNAMO::SimData*, double);

  virtual bool SphereSphereInRoot(CPDData&, const double&, bool p1Dynamic, bool p2Dynamic) const;
  virtual bool SphereSphereOutRoot(CPDData&, const double&, bool p1Dynamic, bool p2Dynamic) const;  
  virtual bool sphereOverlap(const CPDData&, const double&) const;

  virtual void streamParticle(Particle&, const double&) const;

  virtual PairEventData SmoothSpheresColl(const IntEvent&, const double&, const double&, const EEventType&) const;

  virtual PairEventData SphereWellEvent(const IntEvent&, const double&, const double&) const;
  
  virtual Liouvillean* Clone() const { return new LCompression(*this); };

  double getGrowthRate() const { return growthRate; }

  virtual double getPBCSentinelTime(const Particle&, const double&) const;

  virtual bool CubeCubeInRoot(CPDData& pd, const double& d) const { M_throw() << "Not Implemented"; }

  virtual bool CubeCubeOutRoot(CPDData&, const double& d) const { M_throw() << "Not Implemented"; }

  virtual bool cubeOverlap(const CPDData& PD, const double& d) const { M_throw() << "Not Implemented"; }

  virtual PairEventData parallelCubeColl(const IntEvent& event, 
					  const double& e, 
					  const double& d, 
					  const EEventType& eType = CORE
					  ) const;
  
protected:
  virtual void outputXML(xml::XmlStream& ) const;
  
  double growthRate;
};
#endif
