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

#pragma once

#include "NewtonL.hpp"

class LNewtonianGravity: public LNewtonian
{
public:
  LNewtonianGravity(DYNAMO::SimData*, const XMLNode&);

  LNewtonianGravity(DYNAMO::SimData* tmp, double gravity,
		    size_t gravityDim, double eV = 0, double tc = -HUGE_VAL);

  void initialise();

  virtual bool SphereSphereInRoot(CPDData&, const double&, bool p1Dynamic, bool p2Dynamic) const;
  virtual bool SphereSphereOutRoot(CPDData&, const double&, bool p1Dynamic, bool p2Dynamic) const;  

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

  virtual double getPBCSentinelTime(const Particle&, const double&) const;

  virtual double getParabolaSentinelTime(const Particle&, unsigned char&) const;

  virtual void enforceParabola(const Particle&) const;

  virtual double getWallCollision(const Particle&, 
				const Vector &, 
				const Vector &) const;

  virtual double getCylinderWallCollision(const Particle&, 
					  const Vector &, 
					  const Vector &,
					  const double&
					  ) const;

  virtual PairEventData SmoothSpheresColl(const IntEvent&, const double&, 
					  const double&, 
					  const EEventType& eType) const;

  //Cloning
  virtual Liouvillean* Clone() const { return new LNewtonianGravity(*this); }

  size_t getGravityDimension() const { return GravityDim; }
  double getGravity() const { return Gravity; }

protected:
  double Gravity;
  size_t GravityDim;
  double elasticV;
  Vector g;

  mutable std::vector<long double> _tcList;
  double _tc;

  virtual void outputXML(xml::XmlStream&) const;
};
