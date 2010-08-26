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

#pragma once

#include "NewtonL.hpp"

class LNewtonianGravity: public LNewtonian
{
public:
  LNewtonianGravity(DYNAMO::SimData*, const XMLNode&);

  //Pair particle dynamics

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

  virtual Iflt getPBCSentinelTime(const Particle&, const Iflt&) const;

  virtual Iflt getParabolaSentinelTime(const Particle&, unsigned char&) const;

  virtual void enforceParabola(const Particle&) const;

  virtual Iflt getWallCollision(const Particle&, 
				const Vector &, 
				const Vector &) const;

  virtual Iflt getCylinderWallCollision(const Particle&, 
					const Vector &, 
					const Vector &,
					const Iflt&
					) const;

  //Cloning
  virtual Liouvillean* Clone() const { return new LNewtonianGravity(*this); }

  size_t getGravityDimension() const { return GravityDim; }
  Iflt getGravity() const { return Gravity; }

protected:
  Iflt Gravity;
  size_t GravityDim;

  virtual void outputXML(xml::XmlStream&) const;
};
