/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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

#include <dynamo/dynamics/newtonian.hpp>

namespace dynamo {
  /*! \brief A Dynamics which implements standard Newtonian dynamics
    with an additional constant force vector.
  */
  class Dyngravity: public DynNewtonian
  {
  public:
    Dyngravity(dynamo::Simulation*, const magnet::xml::Node&);

    Dyngravity(dynamo::Simulation* tmp, Vector gravity, double eV = 0, double tc = -HUGE_VAL);

    void initialise();

    const Vector& getGravityVector() const { return g; }

    virtual double SphereSphereInRoot(const Particle& p1, const Particle& p2, double d) const;
    virtual double SphereSphereInRoot(const Range& p1, const Range& p2, double d) const;
    virtual double SphereSphereOutRoot(const Particle& p1, const Particle& p2, double d) const;
    virtual double SphereSphereOutRoot(const Range& p1, const Range& p2, double d) const;

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

    virtual double getParabolaSentinelTime(const Particle&) const;

    virtual void enforceParabola(Particle&) const;

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

    virtual std::pair<double, Dynamics::TriangleIntersectingPart> 
    getSphereTriangleEvent(const Particle& part, 
			   const Vector & A, 
			   const Vector & B, 
			   const Vector & C,
			   const double dist
			   ) const;

    virtual ParticleEventData runWallCollision(Particle&, 
					       const Vector &,
					       const double&
					       ) const;
  protected:
    double elasticV;
    Vector g;

    mutable std::vector<long double> _tcList;
    double _tc;

    virtual void outputXML(magnet::xml::XmlStream&) const;
  };
}
