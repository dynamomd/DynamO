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
  class DynCompression: public DynNewtonian
  {
  public:
    DynCompression(dynamo::Simulation*, double);

    virtual double SphereSphereInRoot(const Particle& p1, const Particle& p2, double d) const;
    virtual double SphereSphereOutRoot(const Particle& p1, const Particle& p2, double d) const;  
    virtual double sphereOverlap(const Particle& p1, const Particle& p2,
				 const double& d) const;

    virtual PairEventData SmoothSpheresColl(const IntEvent&, const double&, const double&, const EEventType&) const;

    virtual PairEventData SphereWellEvent(const IntEvent&, const double&, const double&) const;
  
    double getGrowthRate() const { return growthRate; }

    virtual double getPBCSentinelTime(const Particle&, const double&) const;

    virtual double CubeCubeInRoot(const Particle& p1, const Particle& p2, 
				  double d) const 
    { M_throw() << "Not Implemented"; }

    virtual bool cubeOverlap(const Particle& p1, const Particle& p2, 
			     const double d) const 
    { M_throw() << "Not Implemented"; }

    virtual PairEventData parallelCubeColl(const IntEvent& event,
					   const double& e, const double& d,
					   const EEventType& eType = CORE) const;
  
  protected:
    virtual void outputXML(magnet::xml::XmlStream& ) const;
  
    double growthRate;
  };
}
