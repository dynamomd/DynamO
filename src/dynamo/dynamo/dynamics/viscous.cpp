/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>
    Copyright (C) 2011  Sebastian Gonzalez

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

#include <dynamo/dynamics/viscous.hpp>
#include <dynamo/2particleEventData.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/species/species.hpp>
#include <dynamo/units/units.hpp>
#include <magnet/overlap/point_prism.hpp>
#include <magnet/overlap/point_cube.hpp>
#include <magnet/intersection/ray_triangle.hpp>
#include <magnet/intersection/ray_rod.hpp>
#include <magnet/intersection/ray_sphere.hpp>
#include <magnet/intersection/ray_plane.hpp>
#include <magnet/intersection/ray_cube.hpp>
#include <magnet/intersection/line_line.hpp>
#include <magnet/intersection/overlapfuncs/oscillatingplate.hpp>
#include <magnet/math/matrix.hpp>
#include <magnet/xmlwriter.hpp>

#include <magnet/intersection/polynomial.hpp>

namespace dynamo {
  DynViscous::DynViscous(dynamo::Simulation* tmp, const magnet::xml::Node& XML):
    DynNewtonian(tmp), _g({0, -1, 0})
  {
    _g << XML.getNode("g");
    _g *= Sim->units.unitAcceleration();
    _gamma = XML.getAttribute("gamma").as<double>() * Sim->units.unitAcceleration();
  }
 
  void 
  DynViscous::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "Viscous";
    XML << magnet::xml::attr("gamma") << _gamma * Sim->units.unitTime();
    XML << magnet::xml::tag("g") << _g / Sim->units.unitAcceleration() << magnet::xml::endtag("g");
  }

  void
  DynViscous::streamParticle(Particle &particle, const double &dt) const
  {
    const double m = Sim->species(particle)->getMass(particle.getID());

    particle.getPosition()
      += - _g * dt / (_gamma * m)
      + (1 - std::exp(-_gamma * dt)) * (particle.getVelocity() + _g / _gamma) / (_gamma * m);
      
    particle.getVelocity() = -_g / _gamma + (particle.getVelocity() + _g / _gamma) * std::exp(-_gamma * dt);

    //This part is also incorrect, but left here for future work
    if (hasOrientationData())
      {
	orientationData[particle.getID()].orientation = Quaternion::fromRotationAxis(orientationData[particle.getID()].angularVelocity * dt)
	  * orientationData[particle.getID()].orientation ;
	orientationData[particle.getID()].orientation.normalise();
      }
  }

  double
  DynViscous::SphereSphereInRoot(const Particle& p1, const Particle& p2, double sigma) const
  {
    Vector r12 = p1.getPosition() - p2.getPosition();
    Sim->BCs->applyBC(r12);
    const Vector X = r12;
    const double m1 = Sim->species(p1)->getMass(p1.getID());
    const double m2 = Sim->species(p2)->getMass(p2.getID());
    const Vector V = (p1.getVelocity() / m1) - (p2.getVelocity() / m2);

    if (m1!=m2)
      M_throw() << "Not implemented asymmetric particle masses for viscous dynamics";

    const double c = X.nrm2() - sigma * sigma;
    const double b = 2 * (V | X) / _gamma;
    const double a = V.nrm2() / (_gamma * _gamma);

    //As the transform from y=e^{-\gamma t} to t is a monotonic
    //transform we can perform the stable algorithm in y, but we need
    //to ensure that t=0 corresponds to y=0 and t>0 corresponds to y>0
    // Thus the appropriate transformation is y=1-e^{-\gamma\,t}
    magnet::intersection::detail::PolynomialFunction<2> f(c, b, 2 * a);
    const double y = magnet::intersection::detail::nextEvent(f);

    //If y>1 then there's no real roots
    if (y > 1) return HUGE_VAL;

    const double t = - std::log(1-y) / _gamma;
    return t;
  }
  
  double 
  DynViscous::getPBCSentinelTime(const Particle& part, const double& lMax) const
  {
    //This is bad for low densities!
    return HUGE_VAL;
  }
}
