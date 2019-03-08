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

namespace dynamo {
  DynViscous::DynViscous(dynamo::Simulation* tmp, const magnet::xml::Node& XML):
    DynNewtonian(tmp), g({0, -1, 0}), lastAbsoluteClock(-1), lastCollParticle1(0), lastCollParticle2(0)
  {
    g << XML.getNode("g");
    g *= Sim->units.unitAcceleration();
  }

  void 
  DynViscous::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") << "Viscous";
    XML << magnet::xml::tag("g") << g / Sim->units.unitAcceleration() << magnet::xml::endtag("g");
  }

  void
  DynViscous::streamParticle(Particle &particle, const double &dt) const
  {
    particle.getPosition() += particle.getVelocity() * dt;

    if (hasOrientationData())
      {
	orientationData[particle.getID()].orientation = Quaternion::fromRotationAxis(orientationData[particle.getID()].angularVelocity * dt)
	  * orientationData[particle.getID()].orientation ;
	orientationData[particle.getID()].orientation.normalise();
      }
  }

  double
  DynViscous::SphereSphereInRoot(const Particle& p1, const Particle& p2, double d) const
  {
    Vector r12 = p1.getPosition() - p2.getPosition();
    Vector v12 = p1.getVelocity() - p2.getVelocity();
    Sim->BCs->applyBC(r12, v12);
    return magnet::intersection::ray_sphere(r12, v12, d);
  }
  
  double 
  DynViscous::getPBCSentinelTime(const Particle& part, const double& lMax) const
  {
#ifdef DYNAMO_DEBUG
    if (!isUpToDate(part))
      M_throw() << "Particle is not up to date";
#endif

    Vector pos(part.getPosition()), vel(part.getVelocity());

    Sim->BCs->applyBC(pos, vel);

    double retval = std::numeric_limits<float>::infinity();

    for (size_t i(0); i < NDIM; ++i)
      if (vel[i] != 0)
	{
	  double tmp = (0.5 * (0.5 * Sim->primaryCellSize[i] - lMax)) / std::abs(vel[i]);
	  
	  if (tmp < retval)
	    retval = tmp;
	}

    return retval;
  }
}
