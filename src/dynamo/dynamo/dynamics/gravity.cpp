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

#include <dynamo/dynamics/gravity.hpp>
#include <dynamo/interactions/intEvent.hpp>
#include <dynamo/2particleEventData.hpp>
#include <dynamo/NparticleEventData.hpp>

#include <dynamo/BC/BC.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/species/species.hpp>
#include <dynamo/schedulers/sorters/event.hpp>
#include <dynamo/dynamics/shapes/oscillatingplate.hpp>
#include <dynamo/globals/neighbourList.hpp>
#include <dynamo/globals/ParabolaSentinel.hpp>
#include <dynamo/units/units.hpp>
#include <magnet/overlap/point_prism.hpp>
#include <magnet/intersection/parabola_sphere.hpp>
#include <magnet/intersection/parabola_plane.hpp>
#include <magnet/intersection/parabola_triangle.hpp>
#include <magnet/intersection/parabola_rod.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <algorithm>

namespace dynamo {
  DynGravity::DynGravity(dynamo::Simulation* tmp, const magnet::xml::Node& XML):
    DynNewtonian(tmp),
    elasticV(0),
    g(0, -1, 0),
    _tc(-HUGE_VAL)
  {
    if (strcmp(XML.getAttribute("Type"), "NewtonianGravity"))
      M_throw() << "Attempting to load gravity from "
		<< XML.getAttribute("Type")
		<< " entry";
    try 
      {
	if (XML.hasAttribute("ElasticV"))
	  elasticV = XML.getAttribute("ElasticV").as<double>()
	    * Sim->units.unitVelocity();

	if (XML.hasAttribute("tc"))
	  {
	    _tc = XML.getAttribute("tc").as<double>() * Sim->units.unitTime();

	    if (_tc <= 0) 
	      M_throw() << "tc must be positive! (tc = " 
			<< _tc/ Sim->units.unitTime() << ")";
	  }

	g << XML.getNode("g");
      }
    catch (boost::bad_lexical_cast &)
      {
	M_throw() << "Failed a lexical cast in DynGravity";
      }
  
    g *= Sim->units.unitAcceleration();
  }

  DynGravity::DynGravity(dynamo::Simulation* tmp, Vector gravity, double eV, double tc):
    DynNewtonian(tmp), 
    elasticV(eV),
    g(gravity),
    _tc(tc)
  {}

  void
  DynGravity::streamParticle(Particle &particle, const double &dt) const
  {
    bool isDynamic = particle.testState(Particle::DYNAMIC);
    particle.getPosition() += dt * (particle.getVelocity() + 0.5 * dt * g * isDynamic);
    particle.getVelocity() += dt * g * isDynamic;
  }

  double
  DynGravity::SphereSphereInRoot(const Particle& p1, const Particle& p2, double d) const
  {
    bool p1Dynamic = p1.testState(Particle::DYNAMIC);
    bool p2Dynamic = p2.testState(Particle::DYNAMIC);
    
    //If both particles feel gravity, or both don't, the root finding is the same.
    if (p1Dynamic == p2Dynamic)
      return DynNewtonian::SphereSphereInRoot(p1, p2, d);

    Vector r12 = p1.getPosition() - p2.getPosition();
    Vector v12 = p1.getVelocity() - p2.getVelocity();
    Sim->BCs->applyBC(r12, v12);

    //One particle feels gravity and the other does not. Here we get the
    //sign right on the acceleration g12
    Vector g12(g);
    if (p2Dynamic) g12 = -g;

    //Now test for a parabolic ray and sphere intersection
    return magnet::intersection::parabola_sphere_bfc(r12, v12, g12, d);
  }

  double
  DynGravity::SphereSphereInRoot(const Range& p1, const Range& p2, double d) const
  {
    double accel1sum = 0;
    double mass1 = 0;
    BOOST_FOREACH(const size_t ID, p1)
      {
	const Particle& part = Sim->particleList[ID];       
	double mass = Sim->species[part].getMass(ID);

	if (part.testState(Particle::DYNAMIC))
	  accel1sum += mass;

	mass1 += mass;
      }
    accel1sum /= mass1;
    
    double accel2sum = 0;
    double mass2 = 0;
    BOOST_FOREACH(const size_t ID, p2)
      {
	const Particle& part = Sim->particleList[ID];       
	double mass = Sim->species[part].getMass(ID);

	if (part.testState(Particle::DYNAMIC))
	  accel2sum += mass;

	mass2 += mass;
      }
    accel2sum /= mass2;

    std::pair<Vector, Vector> r1data = getCOMPosVel(p1);
    std::pair<Vector, Vector> r2data = getCOMPosVel(p2);
    Vector r12 = r1data.first - r2data.first;
    Vector v12 = r1data.second - r2data.second;
    Vector a12 = g * (accel1sum - accel2sum);
    Sim->BCs->applyBC(r12, v12);
    return magnet::intersection::parabola_sphere_bfc(r12, v12, a12, d);
  }
  
  double
  DynGravity::SphereSphereOutRoot(const Particle& p1, const Particle& p2, double d) const
  {
    bool p1Dynamic = p1.testState(Particle::DYNAMIC);
    bool p2Dynamic = p2.testState(Particle::DYNAMIC);

    if (p1Dynamic == p2Dynamic)
      return DynNewtonian::SphereSphereOutRoot(p1, p2, d);

    Vector r12 = p1.getPosition() - p2.getPosition();
    Vector v12 = p1.getVelocity() - p2.getVelocity();
    Sim->BCs->applyBC(r12, v12);

    //One particle feels gravity and the other does not. Here we get the
    //sign right on the acceleration g12
    Vector g12(g);
    if (p2Dynamic) g12 = -g;

    //return magnet::intersection::parabola_inv_sphere_bfc(r12, v12, a12, d);
    M_throw() << "Function not implemented";
  }

  double
  DynGravity::SphereSphereOutRoot(const Range& p1, const Range& p2, double d) const
  {
    double accel1sum = 0;
    double mass1 = 0;
    BOOST_FOREACH(const size_t ID, p1)
      {
	const Particle& part = Sim->particleList[ID];       
	double mass = Sim->species[part].getMass(ID);

	if (part.testState(Particle::DYNAMIC))
	  accel1sum += mass;

	mass1 += mass;
      }
    accel1sum /= mass1;
    
    double accel2sum = 0;
    double mass2 = 0;
    BOOST_FOREACH(const size_t ID, p2)
      {
	const Particle& part = Sim->particleList[ID];       
	double mass = Sim->species[part].getMass(ID);

	if (part.testState(Particle::DYNAMIC))
	  accel2sum += mass;

	mass2 += mass;
      }
    accel2sum /= mass2;

    std::pair<Vector, Vector> r1data = getCOMPosVel(p1);
    std::pair<Vector, Vector> r2data = getCOMPosVel(p2);
    Vector r12 = r1data.first - r2data.first;
    Vector v12 = r1data.second - r2data.second;
    Vector a12 = g * (accel1sum - accel2sum);
    Sim->BCs->applyBC(r12, v12);

    M_throw() << "Not implemented yet";
    //return magnet::intersection::parabola_inv_sphere_bfc(r12, v12, a12, d);
  }


  double 
  DynGravity::getWallCollision(const Particle &part, 
				      const Vector  &wallLoc, 
				      const Vector  &wallNorm) const
  {
    Vector  rij = part.getPosition() - wallLoc,
      vij = part.getVelocity();

    Sim->BCs->applyBC(rij, vij);

    return magnet::intersection::parabola_plane_bfc(rij, vij, g * part.testState(Particle::DYNAMIC), wallNorm);
  }

  double
  DynGravity::getSquareCellCollision2(const Particle& part, 
					     const Vector & origin, 
					     const Vector & width) const
  {
    Vector rpos(part.getPosition() - origin);
    Vector vel(part.getVelocity());
    Sim->BCs->applyBC(rpos, vel);
  
#ifdef DYNAMO_DEBUG
    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      if ((vel[iDim] == 0) && (std::signbit(vel[iDim])))
	M_throw() << "You have negative zero velocities, dont use them."
		  << "\nPlease think of the neighbour lists.";
#endif 

    double retVal = HUGE_VAL;

    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      if ((g[iDim] != 0) && part.testState(Particle::DYNAMIC))
	{
	  //First check the "upper" boundary that may have no roots
	  double r = (g[iDim] < 0) ? rpos[iDim] - width[iDim] : rpos[iDim];
	  double arg = vel[iDim] * vel[iDim] - 2 * r * g[iDim];
	  double upperRoot1(HUGE_VAL), upperRoot2(HUGE_VAL);
	
	  if (arg >= 0)
	    {
	      double t = -(vel[iDim] + ((vel[iDim]<0) ? -1: 1) * std::sqrt(arg));
	      upperRoot1 = t / g[iDim];
	      upperRoot2 = 2 * r / t;
	      if (upperRoot2 < upperRoot1) std::swap(upperRoot2, upperRoot1);
	    }

	
	  //Now the lower boundary which always has roots
	  r = (g[iDim] < 0) ? rpos[iDim] : rpos[iDim] - width[iDim];
	  arg = vel[iDim] * vel[iDim] - 2 * r * g[iDim];
	  double lowerRoot1(HUGE_VAL), lowerRoot2(HUGE_VAL);
	  if (arg >= 0)
	    {
	      double t = -(vel[iDim] + ((vel[iDim]<0) ? -1: 1) * std::sqrt(arg));
	      lowerRoot1 = t / g[iDim];
	      lowerRoot2 = 2 * r / t;
	      if (lowerRoot2 < lowerRoot1) std::swap(lowerRoot2, lowerRoot1);
	    }

	  double root = HUGE_VAL;
	  //Now, if the velocity is "up", and the upper roots exist,
	  //then pick the shortest one
	  if (!((g[iDim] < 0) - (vel[iDim] > 0))
	      && (upperRoot1 != HUGE_VAL))
	    root = upperRoot1;

	  //Otherwise its usually the latest lowerRoot
	  if (root == HUGE_VAL)
	    root = lowerRoot2;

	  if (root < retVal)
	    retVal = root;
	}
      else
	{
	  double tmpdt((vel[iDim] < 0)
		       ? -rpos[iDim]/vel[iDim] 
		       : (width[iDim]-rpos[iDim]) / vel[iDim]);
	
	  if (tmpdt < retVal)
	    retVal = tmpdt;
	}
  
    return retVal;
  }

  int
  DynGravity::getSquareCellCollision3(const Particle& part, 
					     const Vector & origin, 
					     const Vector & width) const
  {
    Vector  rpos(part.getPosition() - origin);
    Vector  vel(part.getVelocity());

    Sim->BCs->applyBC(rpos, vel);

    int retVal(0);
    double time(HUGE_VAL);
  
#ifdef DYNAMO_DEBUG
    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      if ((vel[iDim] == 0) && (std::signbit(vel[iDim])))
	M_throw() << "You have negative zero velocities, dont use them."
		  << "\nPlease think of the neighbour lists.";
#endif

    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      if ((g[iDim] != 0) && part.testState(Particle::DYNAMIC))
	{
	  //First check the "upper" boundary that may have no roots
	  double rdot = (g[iDim] < 0) ? rpos[iDim] - width[iDim]: rpos[iDim];
	  double arg = vel[iDim] * vel[iDim] - 2 * rdot * g[iDim];
	  double upperRoot1(HUGE_VAL), upperRoot2(HUGE_VAL);
	
	  if (arg >= 0)
	    {
	      double t = -(vel[iDim] + ((vel[iDim]<0) ? -1: 1) * std::sqrt(arg));
	      upperRoot1 = t / g[iDim];
	      upperRoot2 = 2 * rdot / t;
	      if (upperRoot2 < upperRoot1) std::swap(upperRoot2, upperRoot1);
	    }

	
	  //Now the lower boundary which always has roots
	  rdot = (g[iDim] < 0) ? rpos[iDim] : rpos[iDim]-width[iDim];
	  arg = vel[iDim] * vel[iDim] - 2 * rdot * g[iDim];
	  double lowerRoot1(HUGE_VAL), lowerRoot2(HUGE_VAL);
	  if (arg >= 0)
	    {
	      double t = -(vel[iDim] + ((vel[iDim]<0) ? -1: 1) * std::sqrt(arg));
	      lowerRoot1 = t / g[iDim];
	      lowerRoot2 = 2 * rdot / t;
	      if (lowerRoot2 < lowerRoot1) std::swap(lowerRoot2, lowerRoot1);
	    }

	  //Now, if the velocity is "up", and the upper roots exist,
	  //then pick the shortest one
	  if (!((g[iDim] < 0) - (vel[iDim] > 0)))
	    if (upperRoot1 < time)
	      {
		time = upperRoot1;
		retVal = (g[iDim] < 0) ? (iDim + 1) : -(iDim + 1);
	      }

	  //Otherwise its usually the latest lowerRoot
	  if (lowerRoot2 < time)
	    {
	      time = lowerRoot2;
	      retVal = (g[iDim] < 0) ? - (iDim + 1) : (iDim + 1);
	    }
	}  
      else
	{
	  double tmpdt = ((vel[iDim] < 0) 
			  ? -rpos[iDim] / vel[iDim] 
			  : (width[iDim] - rpos[iDim]) / vel[iDim]);

	  if (tmpdt < time)
	    {
	      time = tmpdt;
	      retVal = (vel[iDim] < 0) ? - (iDim + 1) : iDim + 1;
	    }
	}

    return retVal;
  }

  void 
  DynGravity::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") 
	<< "NewtonianGravity";

    if (elasticV)
      XML << magnet::xml::attr("ElasticV") << elasticV / Sim->units.unitVelocity();

    if (_tc > 0)
      XML << magnet::xml::attr("tc") << _tc / Sim->units.unitTime();

    XML << magnet::xml::tag("g") << g / Sim->units.unitAcceleration() << magnet::xml::endtag("g");
  }

  double 
  DynGravity::getPBCSentinelTime(const Particle& part, const double& lMax) const
  {
#ifdef DYNAMO_DEBUG
    if (!isUpToDate(part))
      M_throw() << "Particle is not up to date";
#endif

    if (!part.testState(Particle::DYNAMIC)) return DynNewtonian::getPBCSentinelTime(part, lMax);

    Vector pos(part.getPosition()), vel(part.getVelocity());

    Sim->BCs->applyBC(pos, vel);

    double retval = HUGE_VAL;

    for (size_t i(0); i < NDIM; ++i)
      if (g[i] == 0)
	{
	  if (vel[i] != 0)
	    {
	      double tmp = (0.5 * Sim->primaryCellSize[i] - lMax) / fabs(vel[i]);
	      
	      if (tmp < retval) retval = tmp;
	    }
	}
      else
	{
	  double roots[2];
	  if (magnet::math::quadSolve((0.5 * Sim->primaryCellSize[i] - lMax),
				      vel[i], 0.5 * g[i], roots[0], roots[1]))
	    {
	      if ((roots[0] > 0) && (roots[0] < retval)) retval = roots[0];
	      if ((roots[1] > 0) && (roots[1] < retval)) retval = roots[1];
	    }

	  if (magnet::math::quadSolve(-(0.5 * Sim->primaryCellSize[i] - lMax),
				      vel[i], 0.5 * g[i], roots[0], roots[1]))
	    {
	      if ((roots[0] > 0) && (roots[0] < retval)) retval = roots[0];
	      if ((roots[1] > 0) && (roots[1] < retval)) retval = roots[1];
	    }
	}

    return retval;
  }

  std::pair<bool,double>
  DynGravity::getPointPlateCollision(const Particle& part, const Vector& nrw0,
					    const Vector& nhat, const double& Delta,
					    const double& Omega, const double& Sigma,
					    const double& t, bool lastpart) const
  {
    M_throw() << "Not implemented yet";
  }

  void
  DynGravity::initialise()
  {
    if (_tc > 0) _tcList.resize(Sim->N, -HUGE_VAL);
    DynNewtonian::initialise();

    //Now add the parabola sentinel if there are cell neighbor lists in
    //use.
    bool hasNBlist = false;
    BOOST_FOREACH(shared_ptr<Global>& glob, Sim->globals)
      if (std::tr1::dynamic_pointer_cast<GNeighbourList>(glob))
	{ hasNBlist = true; break; }
  
    if (hasNBlist)
      Sim->globals.push_back(shared_ptr<Global>(new GParabolaSentinel(Sim, "NBListParabolaSentinel")));
  }

  PairEventData 
  DynGravity::SmoothSpheresColl(const IntEvent& event, const double& ne,
				       const double& d2, const EEventType& eType) const
  {
    Particle& particle1 = Sim->particleList[event.getParticle1ID()];
    Particle& particle2 = Sim->particleList[event.getParticle2ID()];

    updateParticlePair(particle1, particle2);  

    Vector rij = particle1.getPosition() - particle2.getPosition(),
      vij = particle1.getVelocity() - particle2.getVelocity();

    Sim->BCs->applyBC(rij, vij);


    //Check if two particles are collapsing
    //First, the elastic V calculation
    double vnrm = std::fabs((rij | vij) / rij.nrm());  
    double e = ne;
    if (vnrm < elasticV) e = 1.0;
  
    //Check if a particle is collapsing on a static particle
    if (!particle1.testState(Particle::DYNAMIC) 
	|| !particle2.testState(Particle::DYNAMIC))
      {
	double gnrm = g.nrm();
	if (gnrm > 0)
	  if (std::fabs((vij | g) / gnrm) < elasticV) e = 1.0;  
      }
  
    //Now the tc model;
    if (_tc > 0)
      {
	if ((Sim->dSysTime - _tcList[particle1.getID()] < _tc)
	    || (Sim->dSysTime - _tcList[particle2.getID()] < _tc))
	  e = 1.0;
      
	_tcList[particle1.getID()] = Sim->dSysTime;
	_tcList[particle2.getID()] = Sim->dSysTime;
      }

    return DynNewtonian::SmoothSpheresColl(event, e, d2, eType);
  }


  double 
  DynGravity::getCylinderWallCollision(const Particle& part, 
					      const Vector& wallLoc, 
					      const Vector& wallNorm,
					      const double& radius) const
  {
    M_throw() << "Not implemented yet";
  }

  double 
  DynGravity::getParabolaSentinelTime(const Particle& part) const
  {
#ifdef DYNAMO_DEBUG
    if (!isUpToDate(part))
      M_throw() << "Particle is not up to date";
#endif
  
    if (!part.testState(Particle::DYNAMIC)) return HUGE_VAL; //Particle is not dynamic (does not feel gravity)

    Vector pos(part.getPosition()), vel(part.getVelocity());  
    Sim->BCs->applyBC(pos, vel);
  
    //We return the time of the next turning point, per dimension
    double time = HUGE_VAL;
    for (size_t iDim(0); iDim < NDIM; ++iDim)
      if (g[iDim] != 0)
	{
	  double tmpTime = - vel[iDim] / g[iDim];
	  if ((tmpTime > 0) && (tmpTime < time)) time = tmpTime;
	}
  
    return time;
  }

  void 
  DynGravity::enforceParabola(Particle& part) const
  {
    updateParticle(part);
    Vector pos(part.getPosition()), vel(part.getVelocity());
    Sim->BCs->applyBC(pos, vel);

    //Find the dimension that is closest to 
    size_t dim = NDIM;
    double time = HUGE_VAL;
    for (size_t iDim(0); iDim < NDIM; ++iDim)
      if (g[iDim] != 0)
	{
	  double tmpTime = std::abs(- vel[iDim] / g[iDim]);
	  if ((std::abs(tmpTime) < time)) 
	    {
	      time = tmpTime;
	      dim = iDim;
	    }
	}

#ifdef DYNAMO_DEBUG
    if (dim >= NDIM) M_throw() << "Could not find a dimension to enforce the parabola in!";
#endif

    part.getVelocity()[dim] = 0;
  }

  std::pair<double, Dynamics::TriangleIntersectingPart>
  DynGravity::getSphereTriangleEvent(const Particle& part, 
					    const Vector & A, 
					    const Vector & B, 
					    const Vector & C,
					    const double dist
					    ) const
  {
    //If the particle doesn't feel gravity, fall back to the standard function
    if (!part.testState(Particle::DYNAMIC)) 
      return DynNewtonian::getSphereTriangleEvent(part, A, B, C, dist);

    typedef std::pair<double, Dynamics::TriangleIntersectingPart> RetType;
    //The Origin, relative to the first vertex
    Vector T = part.getPosition() - A;
    //The ray direction
    Vector D = part.getVelocity();
    Sim->BCs->applyBC(T, D);

    //The edge vectors
    Vector E1 = B - A;
    Vector E2 = C - A;
  
    Vector N = E1 ^ E2;
    double nrm2 = N.nrm2();
#ifdef DYNAMO_DEBUG
    if (!nrm2) M_throw() << "Degenerate triangle detected!";
#endif
    N /= std::sqrt(nrm2);
  
    //First test for intersections with the triangle faces.
    double t1 = magnet::intersection::parabola_triangle_bfc(T - N * dist, D, g, E1, E2);
    
    if (t1 < 0)
      {
	t1 = HUGE_VAL;
	if ((D | N) > 0)
	  if (magnet::overlap::point_prism(T - N * dist, E1, E2, N, dist)) t1 = 0;
      }

    double t2 = magnet::intersection::parabola_triangle_bfc(T + N * dist, D, g, E2, E1);

    if (t2 < 0)
      {
	t2 = HUGE_VAL;
	if ((D | N) < 0)
	  if (magnet::overlap::point_prism(T + N * dist, E2, E1, -N, dist)) t2 = 0;
      }
  
    RetType retval(std::min(t1, t2), T_FACE);
  
    //Early jump out, to make sure that if we have zero time
    //interactions for the triangle faces, we take them.
    if (retval.first == 0) return retval;
  
    //Now test for intersections with the triangle corners
    double t = magnet::intersection::parabola_sphere_bfc(T, D, g, dist);
    if (t < retval.first) retval = RetType(t, T_A_CORNER);
    t = magnet::intersection::parabola_sphere_bfc(T - E1, D, g, dist);
    if (t < retval.first) retval = RetType(t, T_B_CORNER);
    t = magnet::intersection::parabola_sphere_bfc(T - E2, D, g, dist);
    if (t < retval.first) retval = RetType(t, T_C_CORNER);

    //Now for the edge collision detection
    t = magnet::intersection::parabola_rod_bfc(T, D, g, B - A, dist);
    if (t < retval.first) retval = RetType(t, T_AB_EDGE);
    t = magnet::intersection::parabola_rod_bfc(T, D, g, C - A, dist);
    if (t < retval.first) retval = RetType(t, T_AC_EDGE);
    t = magnet::intersection::parabola_rod_bfc(T - E2, D, g, B - C, dist);
    if (t < retval.first) retval = RetType(t, T_BC_EDGE);

    if (retval.first < 0) retval.first = 0;

    return retval;
  }

  ParticleEventData 
  DynGravity::runWallCollision(Particle &part, 
				      const Vector& vNorm,
				      const double& e
				      ) const
  {
    updateParticle(part);

    double e_val = e;
    if (std::fabs(part.getVelocity() | vNorm) < elasticV) e_val = 1;

    //Now the tc model;
    if (_tc > 0)
      {
	if ((Sim->dSysTime - _tcList[part.getID()] < _tc))
	  e_val = 1;
      
	_tcList[part.getID()] = Sim->dSysTime;
      }

    return DynNewtonian::runWallCollision(part, vNorm, e_val);
  }
}

