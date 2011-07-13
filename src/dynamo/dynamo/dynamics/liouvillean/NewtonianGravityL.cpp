/*  dynamo:- Event driven molecular dynamics simulator 
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

#include "NewtonianGravityL.hpp"
#include "../interactions/intEvent.hpp"
#include "../2particleEventData.hpp"
#include "../NparticleEventData.hpp"
#include "../dynamics.hpp"
#include "../BC/BC.hpp"
#include "../../base/is_simdata.hpp"
#include "../species/species.hpp"
#include "../../schedulers/sorters/datastruct.hpp"
#include "shapes/frenkelroot.hpp"
#include "shapes/oscillatingplate.hpp"
#include <dynamo/datatypes/vector.xml.hpp>
#include <dynamo/dynamics/globals/neighbourList.hpp>
#include <dynamo/dynamics/globals/ParabolaSentinel.hpp>
#include <dynamo/dynamics/units/units.hpp>
#include <magnet/overlap/point_prism.hpp>
#include <magnet/intersection/parabola_sphere.hpp>
#include <magnet/intersection/parabola_plane.hpp>
#include <magnet/intersection/parabola_triangle.hpp>
#include <magnet/intersection/parabola_rod.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <algorithm>



LNewtonianGravity::LNewtonianGravity(dynamo::SimData* tmp, const magnet::xml::Node& XML):
  LNewtonian(tmp),
  elasticV(0),
  g(0, -1, 0),
  _tc(-HUGE_VAL)
{
  if (strcmp(XML.getAttribute("Type"), "NewtonianGravity"))
    M_throw() << "Attempting to load NewtonianGravity from "
	      << XML.getAttribute("Type")
	      << " entry";
  try 
    {
      if (XML.hasAttribute("ElasticV"))
	elasticV = XML.getAttribute("ElasticV").as<double>()
	  * Sim->dynamics.units().unitVelocity();

      if (XML.hasAttribute("tc"))
	{
	  _tc = XML.getAttribute("tc").as<double>() * Sim->dynamics.units().unitTime();

	  if (_tc <= 0) 
	    M_throw() << "tc must be positive! (tc = " 
		      << _tc/ Sim->dynamics.units().unitTime() << ")";
	}

      g << XML.getNode("g");
    }
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in LNewtonianGravity";
    }
  
  g *= Sim->dynamics.units().unitAcceleration();
}

LNewtonianGravity::LNewtonianGravity(dynamo::SimData* tmp, Vector gravity, double eV, double tc):
  LNewtonian(tmp), 
  elasticV(eV),
  g(gravity),
  _tc(tc)
{}

void
LNewtonianGravity::streamParticle(Particle &particle, const double &dt) const
{
  bool isDynamic = particle.testState(Particle::DYNAMIC);
  particle.getPosition() += dt * (particle.getVelocity() + 0.5 * dt * g * isDynamic);
  particle.getVelocity() += dt * g * isDynamic;
}

bool 
LNewtonianGravity::SphereSphereInRoot(CPDData& dat, const double& d2, 
				      bool p1Dynamic, bool p2Dynamic) const
{
  //If both particles feel gravity, or both don't, the root finding is the same.
  if (p1Dynamic == p2Dynamic)
    return LNewtonian::SphereSphereInRoot(dat,d2,p1Dynamic,p2Dynamic);

  //One particle feels gravity and the other does not. Here we get the
  //sign right on the acceleration g12
  Vector gij(g);
  if (p2Dynamic) gij = -g;

  //Now test for a parabolic ray and sphere intersection
  dat.dt = magnet::intersection::parabola_sphere_bfc(dat.rij, dat.vij, gij, std::sqrt(d2));
    
  return (dat.dt != HUGE_VAL);
}
  
bool 
LNewtonianGravity::SphereSphereOutRoot(CPDData& dat, const double& d2, bool p1Dynamic, bool p2Dynamic) const
{
  if (p1Dynamic == p2Dynamic)
    return LNewtonian::SphereSphereOutRoot(dat,d2,p1Dynamic,p2Dynamic);

  M_throw() << "Function not implemented";
}


double 
LNewtonianGravity::getWallCollision(const Particle &part, 
				    const Vector  &wallLoc, 
				    const Vector  &wallNorm) const
{
  Vector  rij = part.getPosition() - wallLoc,
    vij = part.getVelocity();

  Sim->dynamics.BCs().applyBC(rij, vij);

  return magnet::intersection::parabola_plane_bfc(rij, vij, g * part.testState(Particle::DYNAMIC), wallNorm);
}

double
LNewtonianGravity::getSquareCellCollision2(const Particle& part, 
					   const Vector & origin, 
					   const Vector & width) const
{
  Vector rpos(part.getPosition() - origin);
  Vector vel(part.getVelocity());
  Sim->dynamics.BCs().applyBC(rpos, vel);
  
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
LNewtonianGravity::getSquareCellCollision3(const Particle& part, 
					   const Vector & origin, 
					   const Vector & width) const
{
  Vector  rpos(part.getPosition() - origin);
  Vector  vel(part.getVelocity());

  Sim->dynamics.BCs().applyBC(rpos, vel);

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
LNewtonianGravity::outputXML(magnet::xml::XmlStream& XML) const
{
  XML << magnet::xml::attr("Type") 
      << "NewtonianGravity";

  if (elasticV)
    XML << magnet::xml::attr("ElasticV") << elasticV / Sim->dynamics.units().unitVelocity();

  if (_tc > 0)
    XML << magnet::xml::attr("tc") << _tc / Sim->dynamics.units().unitTime();

  XML << magnet::xml::tag("g") << g / Sim->dynamics.units().unitAcceleration() << magnet::xml::endtag("g");
}

double 
LNewtonianGravity::getPBCSentinelTime(const Particle& part, const double& lMax) const
{
#ifdef DYNAMO_DEBUG
  if (!isUpToDate(part))
    M_throw() << "Particle is not up to date";
#endif

  if (!part.testState(Particle::DYNAMIC)) return LNewtonian::getPBCSentinelTime(part, lMax);

  Vector pos(part.getPosition()), vel(part.getVelocity());

  Sim->dynamics.BCs().applyBC(pos, vel);

  double retval = HUGE_VAL;

  for (size_t i(0); i < NDIM; ++i)
    if (g[i] == 0)
      {
	double tmp = (0.5 * Sim->primaryCellSize[i] - lMax) / fabs(vel[i]);
	
	if (tmp < retval) retval = tmp;
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
LNewtonianGravity::getPointPlateCollision(const Particle& part, const Vector& nrw0,
				 const Vector& nhat, const double& Delta,
				 const double& Omega, const double& Sigma,
				 const double& t, bool lastpart) const
{
  M_throw() << "Not implemented yet";
}

void
LNewtonianGravity::initialise()
{
  if (_tc > 0) _tcList.resize(Sim->N, -HUGE_VAL);
  LNewtonian::initialise();

  //Now add the parabola sentinel if there are cell neighbor lists in
  //use.
  bool hasNBlist = false;
  BOOST_FOREACH(magnet::ClonePtr<Global>& glob, Sim->dynamics.getGlobals())
    if (glob.typeTest<GNeighbourList>())
      { hasNBlist = true; break; }
  
  if (hasNBlist)
    Sim->dynamics.addGlobal(new CGParabolaSentinel(Sim, "NBListParabolaSentinel"));
}

PairEventData 
LNewtonianGravity::SmoothSpheresColl(const IntEvent& event, const double& ne,
				     const double& d2, const EEventType& eType) const
{
  const Particle& particle1 = Sim->particleList[event.getParticle1ID()];
  const Particle& particle2 = Sim->particleList[event.getParticle2ID()];

  updateParticlePair(particle1, particle2);  

  Vector rij = particle1.getPosition() - particle2.getPosition(),
    vij = particle1.getVelocity() - particle2.getVelocity();

  Sim->dynamics.BCs().applyBC(rij, vij);


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

  return LNewtonian::SmoothSpheresColl(event, e, d2, eType);
}


double 
LNewtonianGravity::getCylinderWallCollision(const Particle& part, 
				   const Vector& wallLoc, 
				   const Vector& wallNorm,
				   const double& radius) const
{
  M_throw() << "Not implemented yet";
}

double 
LNewtonianGravity::getParabolaSentinelTime(const Particle& part) const
{
#ifdef DYNAMO_DEBUG
  if (!isUpToDate(part))
    M_throw() << "Particle is not up to date";
#endif
  
  if (!part.testState(Particle::DYNAMIC)) return HUGE_VAL; //Particle is not dynamic (does not feel gravity)

  Vector pos(part.getPosition()), vel(part.getVelocity());  
  Sim->dynamics.BCs().applyBC(pos, vel);
  
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
LNewtonianGravity::enforceParabola(const Particle& part) const
{
  updateParticle(part);
  Vector pos(part.getPosition()), vel(part.getVelocity());
  Sim->dynamics.BCs().applyBC(pos, vel);

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

  const_cast<Particle&>(part).getVelocity()[dim] = 0;
}

std::pair<double, Liouvillean::TriangleIntersectingPart>
LNewtonianGravity::getSphereTriangleEvent(const Particle& part, 
					  const Vector & A, 
					  const Vector & B, 
					  const Vector & C,
					  const double dist
					  ) const
{
  //If the particle doesn't feel gravity, fall back to the standard function
  if (!part.testState(Particle::DYNAMIC)) 
    return LNewtonian::getSphereTriangleEvent(part, A, B, C, dist);

  typedef std::pair<double, Liouvillean::TriangleIntersectingPart> RetType;
  //The Origin, relative to the first vertex
  Vector T = part.getPosition() - A;
  //The ray direction
  Vector D = part.getVelocity();
  Sim->dynamics.BCs().applyBC(T, D);

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
LNewtonianGravity::runWallCollision(const Particle &part, 
				    const Vector  &vNorm,
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

  return LNewtonian::runWallCollision(part, vNorm, e_val);
}
