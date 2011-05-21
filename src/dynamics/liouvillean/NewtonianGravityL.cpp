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
#include "../datatypes/vector.xml.hpp"
#include "../units/units.hpp"
#include <magnet/math/cubic.hpp>
#include <magnet/math/bisect.hpp>
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

namespace {
  struct QuarticFunc
  {
  public:
    inline double operator()(double t)
    {
      return (((coeffs[0] * t + coeffs[1]) * t + coeffs[2]) * t + coeffs[3]) * t + coeffs[4];
    }

    double coeffs[5];
  };
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

  //This is our lengthscale to which we bisect the roots
  const double rootthreshold = 1e-16 * sqrt(d2);

  magnet::math::Bisect<QuarticFunc> quartic;
  quartic.coeffs[0] = 0.25 * gij.nrm2();
  quartic.coeffs[1] = gij | dat.vij;
  quartic.coeffs[2] = dat.v2 + (gij | dat.rij);
  quartic.coeffs[3] = 2 * dat.rvdot;
  quartic.coeffs[4] = dat.r2 - d2;
  
  //We calculate the roots of the cubic differential of F
  //\f$F=A t^4 + B t^3 + C t^2 + D t + E == 0\f$ taking the differential gives
  //\f$F=4 A t^3 + 3 B t^2 + 2C t + D == 0\f$ and normalizing the cubic term gives
  //\f$F=t^3 + \frac{3 B}{4 A} t^2 + \frac{2C}{4 A} t + \frac{D}{4 A} == 0\f$
  double roots[3];
  size_t rootCount = magnet::math::cubicSolve(quartic.coeffs[1] * 3 / (4 * quartic.coeffs[0]), 
					      quartic.coeffs[2] * 2 / (4 * quartic.coeffs[0]), 
					      quartic.coeffs[3] * 1 / (4 * quartic.coeffs[0]), 
					      roots[0], roots[1], roots[2]);

  //Sort the roots in ascending order
  std::sort(roots, roots + rootCount);
  //The roots are (in order) a minimum, (and if we have 3 roots), a
  //local maximum, then another minimum

  //We need to define the dynamics of overlapping particles
  //
  //We follow the convention that, if particles are overlapped and
  //approaching, we collide them now.
  //
  //But, there are some cases where the particle is overlapped,
  //receeding but will not escape the overlap before turning around
  //and approaching again.
  //
  //Here we make another choice, we collide whenever the particle is
  //touching/overlapping and approaching. We may get some kind of
  //inelastic collapse, which is fair, but it may be avoided by
  //sleeping the particles or by other dynamics
  if (quartic(0.0) <= 0)
    {//We're overlapping! We need the overlapped dynamics!
      if (dat.rvdot < 0)
	{
	  //We're overlapping and approaching! Cause an instantaneous collision
	  dat.dt = 0.0;
	  return true;
	}
      
      //We're overlapping but receeding. If there's only one cubic
      //root, then we're fine, we've passed the local minimum (we're
      //receeding) so we'll escape the particle (never to return) soon
      //enough, return no collision
      if (rootCount == 1) return false;

      //We have a local minimum, maximum and minimum. We must be after
      //the first local minimum as we're overlapping but receeding. If
      //the maximum is in the past, and we're receeding, we must be
      //past the second minimum and we will only get further away! so
      //return no collision
      if (roots[1] < 0) return false;

      //The local maximum is in the future (or now). Check if we will be
      //outside the particle then.
      if (quartic(roots[1]) > 0)
	{ //We will be outside the particle, so check if we re-enter it!
	  if (quartic(roots[2]) < 0)
	    {//We do reenter! So return the time of this reentry as the next event
	      dat.dt = std::max(0.0, quartic.bisectRoot(roots[1], roots[2], rootthreshold));
	      return true;
	    }
	  //We wont reenter, so return false as we'll escape
	  return false;
	}
      
      //Our local maxima is in the future, and it's still an
      //overlapped state!  Return an event at the local maximum, it's
      //the furthest away we're going to get without approaching
      //again. Either we'll have another event before we reach it or
      //we'll have to handle it by sleeping the particle.
      dat.dt = roots[1]; //return the local maximum 
      return true;
    }
      
  //We're not overlapping, so check for standard events  

  //Check the first minimum (we always have one), if it's in the
  //future check that we overlap then, and then bisect the event
  if ((roots[0] >= 0) && (quartic(roots[0]) <= 0))
    {
      dat.dt = std::max(0.0, quartic.bisectRoot(0, roots[0], rootthreshold));
      return true;
    }

  //Check the second minimum if we have one. When bisecting we use
  //either now or the local maximum, depending on which is later as
  //the first point.
  if ((rootCount > 1)
      && (roots[2] > 0)
      && (quartic(roots[2]) < 0)
      && (quartic(std::max(0.0, roots[1])) >= 0))
    {
      dat.dt = std::max(0.0, quartic.bisectRoot(std::max(0.0, roots[1]), roots[2], rootthreshold));
      return true;
    }
  
  return false;
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
  Vector  rij = part.getPosition(),
    vel = part.getVelocity();

  Sim->dynamics.BCs().applyBC(rij, vel);

  double adot = (wallNorm | g) * part.testState(Particle::DYNAMIC);
  double vdot = vel | wallNorm;
  double rdot = (rij - wallLoc) | wallNorm;

  double arg = vdot * vdot - 2 * rdot * adot;
  
  if (arg > 0)
    {
      double t = -(vdot + ((vdot<0) ? -1: 1) * std::sqrt(arg));
      double x1 = t / adot;
      double x2 = 2 * rdot / t;

      if (adot > 0)
	//The particle is arcing under the plate
	return (x1 < x2) ? x1 : x2 ;
      else
	//The particle is arcing over the plate
	return (x1 < x2) ? x2 : x1;
    }

  return HUGE_VAL;
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
LNewtonianGravity::outputXML(xml::XmlStream& XML) const
{
  XML << xml::attr("Type") 
      << "NewtonianGravity";

  if (elasticV)
    XML << xml::attr("ElasticV") << elasticV / Sim->dynamics.units().unitVelocity();

  if (_tc > 0)
    XML << xml::attr("tc") << _tc / Sim->dynamics.units().unitTime();

  XML << xml::tag("g") << g / Sim->dynamics.units().unitAcceleration() << xml::endtag("g");
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
