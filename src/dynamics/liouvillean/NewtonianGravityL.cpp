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

#include "NewtonianGravityL.hpp"
#include "../../extcode/xmlwriter.hpp"
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
#include <boost/math/special_functions/fpclassify.hpp>
#include <magnet/math/cubic.hpp>
#include <magnet/math/bisect.hpp>
#include <algorithm>

LNewtonianGravity::LNewtonianGravity(DYNAMO::SimData* tmp, const XMLNode& XML):
  LNewtonian(tmp),
  Gravity(-1),
  GravityDim(1),
  elasticV(0),
  _tc(-HUGE_VAL)
{
  if (strcmp(XML.getAttribute("Type"),"NewtonianGravity"))
    M_throw() << "Attempting to load NewtonianGravity from "
	      << XML.getAttribute("Type")
	      << " entry";
  try 
    {
      if (XML.isAttributeSet("Gravity"))
	Gravity = boost::lexical_cast<double>(XML.getAttribute("Gravity"));      
      
      if (XML.isAttributeSet("GravityDimension"))
	GravityDim = boost::lexical_cast<double>(XML.getAttribute("GravityDimension"));      

      if (XML.isAttributeSet("ElasticV"))
	elasticV = boost::lexical_cast<double>(XML.getAttribute("ElasticV")) 
	  * Sim->dynamics.units().unitVelocity();

      if (XML.isAttributeSet("tc"))
	{
	  _tc = boost::lexical_cast<double>(XML.getAttribute("tc")) 
	    * Sim->dynamics.units().unitTime();

	  if (_tc <= 0) M_throw() << "tc must be positive! (tc = " << _tc/ Sim->dynamics.units().unitTime() << ")";
	}
      
    }
  catch (boost::bad_lexical_cast &)
    {
      M_throw() << "Failed a lexical cast in LNewtonianGravity";
    }
  
  Gravity *= Sim->dynamics.units().unitAcceleration();

  g[GravityDim] = Gravity;
}

LNewtonianGravity::LNewtonianGravity(DYNAMO::SimData* tmp, double gravity, 
				     size_t gravityDim, double eV, double tc):
  LNewtonian(tmp), 
  Gravity(gravity), 
  GravityDim(gravityDim),
  elasticV(eV),
  _tc(tc)
{
  g[GravityDim] = Gravity;
}

void
LNewtonianGravity::streamParticle(Particle &particle, const double &dt) const
{
  particle.getPosition() += dt * particle.getVelocity();

  bool isDynamic = particle.testState(Particle::DYNAMIC);
  particle.getPosition()[GravityDim] += 0.5 * dt * dt * Gravity * isDynamic;
  particle.getVelocity()[GravityDim] += dt * Gravity * isDynamic;
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
  if (p1Dynamic == p2Dynamic)
    return LNewtonian::SphereSphereInRoot(dat,d2,p1Dynamic,p2Dynamic);

  Vector gLocal = g;
  if (p2Dynamic) gLocal = -g;

  magnet::math::Bisect<QuarticFunc> quartic;
  quartic.coeffs[0] = 0.25 * Gravity * Gravity;
  quartic.coeffs[1] = gLocal | dat.vij;
  quartic.coeffs[2] = dat.v2 + (gLocal | dat.rij);
  quartic.coeffs[3] = 2 * dat.rvdot;
  quartic.coeffs[4] = dat.r2 - d2;
  
  const double rootthreshold = 1e-16 * sqrt(d2);
  double roots[3];

  //First check if we're already at a root
  if ((dat.p1 != NULL) && (dat.p2 != NULL))
    if (((dat.p1->getID() == lastCollParticle1 && dat.p2->getID() == lastCollParticle2)
	 || (dat.p1->getID() == lastCollParticle2 && dat.p2->getID() == lastCollParticle1))
	&& Sim->dSysTime == lastAbsoluteClock)
      {
	//This collision has already happened, we can factor out this
	//root and use the cubic formula to test for more.
	//\f$A t^4 + B t^3 + C t^2 + D t + 0 == 0\f$ is rewritten as
	//\f$t^3 + (B/A) t^2 + (C/A) t + (D/A) == 0\f$
	size_t rootCount = magnet::math::cubicSolve(quartic.coeffs[1] / quartic.coeffs[0], 
						    quartic.coeffs[2] / quartic.coeffs[0], 
						    quartic.coeffs[3] / quartic.coeffs[0], 
						    roots[0], roots[1], roots[2]);
	
	//If there is just one root, it's the entrance root to the
	//current exit root. (assuming our last collision was an exit!)
	if (rootCount != 3) return false;
	
	//Sort all roots
	std::sort(roots, roots + 3);
	
	//t=0 is either the second or fourth root of the quartic (we
	//just had a collision so t=0 is the exit root).  Check the
	//second of the cubics roots, if it's positive it's the
	//re-entry root
	if (roots[1] > 0)
	  { dat.dt = std::max(0.0, 0.999 * roots[1]); return true; }
	
	//There is the chance that this is the 
	return false;
      }
  
  //We calculate the roots of the cubic differential of F
  //\f$F=A t^4 + B t^3 + C t^2 + D t + E == 0\f$ taking the differential gives
  //\f$F=4 A t^3 + 3 B t^2 + 2C t + D == 0\f$ and normalizing the cubic term gives
  //\f$F=t^3 + \frac{3 B}{4 A} t^2 + \frac{2C}{4 A} t + \frac{D}{4 A} == 0\f$
  size_t rootCount = magnet::math::cubicSolve(quartic.coeffs[1] * 3 / (4 * quartic.coeffs[0]), 
					      quartic.coeffs[2] * 2 / (4 * quartic.coeffs[0]), 
					      quartic.coeffs[3] * 1 / (4 * quartic.coeffs[0]), 
					      roots[0], roots[1], roots[2]);
  
  //Sort the roots in ascending order
  std::sort(roots, roots + rootCount);

  //The roots are (in order) a minimum, (and if we have 3 roots), a
  //local maximum, then another minimum

  //Only accept time-positive minimums (otherwise collision was in the
  //past) and check an overlap actually occurs at the minimum

  //Check the first minimum (we always have one)
  if ((roots[0] >= 0) && (quartic(roots[0]) <= 0))
    {
      dat.dt = std::max(0.0, 0.999 * quartic.bisectRoot(0, roots[0], rootthreshold));
      return true;
    }

  //Check the second minimum if we have one
  if ((rootCount > 1)
      && (roots[2] > 0)
      && (quartic(roots[2]) < 0)
      && (quartic(std::max(0.0, roots[1])) > 0))
    {
      dat.dt = std::max(0.0, 0.999 * quartic.bisectRoot(std::max(0.0, roots[1]), roots[2], rootthreshold));
      return true;
    }
  
  return false;
}
  
bool 
LNewtonianGravity::SphereSphereOutRoot(CPDData& dat, const double& d2, bool p1Dynamic, bool p2Dynamic) const
{
  if (p1Dynamic == p2Dynamic)
    return LNewtonian::SphereSphereOutRoot(dat,d2,p1Dynamic,p2Dynamic);

  M_throw() << "Unsupported";
}


double 
LNewtonianGravity::getWallCollision(const Particle &part, 
				    const Vector  &wallLoc, 
				    const Vector  &wallNorm) const
{
  Vector  rij = part.getPosition(),
    vel = part.getVelocity();

  Sim->dynamics.BCs().applyBC(rij, vel);

  double adot = wallNorm[GravityDim] * Gravity * part.testState(Particle::DYNAMIC);
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
    if ((iDim == GravityDim) && part.testState(Particle::DYNAMIC))
      {
	double adot = Gravity;
	double vdot = vel[GravityDim];

	//First check the "upper" boundary that may have no roots
	double rdot = (Gravity < 0) ? rpos[iDim]-width[iDim] : rpos[iDim];
	double arg = vdot * vdot - 2 * rdot * adot;
	double upperRoot1(HUGE_VAL), upperRoot2(HUGE_VAL);
	
	if (arg >= 0)
	  {
	    double t = -(vdot + ((vdot<0) ? -1: 1) * std::sqrt(arg));
	    upperRoot1 = t / adot;
	    upperRoot1 = 2 * rdot / t;
	    if (upperRoot2 < upperRoot1) std::swap(upperRoot2, upperRoot1);
	  }

	
	//Now the lower boundary which always has roots
	rdot = (Gravity < 0) ? rpos[iDim] : rpos[iDim] - width[iDim];
	arg = vdot * vdot - 2 * rdot * adot;
	double lowerRoot1(HUGE_VAL), lowerRoot2(HUGE_VAL);
	if (arg >= 0)
	  {
	    double t = -(vdot + ((vdot<0) ? -1: 1) * std::sqrt(arg));
	    lowerRoot1 = t / adot;
	    lowerRoot2 = 2 * rdot / t;
	    if (lowerRoot2 < lowerRoot1) std::swap(lowerRoot2, lowerRoot1);
	  }

	double root = HUGE_VAL;
	//Now, if the velocity is "up", and the upper roots exist,
	//then pick the shortest one
	if (!((Gravity < 0) - (vel[GravityDim] > 0))
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
    if ((iDim == GravityDim) && part.testState(Particle::DYNAMIC))
      {
	double adot = Gravity;
	double vdot = vel[GravityDim];
	
	//First check the "upper" boundary that may have no roots
	double rdot = (Gravity < 0) ? rpos[iDim]-width[iDim]: rpos[iDim];
	double arg = vdot * vdot - 2 * rdot * adot;
	double upperRoot1(HUGE_VAL), upperRoot2(HUGE_VAL);
	
	if (arg >= 0)
	  {
	    double t = -(vdot + ((vdot<0) ? -1: 1) * std::sqrt(arg));
	    upperRoot1 = t / adot;
	    upperRoot1 = 2 * rdot / t;
	    if (upperRoot2 < upperRoot1) std::swap(upperRoot2, upperRoot1);
	  }

	
	//Now the lower boundary which always has roots
	rdot = (Gravity < 0) ? rpos[iDim] : rpos[iDim]-width[iDim];
	arg = vdot * vdot - 2 * rdot * adot;
	double lowerRoot1(HUGE_VAL), lowerRoot2(HUGE_VAL);
	if (arg >= 0)
	  {
	    double t = -(vdot + ((vdot<0) ? -1: 1) * std::sqrt(arg));
	    lowerRoot1 = t / adot;
	    lowerRoot2 = 2 * rdot / t;
	    if (lowerRoot2 < lowerRoot1) std::swap(lowerRoot2, lowerRoot1);
	  }

	//Now, if the velocity is "up", and the upper roots exist,
	//then pick the shortest one
	if (!((Gravity < 0) - (vel[GravityDim] > 0)))
	  if (upperRoot1 < time)
	    {
	      time = upperRoot1;
	      retVal = (Gravity < 0) ? (iDim + 1) : -(iDim + 1);
	    }

	//Otherwise its usually the latest lowerRoot
	if (lowerRoot2 < time)
	  {
	    time = lowerRoot2;
	    retVal = (Gravity < 0) ? - (iDim + 1) : (iDim + 1);
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
      << "NewtonianGravity"
      << xml::attr("Gravity") 
      << Gravity / Sim->dynamics.units().unitAcceleration()
      << xml::attr("GravityDimension") 
      << GravityDim
    ;

  if (elasticV)
    XML << xml::attr("ElasticV") << elasticV / Sim->dynamics.units().unitVelocity();

  if (_tc > 0)
    XML << xml::attr("tc") << _tc / Sim->dynamics.units().unitTime();
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
    if (i != GravityDim)
      {
	double tmp = (0.5 * Sim->aspectRatio[i] - lMax) / fabs(vel[i]);
	
	if (tmp < retval) retval = tmp;
      }
    else
      {
	double roots[2];
	if (magnet::math::quadSolve((0.5 * Sim->aspectRatio[i] - lMax),
				    vel[i], 0.5 * Gravity, roots[0], roots[1]))
	  {
	    if ((roots[0] > 0) && (roots[0] < retval)) retval = roots[0];
	    if ((roots[1] > 0) && (roots[1] < retval)) retval = roots[1];
	  }

	if (magnet::math::quadSolve(-(0.5 * Sim->aspectRatio[i] - lMax),
				    vel[i], 0.5 * Gravity, roots[0], roots[1]))
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
LNewtonianGravity::getParabolaSentinelTime(const Particle& part, 
					   unsigned char& passed) const
{
#ifdef DYNAMO_DEBUG
  if (!isUpToDate(part))
    M_throw() << "Particle is not up to date";
#endif
  
  Vector pos(part.getPosition()), vel(part.getVelocity());
  
  Sim->dynamics.BCs().applyBC(pos, vel);
  
  double turningPoint = - vel[GravityDim] / Gravity;
  
  if (turningPoint <= 0)
    {
      passed = true;
      return HUGE_VAL;
    }
  
  return turningPoint;
}

void 
LNewtonianGravity::enforceParabola(const Particle& part) const
{
  updateParticle(part);

  const_cast<Particle&>(part).getVelocity()[GravityDim] = 0.0;
}
