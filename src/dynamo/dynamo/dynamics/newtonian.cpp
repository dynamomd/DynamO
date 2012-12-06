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

#include <dynamo/dynamics/newtonian.hpp>
#include <dynamo/interactions/intEvent.hpp>
#include <dynamo/2particleEventData.hpp>
#include <dynamo/NparticleEventData.hpp>

#include <dynamo/BC/BC.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/species/species.hpp>
#include <dynamo/schedulers/sorters/event.hpp>
#include <dynamo/dynamics/shapes/oscillatingplate.hpp>
#include <dynamo/dynamics/shapes/lines.hpp>
#include <dynamo/dynamics/shapes/offcentre_sphere.hpp>
#include <dynamo/units/units.hpp>
#include <magnet/math/frenkelroot.hpp>
#include <magnet/overlap/point_prism.hpp>
#include <magnet/overlap/point_cube.hpp>
#include <magnet/intersection/ray_triangle.hpp>
#include <magnet/intersection/ray_rod.hpp>
#include <magnet/intersection/ray_sphere.hpp>
#include <magnet/intersection/ray_plane.hpp>
#include <magnet/intersection/ray_cube.hpp>
#include <magnet/math/matrix.hpp>
#include <magnet/xmlwriter.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

namespace dynamo {
  double
  DynNewtonian::CubeCubeInRoot(const Particle& p1, const Particle& p2, 
			     double d) const 
  {
    Vector r12 = p1.getPosition() - p2.getPosition();
    Vector v12 = p1.getVelocity() - p2.getVelocity();
    Sim->BCs->applyBC(r12, v12);
    return magnet::intersection::ray_AAcube_bfc(r12, v12, 2 * Vector(d, d, d));
  }

  bool 
  DynNewtonian::cubeOverlap(const Particle& p1, const Particle& p2, 
			    const double d) const
  {
    Vector r12 = p1.getPosition() - p2.getPosition();
    Sim->BCs->applyBC(r12);
    return magnet::overlap::point_cube(r12, 2 * Vector(d, d, d));
  }

  double
  DynNewtonian::SphereSphereInRoot(const Particle& p1, const Particle& p2, double d) const
  {
    Vector r12 = p1.getPosition() - p2.getPosition();
    Vector v12 = p1.getVelocity() - p2.getVelocity();
    Sim->BCs->applyBC(r12, v12);
    return magnet::intersection::ray_sphere_bfc(r12, v12, d);
  }

  double
  DynNewtonian::SphereSphereInRoot(const IDRange& p1, const IDRange& p2, double d) const
  {
    std::pair<Vector, Vector> r1data = getCOMPosVel(p1);
    std::pair<Vector, Vector> r2data = getCOMPosVel(p2);
    Vector r12 = r1data.first - r2data.first;
    Vector v12 = r1data.second - r2data.second;
    Sim->BCs->applyBC(r12, v12);
    return magnet::intersection::ray_sphere_bfc(r12, v12, d);
  }
  
  double
  DynNewtonian::SphereSphereOutRoot(const Particle& p1, const Particle& p2, double d) const
  {
    Vector r12 = p1.getPosition() - p2.getPosition();
    Vector v12 = p1.getVelocity() - p2.getVelocity();
    Sim->BCs->applyBC(r12, v12);
    return magnet::intersection::ray_inv_sphere_bfc(r12, v12, d);
  }

  double
  DynNewtonian::SphereSphereOutRoot(const IDRange& p1, const IDRange& p2, double d) const
  {
    std::pair<Vector, Vector> r1data = getCOMPosVel(p1);
    std::pair<Vector, Vector> r2data = getCOMPosVel(p2);
    Vector r12 = r1data.first - r2data.first;
    Vector v12 = r1data.second - r2data.second;
    Sim->BCs->applyBC(r12, v12);
    return magnet::intersection::ray_inv_sphere_bfc(r12, v12, d);
  }

  ParticleEventData 
  DynNewtonian::randomGaussianEvent(Particle& part, const double& sqrtT, 
				  const size_t dimensions) const
  {
#ifdef DYNAMO_DEBUG
    if (dimensions > NDIM)
      M_throw() << "Number of dimensions passed larger than NDIM!";
#endif

    //See http://mathworld.wolfram.com/SpherePointPicking.html

    if (hasOrientationData())
      M_throw() << "Need to implement thermostating of the rotational degrees"
	" of freedom";  

    //Ensure the particle is free streamed first
    updateParticle(part);

    //Collect the precoll data
    ParticleEventData tmpDat(part, *Sim->species[part], GAUSSIAN);

    double mass = Sim->species[tmpDat.getSpeciesID()]->getMass(part.getID());
    double factor = sqrtT / std::sqrt(mass);

    //Assign the new velocities
    for (size_t iDim = 0; iDim < dimensions; iDim++)
      part.getVelocity()[iDim] = Sim->normal_sampler() * factor;

    tmpDat.setDeltaKE(0.5 * mass * (part.getVelocity().nrm2() - tmpDat.getOldVel().nrm2()));
  
    return tmpDat;
  }

  DynNewtonian::DynNewtonian(dynamo::Simulation* tmp):
    Dynamics(tmp),
    lastAbsoluteClock(-1),
    lastCollParticle1(0),
    lastCollParticle2(0)  
  {}

  void
  DynNewtonian::streamParticle(Particle &particle, const double &dt) const
  {
    particle.getPosition() += particle.getVelocity() * dt;

    //The Vector copy is required to make sure that the cached
    //orientation doesn't change during calculation
    if (hasOrientationData())
      orientationData[particle.getID()].orientation 
	= Rodrigues(orientationData[particle.getID()].angularVelocity * dt)
	* Vector(orientationData[particle.getID()].orientation); 
  }

  double 
  DynNewtonian::getPlaneEvent(const Particle& part, const Vector& wallLoc, const Vector& wallNorm, double diameter) const
  {
    Vector rij = part.getPosition() - wallLoc,
      vel = part.getVelocity();

    Sim->BCs->applyBC(rij, vel);

    return magnet::intersection::ray_plane(rij, vel, wallNorm, diameter);
  }

  std::pair<double, Dynamics::TriangleIntersectingPart>
  DynNewtonian::getSphereTriangleEvent(const Particle& part, const Vector & A, const Vector & B, const Vector & C, const double dist) const
  {
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
    double t1 = magnet::intersection::ray_triangle<true, true>(T - N * dist, D, E1, E2);
    
    if (t1 < 0)
      {
	t1 = HUGE_VAL;
	if (magnet::overlap::point_prism(T - N * dist, E1, E2, N, dist)) t1 = 0;
      }

    double t2 = magnet::intersection::ray_triangle<true, true>(T + N * dist, D, E2, E1);

    if (t2 < 0)
      {
	t2 = HUGE_VAL;
	if (magnet::overlap::point_prism(T + N * dist, E2, E1, -N, dist)) t2 = 0;
      }
  
    RetType retval(std::min(t1, t2), T_FACE);
  
    //Early jump out, to make sure that if we have zero time
    //interactions for the triangle faces, we take them.
    if (retval.first == 0) return retval;
  
    //Now test for intersections with the triangle corners
    double t = magnet::intersection::ray_sphere_bfc(T, D, dist);
    if (t < retval.first) retval = RetType(t, T_A_CORNER);
    t = magnet::intersection::ray_sphere_bfc(T - E1, D, dist);
    if (t < retval.first) retval = RetType(t, T_B_CORNER);
    t = magnet::intersection::ray_sphere_bfc(T - E2, D, dist);
    if (t < retval.first) retval = RetType(t, T_C_CORNER);

    //Now for the edge collision detection
    t = magnet::intersection::ray_rod_bfc(T, D, B - A, dist);
    if (t < retval.first) retval = RetType(t, T_AB_EDGE);
    t = magnet::intersection::ray_rod_bfc(T, D, C - A, dist);
    if (t < retval.first) retval = RetType(t, T_AC_EDGE);
    t = magnet::intersection::ray_rod_bfc(T - E2, D, B - C, dist);
    if (t < retval.first) retval = RetType(t, T_BC_EDGE);

    if (retval.first < 0) retval.first = 0;

    return retval;
  }

  ParticleEventData 
  DynNewtonian::runPlaneEvent(Particle& part, const Vector& vNorm, const double e, const double diameter) const
  {
    updateParticle(part);

    ParticleEventData retVal(part, *Sim->species[part], WALL);
  
    part.getVelocity() -= (1+e) * (vNorm | part.getVelocity()) * vNorm;
  
    retVal.setDeltaKE(0.5 * Sim->species[retVal.getSpeciesID()]->getMass(part.getID())
		      * (part.getVelocity().nrm2() - retVal.getOldVel().nrm2()));
    return retVal; 
  }

  ParticleEventData 
  DynNewtonian::runAndersenWallCollision(Particle& part, const Vector & vNorm, const double& sqrtT, const double) const
  {  
    updateParticle(part);

    if (hasOrientationData())
      M_throw() << "Need to implement thermostating of the rotational degrees"
	" of freedom";

    //This gives a completely new random unit vector with a properly
    //distributed Normal component. See Granular Simulation Book
    ParticleEventData tmpDat(part, *Sim->species[part], WALL);
 
    double mass = Sim->species[tmpDat.getSpeciesID()]->getMass(part.getID());

    for (size_t iDim = 0; iDim < NDIM; iDim++)
      part.getVelocity()[iDim] = Sim->normal_sampler() * sqrtT / std::sqrt(mass);
  
    part.getVelocity() 
      //This first line adds a component in the direction of the normal
      += vNorm * (sqrtT * sqrt(-2.0*log(1.0-Sim->uniform_sampler()) / mass)
		  //This removes the original normal component
		  -(part.getVelocity() | vNorm));

    tmpDat.setDeltaKE(0.5 * mass * (part.getVelocity().nrm2() - tmpDat.getOldVel().nrm2()));
  
    return tmpDat; 
  }

  double
  DynNewtonian::getSquareCellCollision2(const Particle& part, const Vector & origin, const Vector & width) const
  {
    Vector  rpos(part.getPosition() - origin);
    Vector  vel(part.getVelocity());
    Sim->BCs->applyBC(rpos, vel);
  
#ifdef DYNAMO_DEBUG
    for (size_t iDim = 0; iDim < NDIM; ++iDim)
      if ((vel[iDim] == 0) && (std::signbit(vel[iDim])))
	M_throw() << "You have negative zero velocities, don't use them.";
#endif 

    double retVal;
    if (vel[0] < 0)
      retVal = -rpos[0] / vel[0];
    else
      retVal = (width[0]-rpos[0]) / vel[0];

    for (size_t iDim = 1; iDim < NDIM; ++iDim)
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
  DynNewtonian::getSquareCellCollision3(const Particle& part, const Vector & origin, const Vector & width) const
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
      {
	double tmpdt = ((vel[iDim] < 0) 
			? -rpos[iDim]/vel[iDim] 
			: (width[iDim]-rpos[iDim]) / vel[iDim]);

	if (tmpdt < time)
	  {
	    time = tmpdt;
	    retVal = (vel[iDim] < 0) ? -(iDim+1) : (iDim+1);
	  }
      }

    if (((retVal < 0) && (vel[abs(retVal)-1] > 0))
	|| ((retVal > 0) && (vel[abs(retVal)-1] < 0)))
      M_throw() << "Found an error! retVal " << retVal
		<< " vel is " << vel[abs(retVal)-1];

    return retVal;
  }

  bool 
  DynNewtonian::DSMCSpheresTest(Particle& p1, Particle& p2, double& maxprob, const double& factor, Vector rij) const
  {
    updateParticlePair(Sim->particles[p1.getID()], Sim->particles[p2.getID()]);

    Vector vij = p1.getVelocity() - p2.getVelocity();
    Sim->BCs->applyBC(rij, vij);

    //Sim->BCs->applyBC(pdat.rij, pdat.vij);
    double rvdot = (rij | vij);
  
    if (rvdot > 0)
      return false; //Positive rvdot

    double prob = factor * (-rvdot);

    if (prob > maxprob)
      maxprob = prob;

    return prob > Sim->uniform_sampler() * maxprob;
  }

  PairEventData
  DynNewtonian::DSMCSpheresRun(Particle& p1, Particle& p2, const double& e, Vector rij) const
  {
    updateParticlePair(Sim->particles[p1.getID()], Sim->particles[p2.getID()]);

    Vector vij = p1.getVelocity() - p2.getVelocity();
    Sim->BCs->applyBC(rij, vij);

    double rvdot = (rij | vij);

    PairEventData retVal(p1, p2,
			 *Sim->species[p1],
			 *Sim->species[p2],
			 CORE);
  
    retVal.rij = rij;
    retVal.rvdot = rvdot;

    double p1Mass = Sim->species[retVal.particle1_.getSpeciesID()]->getMass(p1.getID());
    double p2Mass = Sim->species[retVal.particle2_.getSpeciesID()]->getMass(p2.getID());
    double mu = p1Mass * p2Mass/(p1Mass+p2Mass);

    retVal.dP = rij * ((1.0 + e) * mu * rvdot / rij.nrm2());  

    //This function must edit particles so it overrides the const!
    p1.getVelocity() -= retVal.dP / p1Mass;
    p2.getVelocity() += retVal.dP / p2Mass;

    retVal.particle1_.setDeltaKE(0.5 * p1Mass * (p1.getVelocity().nrm2() 
						 - retVal.particle1_.getOldVel().nrm2()));
  
    retVal.particle2_.setDeltaKE(0.5 * p2Mass * (p2.getVelocity().nrm2() 
						 - retVal.particle2_.getOldVel().nrm2()));
    return retVal;
  }

  PairEventData 
  DynNewtonian::SmoothSpheresColl(const IntEvent& event, const double& e, const double&, const EEventType& eType) const
  {
    Particle& particle1 = Sim->particles[event.getParticle1ID()];
    Particle& particle2 = Sim->particles[event.getParticle2ID()];

    updateParticlePair(particle1, particle2);  

    PairEventData retVal(particle1, particle2,
			 *Sim->species[particle1],
			 *Sim->species[particle2],
			 eType);
    
    Sim->BCs->applyBC(retVal.rij, retVal.vijold);
  
    double p1Mass = Sim->species[retVal.particle1_.getSpeciesID()]->getMass(particle1.getID()); 
    double p2Mass = Sim->species[retVal.particle2_.getSpeciesID()]->getMass(particle2.getID());
 
    retVal.rvdot = (retVal.rij | retVal.vijold);

    //Treat special cases if one particle has infinite mass
    if ((p1Mass == 0) && (p2Mass != 0))
      //if (!particle1.testState(Particle::DYNAMIC) && particle2.testState(Particle::DYNAMIC))
      {
	retVal.dP = p2Mass * retVal.rij * ((1.0 + e) * retVal.rvdot / retVal.rij.nrm2());  
	//This function must edit particles so it overrides the const!
	particle2.getVelocity() += retVal.dP / p2Mass;
      }
    else 
      if ((p1Mass != 0) && (p2Mass == 0))
	//if (particle1.testState(Particle::DYNAMIC) && !particle2.testState(Particle::DYNAMIC))
	{
	  retVal.dP = p1Mass * retVal.rij * ((1.0 + e) * retVal.rvdot / retVal.rij.nrm2());
	  //This function must edit particles so it overrides the const!
	  particle1.getVelocity() -= retVal.dP / p1Mass;
	}
      else
	{
	  bool isInfInf = (p1Mass == 0) && (p2Mass == 0);

	  //If both particles have infinite mass we just collide them as identical masses
	  if (isInfInf) p1Mass = p2Mass = 1;

	  double mu = p1Mass * p2Mass / (p1Mass + p2Mass);

	  retVal.dP = retVal.rij * ((1.0 + e) * mu * retVal.rvdot / retVal.rij.nrm2());  

	  //This function must edit particles so it overrides the const!
	  particle1.getVelocity() -= retVal.dP / p1Mass;
	  particle2.getVelocity() += retVal.dP / p2Mass;

	  //If both particles have infinite mass we pretend no momentum was transferred
	  retVal.dP *= !isInfInf;
	}

    retVal.particle1_.setDeltaKE(0.5 * p1Mass * (particle1.getVelocity().nrm2()
						 - retVal.particle1_.getOldVel().nrm2()));
  
    retVal.particle2_.setDeltaKE(0.5 * p2Mass * (particle2.getVelocity().nrm2() 
						 - retVal.particle2_.getOldVel().nrm2()));

    lastCollParticle1 = particle1.getID();
    lastCollParticle2 = particle2.getID();
    lastAbsoluteClock = Sim->systemTime;

    return retVal;
  }

  PairEventData 
  DynNewtonian::parallelCubeColl(const IntEvent& event, const double& e, const double&, const EEventType& eType) const
  {
    Particle& particle1 = Sim->particles[event.getParticle1ID()];
    Particle& particle2 = Sim->particles[event.getParticle2ID()];

    updateParticlePair(particle1, particle2);

    PairEventData retVal(particle1, particle2,
			 *Sim->species[particle1],
			 *Sim->species[particle2],
			 eType);
    
    Sim->BCs->applyBC(retVal.rij, retVal.vijold);
  
    size_t dim(0);
   
    for (size_t iDim(1); iDim < NDIM; ++iDim)
      if (fabs(retVal.rij[dim]) < fabs(retVal.rij[iDim])) dim = iDim;

    double p1Mass = Sim->species[retVal.particle1_.getSpeciesID()]->getMass(particle1.getID()); 
    double p2Mass = Sim->species[retVal.particle2_.getSpeciesID()]->getMass(particle2.getID());
    double mu = p1Mass * p2Mass/ (p1Mass + p2Mass);
  
    Vector collvec(0,0,0);

    if (retVal.rij[dim] < 0)
      collvec[dim] = -1;
    else
      collvec[dim] = 1;

    retVal.rvdot = (retVal.rij | retVal.vijold);

    retVal.dP = collvec * (1.0 + e) * mu * (collvec | retVal.vijold);  

    //This function must edit particles so it overrides the const!
    particle1.getVelocity() -= retVal.dP / p1Mass;
    particle2.getVelocity() += retVal.dP / p2Mass;

    retVal.particle1_.setDeltaKE(0.5 * p1Mass * (particle1.getVelocity().nrm2() 
						 - retVal.particle1_.getOldVel().nrm2()));
  
    retVal.particle2_.setDeltaKE(0.5 * p2Mass * (particle2.getVelocity().nrm2() 
						 - retVal.particle2_.getOldVel().nrm2()));

    return retVal;
  }

  NEventData 
  DynNewtonian::multibdyCollision(const IDRange& range1, const IDRange& range2, const double&, const EEventType& eType) const
  {
    Vector COMVel1(0,0,0), COMVel2(0,0,0), COMPos1(0,0,0), COMPos2(0,0,0);
  
    double structmass1(0), structmass2(0);
  
    BOOST_FOREACH(const size_t& ID, range1)
      {
	updateParticle(Sim->particles[ID]);
      
	double mass = Sim->species[Sim->particles[ID]]->getMass(ID);
	structmass1 += mass;
      
	Vector pos(Sim->particles[ID].getPosition()),
	  vel(Sim->particles[ID].getVelocity());

	Sim->BCs->applyBC(pos, vel);

	COMVel1 += vel * mass;

	COMPos1 += pos * mass;
      }
  
    BOOST_FOREACH(const size_t& ID, range2)
      {
	updateParticle(Sim->particles[ID]);

	double mass = Sim->species[Sim->particles[ID]]->getMass(ID);
	structmass2 += mass;
      
	Vector pos(Sim->particles[ID].getPosition()),
	  vel(Sim->particles[ID].getVelocity());

	Sim->BCs->applyBC(pos, vel);

	COMVel2 += vel * mass;      

	COMPos2 += pos * mass;
      }
  
    COMVel1 /= structmass1;
    COMVel2 /= structmass2;
  
    COMPos1 /= structmass1;
    COMPos2 /= structmass2;
  
    Vector  rij = COMPos1 - COMPos2, vij = COMVel1 - COMVel2;
    Sim->BCs->applyBC(rij, vij);
    double rvdot = (rij | vij);

    double mu = structmass1 * structmass2 / (structmass1 + structmass2);

    static const double e = 1.0;
    Vector  dP = rij * ((1.0 + e) * mu * rvdot / rij.nrm2());

    NEventData retVal;
    BOOST_FOREACH(const size_t& ID, range1)
      {
	ParticleEventData tmpval(Sim->particles[ID],
				 *Sim->species[Sim->particles[ID]],
				 eType);

	Sim->particles[ID].getVelocity() -= dP / structmass1;
      
	tmpval.setDeltaKE(0.5 * Sim->species[tmpval.getSpeciesID()]->getMass(ID)
			  * (Sim->particles[ID].getVelocity().nrm2() 
			     - tmpval.getOldVel().nrm2()));
      
	retVal.L1partChanges.push_back(tmpval);
      }

    BOOST_FOREACH(const size_t& ID, range2)
      {
	ParticleEventData tmpval
	  (Sim->particles[ID],
	   *Sim->species[Sim->particles[ID]],
	   eType);

	Sim->particles[ID].getVelocity() += dP / structmass2;
      
	tmpval.setDeltaKE(0.5 * Sim->species[tmpval.getSpeciesID()]->getMass(ID)
			  * (Sim->particles[ID].getVelocity().nrm2() 
			     - tmpval.getOldVel().nrm2()));
  
	retVal.L1partChanges.push_back(tmpval);
      }
  
    return retVal;
  }

  NEventData 
  DynNewtonian::multibdyWellEvent(const IDRange& range1, const IDRange& range2, const double&, const double& deltaKE, EEventType& eType) const
  {
    Vector  COMVel1(0,0,0), COMVel2(0,0,0), COMPos1(0,0,0), COMPos2(0,0,0);
  
    double structmass1(0), structmass2(0);
  
    BOOST_FOREACH(const size_t& ID, range1)
      {
	updateParticle(Sim->particles[ID]);
	double mass = Sim->species[Sim->particles[ID]]->getMass(ID);

	structmass1 += mass;

	Vector pos(Sim->particles[ID].getPosition()),
	  vel(Sim->particles[ID].getVelocity());
      
	Sim->BCs->applyBC(pos, vel);

	COMVel1 += vel * mass;

	COMPos1 += pos * mass;
      }
  
    BOOST_FOREACH(const size_t& ID, range2)
      {
	updateParticle(Sim->particles[ID]);

	double mass = Sim->species[Sim->particles[ID]]->getMass(ID);
      
	structmass2 += mass;
      
	Vector pos(Sim->particles[ID].getPosition()),
	  vel(Sim->particles[ID].getVelocity());

	Sim->BCs->applyBC(pos, vel);
	COMVel2 += vel * mass;
	COMPos2 += pos * mass;
      }
  
    COMVel1 /= structmass1;
    COMVel2 /= structmass2;
  
    COMPos1 /= structmass1;
    COMPos2 /= structmass2;
  
    Vector  rij = COMPos1 - COMPos2, vij = COMVel1 - COMVel2;
    Sim->BCs->applyBC(rij, vij);
    double rvdot = (rij | vij);

    double mu = structmass1 * structmass2 / (structmass1 + structmass2);

    double R2 = rij.nrm2();
    double sqrtArg = rvdot * rvdot + 2.0 * R2 * deltaKE / mu;

    Vector  dP(0,0,0);

    if ((deltaKE < 0) && (sqrtArg < 0))
      {
	eType = BOUNCE;
	dP = rij * 2.0 * mu * rvdot / R2;
      }
    else
      {
	if (deltaKE < 0)
	  eType = WELL_KEDOWN;
	else
	  eType = WELL_KEUP;
	  
	if (rvdot < 0)
	  dP = rij 
	    * (2.0 * deltaKE / (std::sqrt(sqrtArg) - rvdot));
	else
	  dP = rij 
	    * (-2.0 * deltaKE / (rvdot + std::sqrt(sqrtArg)));
      }
  
    NEventData retVal;
    BOOST_FOREACH(const size_t& ID, range1)
      {
	ParticleEventData tmpval
	  (Sim->particles[ID],
	   *Sim->species[Sim->particles[ID]],
	   eType);

	Sim->particles[ID].getVelocity() -= dP / structmass1;
      
	tmpval.setDeltaKE(0.5 * Sim->species[tmpval.getSpeciesID()]->getMass(ID)
			  * (Sim->particles[ID].getVelocity().nrm2() 
			     - tmpval.getOldVel().nrm2()));
        
	retVal.L1partChanges.push_back(tmpval);
      }

    BOOST_FOREACH(const size_t& ID, range2)
      {
	ParticleEventData tmpval
	  (Sim->particles[ID],
	   *Sim->species[Sim->particles[ID]],
	   eType);

	Sim->particles[ID].getVelocity() += dP / structmass2;
      
	tmpval.setDeltaKE(0.5 * Sim->species[tmpval.getSpeciesID()]->getMass(ID)
			  * (Sim->particles[ID].getVelocity().nrm2() 
			     - tmpval.getOldVel().nrm2()));
      
	retVal.L1partChanges.push_back(tmpval);
      }
  
    return retVal;
  }

  PairEventData 
  DynNewtonian::SphereWellEvent(const IntEvent& event, const double& deltaKE, const double &) const
  {
    Particle& particle1 = Sim->particles[event.getParticle1ID()];
    Particle& particle2 = Sim->particles[event.getParticle2ID()];

    updateParticlePair(particle1, particle2);  

    PairEventData retVal(particle1, particle2,
			 *Sim->species[particle1],
			 *Sim->species[particle2],
			 event.getType());
    
    Sim->BCs->applyBC(retVal.rij,retVal.vijold);
  
    retVal.rvdot = (retVal.rij | retVal.vijold);
  
    double p1Mass = Sim->species[retVal.particle1_.getSpeciesID()]->getMass(particle1.getID());
    double p2Mass = Sim->species[retVal.particle2_.getSpeciesID()]->getMass(particle2.getID());
    double mu = p1Mass * p2Mass / (p1Mass + p2Mass);  
    double R2 = retVal.rij.nrm2();
    double sqrtArg = retVal.rvdot * retVal.rvdot + 2.0 * R2 * deltaKE / mu;
  
    if ((deltaKE < 0) && (sqrtArg < 0))
      {
	event.setType(BOUNCE);
	retVal.setType(BOUNCE);
	retVal.dP = retVal.rij * 2.0 * mu * retVal.rvdot / R2;
      }
    else if (deltaKE==0)
      {
	event.setType(NON_EVENT);
	retVal.setType(NON_EVENT);
	retVal.dP = Vector(0,0,0);
      }
    else
      {
	if (deltaKE < 0)
	  {
	    event.setType(WELL_KEDOWN);
	    retVal.setType(WELL_KEDOWN);
	  }
	else
	  {
	    event.setType(WELL_KEUP);
	    retVal.setType(WELL_KEUP);	  
	  }
	  
	retVal.particle1_.setDeltaU(-0.5 * deltaKE);
	retVal.particle2_.setDeltaU(-0.5 * deltaKE);	  
      
	if (retVal.rvdot < 0)
	  retVal.dP = retVal.rij 
	    * (2.0 * deltaKE / (std::sqrt(sqrtArg) - retVal.rvdot));
	else
	  retVal.dP = retVal.rij 
	    * (-2.0 * deltaKE / (retVal.rvdot + std::sqrt(sqrtArg)));
      }
  
#ifdef DYNAMO_DEBUG
    if (boost::math::isnan(retVal.dP[0]))
      M_throw() << "A nan dp has ocurred";
#endif
  
    //This function must edit particles so it overrides the const!
    particle1.getVelocity() -= retVal.dP / p1Mass;
    particle2.getVelocity() += retVal.dP / p2Mass;
  
    retVal.particle1_.setDeltaKE(0.5 * p1Mass
				 * (particle1.getVelocity().nrm2() 
				    - retVal.particle1_.getOldVel().nrm2()));
  
    retVal.particle2_.setDeltaKE(0.5 * p2Mass
				 * (particle2.getVelocity().nrm2() 
				    - retVal.particle2_.getOldVel().nrm2()));

    return retVal;
  }

  void 
  DynNewtonian::outputXML(magnet::xml::XmlStream& XML) const
  {
    XML << magnet::xml::attr("Type") 
	<< "Newtonian";
  }

  double 
  DynNewtonian::getPBCSentinelTime(const Particle& part, const double& lMax) const
  {
#ifdef DYNAMO_DEBUG
    if (!isUpToDate(part))
      M_throw() << "Particle is not up to date";
#endif

    Vector pos(part.getPosition()), vel(part.getVelocity());

    Sim->BCs->applyBC(pos, vel);

    double retval = HUGE_VAL;

    for (size_t i(0); i < NDIM; ++i)
      if (vel[i] != 0)
	{
	  double tmp = (0.5 * (0.5 * Sim->primaryCellSize[i] - lMax)) / fabs(vel[i]);
	  
	  if (tmp < retval)
	    retval = tmp;
	}

    return retval;
  }

  std::pair<bool,double>
  DynNewtonian::getPointPlateCollision(const Particle& part, const Vector& nrw0, const Vector& nhat, const double& Delta, const double& Omega, const double& Sigma, const double& t, bool lastpart) const
  {
#ifdef DYNAMO_DEBUG
    if (!isUpToDate(part))
      M_throw() << "Particle1 " << part.getID() << " is not up to date";
#endif
  
    Vector pos(part.getPosition() - nrw0), vel(part.getVelocity());
    Sim->BCs->applyBC(pos, vel);

    double t_high;
    double surfaceOffset = pos | nhat;
    double surfaceVel = vel | nhat;
  
    if (surfaceVel > 0)
      t_high = (Sigma + Delta - surfaceOffset) / surfaceVel;
    else
      t_high = -(Sigma + Delta + surfaceOffset) / surfaceVel;
  

    SFOscillatingPlate fL(vel, nhat, pos, t, Delta, Omega, Sigma);
  
#ifdef DYNAMO_DEBUG
    if (Sigma < 0) M_throw() << "Assuming a positive Sigma here";
#endif

    //A particle has penetrated the plate, probably due to some small numerical error
    //We can just adjust the seperation vector till the particle is on the surface of the plate
    if (fL.F_zeroDeriv() > 0)
      {
#ifdef DYNAMO_DEBUG
	derr << "Particle is penetrating the \"upper\" plate"
	     << "\nTo avoid rediscovering the root we're adjusting the relative position vector to just touching."
	     << "\nThis is fine if it is a rare event." << std::endl;
#endif
	fL.fixFZeroSign(false);

#ifdef DYNAMO_DEBUG
	//This is just incase the oscillating plate shape function is broken
	if (fL.F_zeroDeriv() > 0)
	  M_throw() << "Failed to adjust the plate position";
#endif
      }
      
    double t_low1 = 0, t_low2 = 0;
    if (lastpart)
      {
	if (-fL.F_zeroDeriv() < fL.F_zeroDerivFlip())
	  //Shift the lower bound up so we don't find the same root again
	  t_low1 = fabs(2.0 * fL.F_firstDeriv())
	    / fL.F_secondDeriv_max();
	else
	  t_low2 = fabs(2.0 * fL.F_firstDeriv())
	    / fL.F_secondDeriv_max();
      }


    //Must be careful with collisions at the end of the interval
    t_high *= 1.01;
  
    std::pair<bool,double> root1 
      = magnet::math::frenkelRootSearch(fL, t_low1, t_high, 1e-12 * Sigma);

    fL.flipSigma();

    if (fL.F_zeroDeriv() < 0)
      {
#ifdef DYNAMO_DEBUG
	derr << "Particle is penetrating the \"lower\" plate"
	     << "\nTo avoid rediscovering the root we're adjusting the relative position vector to just touching."
	     << "\nThis is fine if it is a rare event." << std::endl;
#endif
	fL.fixFZeroSign(true);

#ifdef DYNAMO_DEBUG
	//This is just incase the oscillating plate shape function is broken
	if (fL.F_zeroDeriv() < 0)
	  M_throw() << "Failed to adjust the plate position";
#endif
      }

    std::pair<bool,double> root2 
      = magnet::math::frenkelRootSearch(fL, t_low2, t_high, 1e-12 * Sigma);

    //Check if the particle is penetrating a wall
    //Or if no roots are found at all
    if ((fabs(surfaceOffset - (nhat | fL.wallPosition())) > Sigma)
	|| ((root1.second == HUGE_VAL) && (root2.second == HUGE_VAL))
	|| ((t_low1 > t_high) && (t_low2 > t_high)))
      {
	//This can be a problem
#ifdef DYNAMO_DEBUG      
	derr << "Particle " << part.getID() 
	     << " may be outside/heading out of the plates"
	     << "\nerror = "
	     << (fabs(surfaceOffset - (nhat | fL.wallPosition())) - Sigma) 
	  / Sim->units.unitLength()
	     << "\n Root1 = " 
	     << root1.second / Sim->units.unitTime()
	     << "\n Root2 = " 
	     << root2.second / Sim->units.unitTime() << std::endl;
#endif
      
	//If the particle is going out of bounds, collide now
	if (fL.test_root())
	  {
#ifdef DYNAMO_DEBUG
	    {
	      SFOscillatingPlate ftmp(fL);
	      SFOscillatingPlate ftmp2(fL);
	      ftmp.flipSigma();
	    
	      double fl01(ftmp.F_zeroDeriv());
	      ftmp.stream(t_low1);
	      double flt_low1(ftmp.F_zeroDeriv());
	      ftmp.stream(t_high - t_low1);
	      double flt_high1(ftmp.F_zeroDeriv());
	    
	      double fl02(ftmp2.F_zeroDeriv());
	      ftmp2.stream(t_low2);
	      double flt_low2(ftmp2.F_zeroDeriv());
	      ftmp2.stream(t_high - t_low2);
	      double flt_high2(ftmp2.F_zeroDeriv());
	    
	      derr << "****Forcing collision"
		   << "\nsystemTime = " << Sim->systemTime
		   << "\nlNColl = " << Sim->eventCount
		   << "\nlast part = " << (lastpart ? (std::string("True")) : (std::string("False")))
		   << "\nVel = " << part.getVelocity()[0]
		   << "\nPos = " << part.getPosition()[0]
		   << "\nVwall[0] = " << fL.wallVelocity()[0]
		   << "\nRwall[0] = " << fL.wallPosition()[0]
		   << "\nRwall[0]+Sigma = " << fL.wallPosition()[0] + Sigma
		   << "\nRwall[0]-Sigma = " << fL.wallPosition()[0] - Sigma
		   << "\nSigma + Del = " << Sigma+Delta
		   << "\nGood root = " << fL.test_root()
		   << "\nt_low1 = " << t_low1
		   << "\nt_low2 = " << t_low2
		   << "\nt_high = " << t_high
		   << "\nroot1 = " << root1.second
		   << "\nroot2 = " << root2.second
		   << "\nf1(0) = " << fl01
		   << "\nf1(t_low1) = " << flt_low1
		   << "\nf1(t_high) = " << flt_high1
		   << "\nf2(0)_1 = " << fl02
		   << "\nf2(t_low2) = " << flt_low2
		   << "\nf2(t_high) = " << flt_high2
		   << "\nf'(0) =" << fL.F_firstDeriv()
		   << "\nf''(Max) =" << fL.F_secondDeriv_max()
		   << "\nf(x)=" << (pos | nhat)
		   << "+" << (part.getVelocity() | nhat)
		   << " * x - "
		   << Delta 
		   << " * cos(("
		   << t + Sim->systemTime << "+ x) * "
		   << Omega << ") - "
		   << Sigma
		   << " << std::endl; set xrange[0:" << t_high << "]; plot f(x)";
	      ;
	    }
#endif
	    return std::pair<bool, double>(true, 0);
	  }
	else
	  {
	    //The particle and plate are approaching but might not be
	    //before the overlap is fixed, schedule another test later
	    //on
	    double currRoot = HUGE_VAL;

	    if (root1.first)
	      currRoot = root1.second;

	    if (root2.first && (currRoot > root2.second))
	      currRoot = root2.second;
	    //
	    double tmpt = fabs(surfaceVel - fL.velnHatWall());
	    //This next line sets what the recoil velocity should be
	    //We choose the velocity that gives elastic collisions!
	    tmpt += fL.maxWallVel() * 0.002;
	    tmpt /= fL.F_secondDeriv_max();
	    if (tmpt < currRoot)
	      {
#ifdef DYNAMO_DEBUG
		dout << "Making a fake collision at " << tmpt << "for particle " << part.getID() << std::endl;
#endif

		return std::pair<bool,double>(true, tmpt);
	      }
#ifdef DYNAMO_DEBUG
	    else
	      dout << "The current root is lower than the fake one" << std::endl;	    
#endif
	  }
      }
  
    return (root1.second < root2.second) ? root1 : root2;
  }

  ParticleEventData 
  DynNewtonian::runOscilatingPlate(Particle& part, const Vector& rw0, const Vector& nhat, double& delta, const double& omega0, const double& sigma, const double& mass, const double& e, double& t, bool strongPlate) const
  {
    std::cout.flush();
    updateParticle(part);

    ParticleEventData retVal(part, *Sim->species[part], WALL);

    SFOscillatingPlate fL(part.getVelocity(), nhat, part.getPosition(),
			     t + Sim->systemTime, delta, omega0, sigma);

    //Should force the particle to the plate surface

    Vector pos(part.getPosition() - fL.wallPosition()), vel(part.getVelocity());

    Sim->BCs->applyBC(pos, vel);
  
    double pmass = Sim->species[retVal.getSpeciesID()]->getMass(part.getID());
    double mu = (pmass * mass) / (mass + pmass);

    Vector vwall(fL.wallVelocity());

    //  derr << "Running event for part " << part.getID() <<
    //	   "\nsystemTime = " << Sim->systemTime << "\nlNColl = " <<
    //	   Sim->lNColl << "\nVel = " << part.getVelocity()[0] <<
    //	   "\nPos = " << part.getPosition()[0] << "\nVwall[0] = " <<
    //	   fL.wallVelocity()[0] << "\nRwall[0] = " <<
    //	   fL.wallPosition()[0] << "\nRwall[0]+sigma = " <<
    //	   fL.wallPosition()[0] + sigma << "\nRwall[0]-sigma = " <<
    //	   fL.wallPosition()[0] - sigma << "\nsigma + Del = " <<
    //	   sigma+delta << "\nf(0)* = " << fL.F_zeroDeriv() << "\nf'(0)
    //	   =" << fL.F_firstDeriv() << "\nf''(Max) =" <<
    //	   fL.F_secondDeriv_max(0) << "\nf(x)=" << pos[0] << "+" <<
    //	   part.getVelocity()[0] << " * x - " << delta << " * cos(("
    //	   << t << "+ x) * " << omega0 << ") - " << sigma << std::endl;

  
    //Check the root is valid
    if (!fL.test_root())
      {
	double f0 = fL.F_zeroDeriv(), f1 = fL.F_firstDeriv(),
	  f2 = fL.F_secondDeriv_max();
	fL.flipSigma();
      
	derr <<"Particle " << part.getID()
	     << ", is pulling on the oscillating plate!"
	     << "\nRunning event for part " << part.getID()
	     << "\nsystemTime = " << Sim->systemTime
	     << "\nlNColl = " << Sim->eventCount
	     << "\nVel = " << part.getVelocity()[0]
	     << "\nPos = " << part.getPosition()[0]
	     << "\nVwall[0] = " << fL.wallVelocity()[0]
	     << "\nRwall[0] = " << fL.wallPosition()[0]
	     << "\nRwall[0]+sigma = " << fL.wallPosition()[0] + sigma
	     << "\nRwall[0]-sigma = " << fL.wallPosition()[0] - sigma
	     << "\nGood root " << fL.test_root()
	     << "\nsigma + Del = " << sigma+delta
	     << "\nf1(0)* = " << fL.F_zeroDeriv()
	     << "\nf1'(0) =" << fL.F_firstDeriv()
	     << "\nf1''(Max) =" << fL.F_secondDeriv_max()
	     << "\nf2(0)* = " << f0
	     << "\nf2'(0) =" << f1
	     << "\nf2''(Max) =" << f2
	     << "\nf(x)=" << (pos | nhat)
	     << "+" << (part.getVelocity() | nhat)
	     << " * x - "
	     << delta
	     << " * cos(("
	     << t + Sim->systemTime << "+ x) * "
	     << omega0 << ") - "
	     << sigma << std::endl;
      
	return retVal;
      }

    //static size_t elascount(0);

    double inelas = e;

    double rvdot = ((vel - vwall) | nhat);
    if (fabs(rvdot / fL.maxWallVel()) < 0.002)
      {
	/*
	  derr <<"<<!!!>>Particle " << part.getID() 
	  << " gone elastic!\nratio is " << fabs(rvdot / vwall.nrm())
	  << "\nCount is " << ++elascount << std::endl;
	*/
	inelas = 1.0;
	if (fabs(rvdot / fL.maxWallVel()) < 0.001)
	  {
	    if (rvdot < 0)
	      rvdot = -fL.maxWallVel() * 0.01;
	    else
	      rvdot = fL.maxWallVel() * 0.01;
	  }
      }

    Vector delP =  nhat * mu * (1.0 + inelas) * rvdot;

    part.getVelocity() -=  delP / pmass;

    retVal.setDeltaKE(0.5 * pmass
		      * (part.getVelocity().nrm2() 
			 - retVal.getOldVel().nrm2()));
  
    //Don't progress if you want to not change the plate data
    if (strongPlate) return retVal;

    double numerator = -nhat | ((delP / mass) + vwall);

    double reducedt = Sim->systemTime 
      - 2.0 * M_PI * int(Sim->systemTime * omega0 / (2.0*M_PI)) / omega0;
  
    double denominator = omega0 * delta * std::cos(omega0 * (reducedt + t));
  

    double newt = std::atan2(numerator, denominator)/ omega0 
      - Sim->systemTime;
  
    delta *= std::cos(omega0 * (Sim->systemTime + t)) 
      / std::cos(omega0 * (Sim->systemTime + newt));
  
    t = newt;

    t -= 2.0 * M_PI * int(t * omega0 / (2.0*M_PI)) / omega0;

    return retVal; 
  }

  double 
  DynNewtonian::getCylinderWallCollision(const Particle& part, const Vector& wallLoc, const Vector& wallNorm, const double& radius) const
  {
    Vector  rij = part.getPosition() - wallLoc,
      vel = part.getVelocity();

    Sim->BCs->applyBC(rij, vel);

    rij -= Vector((rij | wallNorm) * wallNorm);

    vel -= Vector((vel | wallNorm) * wallNorm);

    double B = (vel | rij),
      A = vel.nrm2(),
      C = rij.nrm2() - radius * radius;

    double t = (std::sqrt(B*B - A*C) - B) / A;

    if (boost::math::isnan(t))
      return HUGE_VAL;
    else
      return t;
  }

  ParticleEventData 
  DynNewtonian::runCylinderWallCollision(Particle& part, const Vector& origin, const Vector& vNorm, const double& e) const
  {
    updateParticle(part);

    ParticleEventData retVal(part, *Sim->species[part], WALL);
  
    Vector rij =  origin - part.getPosition();

    Sim->BCs->applyBC(rij);

    rij -= Vector((rij | vNorm) * vNorm);

    rij /= rij.nrm();

    part.getVelocity() -= (1+e) * (rij | part.getVelocity()) * rij;
  
    retVal.setDeltaKE(0.5 * Sim->species[retVal.getSpeciesID()]->getMass(part.getID())
		      * (part.getVelocity().nrm2() 
			 - retVal.getOldVel().nrm2()));
  
    return retVal; 
  }

  std::pair<bool, double> 
  DynNewtonian::getLineLineCollision(const double length, const Particle& p1, const Particle& p2, double t_high) const
  {  
#ifdef DYNAMO_DEBUG
    if (!hasOrientationData())
      M_throw() << "Cannot use this function without orientational data";

    if (!isUpToDate(p1))
      M_throw() << "Particle1 " << p1.getID() << " is not up to date";

    if (!isUpToDate(p2))
      M_throw() << "Particle2 " << p2.getID() << " is not up to date";
#endif

    Vector r12 = p1.getPosition() - p2.getPosition();
    Vector v12 = p1.getVelocity() - p2.getVelocity();
    Sim->BCs->applyBC(r12, v12);

    double t_low = 0.0;
  
    SFLines fL(r12, v12,
		  orientationData[p1.getID()].angularVelocity,
		  orientationData[p2.getID()].angularVelocity,
		  orientationData[p1.getID()].orientation,
		  orientationData[p2.getID()].orientation,
		  length);
  
    if (((p1.getID() == lastCollParticle1 && p2.getID() == lastCollParticle2)
	 || (p1.getID() == lastCollParticle2 && p2.getID() == lastCollParticle1))
	&& Sim->systemTime == lastAbsoluteClock)
      //Shift the lower bound up so we don't find the same root again
      t_low += fabs(2.0 * fL.F_firstDeriv())
	/ fL.F_secondDeriv_max();
  
    //Find window delimited by discs
    std::pair<double,double> dtw = fL.discIntersectionWindow();
  
    if(dtw.first > t_low)
      t_low = dtw.first;
  
    if(dtw.second < t_high)
      t_high = dtw.second;
  
    return magnet::math::frenkelRootSearch(fL, t_low, t_high, length * 1e-10);
  }

  std::pair<bool, double> 
  DynNewtonian::getOffcentreSpheresCollision(const double offset1, const double diameter1, const double offset2, const double diameter2, const Particle& p1, const Particle& p2, double t_max) const
  {  
#ifdef DYNAMO_DEBUG
    if (!hasOrientationData())
      M_throw() << "Cannot use this function without orientational data";

    if (!isUpToDate(p1))
      M_throw() << "Particle1 " << p1.getID() << " is not up to date";

    if (!isUpToDate(p2))
      M_throw() << "Particle2 " << p2.getID() << " is not up to date";
#endif

    Vector r12 = p1.getPosition() - p2.getPosition();
    Vector v12 = p1.getVelocity() - p2.getVelocity();
    Sim->BCs->applyBC(r12, v12);
  
    SFOffcentre_Spheres fL(r12, v12,
			   orientationData[p1.getID()].angularVelocity,
			   orientationData[p2.getID()].angularVelocity,
			   orientationData[p1.getID()].orientation,
			   orientationData[p2.getID()].orientation,
			   offset1);
    
    double t_low = 0;
    if (((p1.getID() == lastCollParticle1 && p2.getID() == lastCollParticle2)
	 || (p1.getID() == lastCollParticle2 && p2.getID() == lastCollParticle1))
	&& Sim->systemTime == lastAbsoluteClock)
      //Shift the lower bound up so we don't find the same root again
      t_low += fabs(2.0 * fL.F_firstDeriv())
	/ fL.F_secondDeriv_max();
    
    return magnet::math::frenkelRootSearch(fL, t_low, t_max, std::min(diameter1, diameter2) * 1e-10);
  }


  PairEventData 
  DynNewtonian::runLineLineCollision(const IntEvent& eevent, const double& elasticity, const double& length) const
  {
#ifdef DYNAMO_DEBUG
    if (!hasOrientationData())
      M_throw() << "Cannot use this function without orientational data";
#endif

    Particle& particle1 = Sim->particles[eevent.getParticle1ID()];
    Particle& particle2 = Sim->particles[eevent.getParticle2ID()];

    updateParticlePair(particle1, particle2);  

    PairEventData retVal(particle1, particle2,
			 *Sim->species[particle1],
			 *Sim->species[particle2],
			 CORE);
  
    Sim->BCs->applyBC(retVal.rij, retVal.vijold);

    retVal.rvdot = (retVal.rij | retVal.vijold);

    double KE1before = getParticleKineticEnergy(particle1);
    double KE2before = getParticleKineticEnergy(particle2);

    SFLines fL(retVal.rij, retVal.vijold,
		  orientationData[particle1.getID()].angularVelocity,
		  orientationData[particle2.getID()].angularVelocity,
		  orientationData[particle1.getID()].orientation,
		  orientationData[particle2.getID()].orientation,
		  length);

    Vector uPerp = fL.getu1() ^ fL.getu2();

    uPerp /= uPerp.nrm();

    std::pair<double, double> cp = fL.getCollisionPoints();

    // \Delta {\bf v}_{imp}
    Vector  vr = retVal.vijold
      + (cp.first * fL.getw1() ^ fL.getu1()) 
      - (cp.second * fL.getw2() ^ fL.getu2());
  
    double mass = Sim->species[retVal.particle1_.getSpeciesID()]->getMass(particle1.getID());
    double inertia = Sim->species[retVal.particle1_.getSpeciesID()]->getScalarMomentOfInertia(particle1.getID());

    retVal.dP = uPerp
      * (((vr | uPerp) * (1.0 + elasticity))
	 / ((2.0 / mass) + ((cp.first * cp.first + cp.second * cp.second) / inertia)));
  
    particle1.getVelocity() -= retVal.dP / mass;
    particle2.getVelocity() += retVal.dP / mass;

    orientationData[particle1.getID()].angularVelocity 
      -= (cp.first / inertia) * (fL.getu1() ^ retVal.dP);

    orientationData[particle2.getID()].angularVelocity 
      += (cp.second / inertia) * (fL.getu2() ^ retVal.dP);

    retVal.particle1_.setDeltaKE(getParticleKineticEnergy(particle1) - KE1before);
    retVal.particle2_.setDeltaKE(getParticleKineticEnergy(particle2) - KE2before);

    lastCollParticle1 = particle1.getID();
    lastCollParticle2 = particle2.getID();
    lastAbsoluteClock = Sim->systemTime;

    return retVal;
  }

  double 
  DynNewtonian::sphereOverlap(const Particle& p1, const Particle& p2, const double& d) const
  {
    Vector r12 = p1.getPosition() - p2.getPosition();
    Sim->BCs->applyBC(r12);

    return std::sqrt(std::max(d * d - (r12 | r12), 0.0));
  }

  PairEventData 
  DynNewtonian::RoughSpheresColl(const IntEvent& event, const double& e, const double& et, const double& d2, const EEventType& eType) const
  {
#ifdef DYNAMO_DEBUG
    if (!hasOrientationData())
      M_throw() << "Cannot use this function without orientational data";
#endif

    Particle& particle1 = Sim->particles[event.getParticle1ID()];
    Particle& particle2 = Sim->particles[event.getParticle2ID()];

    updateParticlePair(particle1, particle2);  

    PairEventData retVal(particle1, particle2,
			 *Sim->species[particle1],
			 *Sim->species[particle2],
			 eType);
    
    Sim->BCs->applyBC(retVal.rij, retVal.vijold);
  
    double p1Mass = Sim->species[retVal.particle1_.getSpeciesID()]->getMass(particle1.getID()); 
    double p2Mass = Sim->species[retVal.particle2_.getSpeciesID()]->getMass(particle2.getID());
    double mu = p1Mass * p2Mass/(p1Mass+p2Mass);
  
    retVal.rvdot = (retVal.rij | retVal.vijold);

    //The normal impulse
    retVal.dP = retVal.rij * ((1.0 + e) * mu * retVal.rvdot / retVal.rij.nrm2());

    Vector eijn = retVal.rij / retVal.rij.nrm();

    //Now the tangential impulse
    Vector gij = retVal.vijold - std::sqrt(d2) * 0.5 
      * ((orientationData[particle1.getID()].angularVelocity
	  + orientationData[particle2.getID()].angularVelocity) ^ eijn);
  
    Vector gijt = (eijn ^ gij) ^ eijn;

    double Jbar = Sim->species[retVal.particle1_.getSpeciesID()]->getScalarMomentOfInertia(particle1.getID()) 
      / (p1Mass * d2 * 0.25);
  
    retVal.dP += (Jbar * (1-et) / (2*(Jbar + 1))) * gijt;

    double KE1before = getParticleKineticEnergy(particle1);
    double KE2before = getParticleKineticEnergy(particle2);

    //This function must edit particles so it overrides the const!
    particle1.getVelocity() -= retVal.dP / p1Mass;
    particle2.getVelocity() += retVal.dP / p2Mass;

    Vector angularVchange = (1-et) / (std::sqrt(d2) * (Jbar+1)) * (eijn ^ gijt);
 
    orientationData[particle1.getID()].angularVelocity
      += angularVchange;
    orientationData[particle2.getID()].angularVelocity 
      += angularVchange;

    retVal.particle1_.setDeltaKE(getParticleKineticEnergy(particle1) - KE1before);
    retVal.particle2_.setDeltaKE(getParticleKineticEnergy(particle2) - KE2before);

    return retVal;
  }

  ParticleEventData 
  DynNewtonian::runRoughWallCollision(Particle& part, const Vector & vNorm, const double& e, const double& et, const double& r) const
  {
#ifdef DYNAMO_DEBUG
    if (!hasOrientationData())
      M_throw() << "Cannot use this function without orientational data";
#endif

    updateParticle(part);

    ParticleEventData retVal(part, *Sim->species[part], WALL);

    double KE1before = getParticleKineticEnergy(part);

    double p1Mass = Sim->species[retVal.getSpeciesID()]->getMass(part.getID()); 

    double Jbar = Sim->species[retVal.getSpeciesID()]->getScalarMomentOfInertia(part.getID())
      / (p1Mass * r * r);

    Vector gij = part.getVelocity() - r
      * (orientationData[part.getID()].angularVelocity ^ vNorm);
  
    Vector gijt = (vNorm ^ gij) ^ vNorm;
  
    part.getVelocity() -= (1+e) * (vNorm | part.getVelocity()) * vNorm
      + (Jbar * (1-et) / (Jbar + 1)) * gijt;

    Vector angularVchange = (1-et) / (r * (Jbar+1)) * (vNorm ^ gijt);
 
    orientationData[part.getID()].angularVelocity
      += angularVchange;

    retVal.setDeltaKE(getParticleKineticEnergy(part) - KE1before);
    return retVal; 
  }
}
