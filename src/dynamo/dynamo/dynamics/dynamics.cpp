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

#include <dynamo/dynamics/include.hpp>
#include <dynamo/species/inertia.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/NparticleEventData.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/BC/LEBC.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <cstring>

namespace dynamo {
  magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream& XML, const Dynamics& g)
  {
    g.outputXML(XML);
    return XML;
  }

  shared_ptr<Dynamics>
  Dynamics::getClass(const magnet::xml::Node& XML, dynamo::Simulation* tmp)
  {
    if (!XML.getAttribute("Type").getValue().compare("Newtonian"))
      return shared_ptr<Dynamics>(new DynNewtonian(tmp));
    if (!XML.getAttribute("Type").getValue().compare("NewtonianGravity"))
      return shared_ptr<Dynamics>(new DynGravity(tmp, XML));
    else if (!XML.getAttribute("Type").getValue().compare("NewtonianMC"))
      return shared_ptr<Dynamics>(new DynNewtonianMC(tmp, XML));
    else if (!XML.getAttribute("Type").getValue().compare("NewtonianMCCMap"))
      return shared_ptr<Dynamics>(new DynNewtonianMCCMap(tmp, XML));
    else
      M_throw() << XML.getAttribute("Type").getValue()
		<< ", Unknown type of Dynamics encountered";
  }


  void 
  Dynamics::initialise()
  {
    streamFreq = 10 * Sim->N();
    
    if (hasOrientationData())
      {
	double sumEnergy(0.0);
	for (const Particle& part : Sim->particles)
	  sumEnergy += Sim->species(part)->getScalarMomentOfInertia(part.getID())
	  * orientationData[part.getID()].angularVelocity.nrm2();
      
	//Check if any of the species are overridden
	bool hasInertia(false);
	for (const shared_ptr<Species>& spec : Sim->species)
	  if (std::dynamic_pointer_cast<SpInertia>(spec))
	    hasInertia = true;

	if (!hasInertia)
	  M_throw() << "No species have inertia, yet the particles have orientational degrees of freedom set!";

	sumEnergy *= 0.5 / Sim->units.unitEnergy();
  
	dout << "System Rotational Energy " << sumEnergy
	     << "\nRotational kT " << sumEnergy / Sim->N() << std::endl;
      }
  }

  NEventData
  Dynamics::enforceParabola(Particle&) const
  {
    M_throw() << "This is not needed for this type of Dynamics";
  }

  PairEventData 
  Dynamics::runLineLineCollision(Event&, const double&, const double&) const
  { M_throw() << "Not implemented for this Dynamics."; }
  
  std::pair<bool, double>
  Dynamics::getLineLineCollision(const double, const Particle&, const Particle&, double) const
  { M_throw() << "Not implemented for this Dynamics."; }

  std::pair<bool, double> 
  Dynamics::getOffcentreSpheresCollision(const double offset1, const double diameter1, 
					 const double offset2, const double diameter2,
					 const Particle& p1, const Particle& p2,
					 double t_max, double maxdist) const
  { M_throw() << "Not implemented for this Dynamics."; }

  double 
  Dynamics::getPBCSentinelTime(const Particle&, const double&) const
  { M_throw() << "Not implemented for this Dynamics."; }

  void 
  Dynamics::loadParticleXMLData(const magnet::xml::Node& XML)
  {
    dout << "Loading Particle Data" << std::endl;

    bool outofsequence = false;  
  
    for (magnet::xml::Node node = XML.getNode("ParticleData").findNode("Pt"); 
	 node.valid(); ++node)
      {
	if (!node.hasAttribute("ID")
	    || node.getAttribute("ID").as<size_t>() != Sim->particles.size())
	  outofsequence = true;
      
	Particle part(node, Sim->particles.size());
	part.getVelocity() *= Sim->units.unitVelocity();
	part.getPosition() *= Sim->units.unitLength();
	Sim->particles.push_back(part);
      }

    if (outofsequence)
      dout << "Particle ID's out of sequence!\n"
	   << "This can result in incorrect capture map loads etc.\n"
	   << "Erase any capture maps in the configuration file so they are regenerated." << std::endl;

    dout << "Particle count " << Sim->N() << std::endl;

    if (XML.getNode("ParticleData").hasAttribute("OrientationData"))
      {
	orientationData.resize(Sim->N());
	size_t i(0);
	for (magnet::xml::Node node = XML.getNode("ParticleData").findNode("Pt"); node.valid(); ++node, ++i)
	  {
	    orientationData[i].orientation << node.getNode("U");
	    orientationData[i].angularVelocity << node.getNode("O");
      
	    //Makes the vector a unit vector
	    orientationData[i].orientation.normalise();
	    if (orientationData[i].orientation.nrm() == 0)
	      M_throw() << "Particle " << i << " has an invalid zero orientation quaternion";
	  }
      }
  }

  void 
  Dynamics::outputParticleXMLData(magnet::xml::XmlStream& XML, bool applyBC) const
  {
    XML << magnet::xml::tag("ParticleData");
  
    if (hasOrientationData())
      XML << magnet::xml::attr("OrientationData") << "Y";

    for (size_t i = 0; i < Sim->N(); ++i)
      {
	Particle tmp(Sim->particles[i]);
	if (applyBC) 
	  Sim->BCs->applyBC(tmp.getPosition(), tmp.getVelocity());
      
	tmp.getVelocity() *= (1.0 / Sim->units.unitVelocity());
	tmp.getPosition() *= (1.0 / Sim->units.unitLength());
      
	XML << magnet::xml::tag("Pt");
	Sim->_properties.outputParticleXMLData(XML, i);
	XML << tmp;

	if (hasOrientationData())
	  XML << magnet::xml::tag("O")
	      << orientationData[i].angularVelocity
	      << magnet::xml::endtag("O")
	      << magnet::xml::tag("U")
	      << orientationData[i].orientation
	      << magnet::xml::endtag("U") ;

	XML << magnet::xml::endtag("Pt");
      }
  
    XML << magnet::xml::endtag("ParticleData");
  }

  double 
  Dynamics::getParticleKineticEnergy(const Particle& part) const
  {
    const double mass = Sim->species(part)->getMass(part.getID());

    double energy(0);
    if (!std::isinf(mass))
      {
	if (std::dynamic_pointer_cast<BCLeesEdwards>(Sim->BCs))
	  energy += static_cast<const BCLeesEdwards&>(*Sim->BCs).getPeculiarVelocity(part).nrm2() * mass;
	else
	  energy += part.getVelocity().nrm2() * mass;
      }

    if (hasOrientationData())
      {
	const double I = Sim->species(part)->getScalarMomentOfInertia(part.getID());
	if (!std::isinf(I))
	  energy += I * orientationData[part.getID()].angularVelocity.nrm2();
      }

    return 0.5 * energy;
  }

  double 
  Dynamics::getSystemKineticEnergy() const
  {
    double sumEnergy(0);

    for (const Particle& part : Sim->particles)
      sumEnergy += getParticleKineticEnergy(part);

    return sumEnergy;
  }

  void
  Dynamics::rescaleSystemKineticEnergy(const double& scale)
  {
    double scalefactor(sqrt(scale));

    if (std::dynamic_pointer_cast<BCLeesEdwards>(Sim->BCs))
      {
	const BCLeesEdwards& bc = static_cast<const BCLeesEdwards&>(*Sim->BCs);
	for (Particle& part : Sim->particles)
	  {
	    const double mass = Sim->species(part)->getMass(part.getID());
	    if (!std::isinf(mass))
	      part.getVelocity() = Vector(bc.getPeculiarVelocity(part) * scalefactor + bc.getStreamVelocity(part));
	  }
      }
    else
      for (Particle& part : Sim->particles)
	{
	  const double mass = Sim->species(part)->getMass(part.getID());
	  if (!std::isinf(mass))
	    part.getVelocity() *= scalefactor;
	}

    if (hasOrientationData())
      for (Particle& part : Sim->particles)
	{
	  const double I = Sim->species(part)->getScalarMomentOfInertia(part.getID());
	  if (!std::isinf(I))
	    orientationData[part.getID()].angularVelocity *= scalefactor;
	}
  }

  PairEventData 
  Dynamics::parallelCubeColl(Event& event, const double& e, const double& d, const EEventType& eType) const
  { M_throw() << "Not Implemented"; }

  std::pair<bool,double>
  Dynamics::getPointPlateCollision(const Particle& np1, const Vector& nrw0,
				      const Vector& nhat, const double& Delta,
				      const double& Omega, const double& Sigma,
				      const double& t, bool) const
  {
    M_throw() << "Not Implemented";
  }

  ParticleEventData 
  Dynamics::runOscilatingPlate
  (Particle& part, const Vector& rw0, const Vector& nhat, double& delta, 
   const double& omega0, const double& sigma, const double& mass, const double& e,
   double& t, bool strongPlate) const
  {
    M_throw() << "Not Implemented";
  }


  double 
  Dynamics::getCylinderWallCollision(const Particle& part, 
					const Vector& origin, 
					const Vector& norm,
					const double&
					) const
  {
    M_throw() << "Not Implemented";
  }

  ParticleEventData 
  Dynamics::runCylinderWallCollision(Particle&, 
					const Vector &,
					const Vector &,
					const double&
					) const
  {
    M_throw() << "Not Implemented";
  }

  PairEventData 
  Dynamics::RoughSpheresColl(Event& event, const double& e, const double& et, const double& d1, const double& d2, const EEventType& eType) const
  {
    M_throw() << "Not Implemented, you need rotational dynamics";
  }

  ParticleEventData 
  Dynamics::runRoughWallCollision(Particle& part, const Vector & vNorm, const double& e, const double& et, const double& r) const
  {
    M_throw() << "Not Implemented, you need rotational dynamics";
  }

  void 
  Dynamics::initOrientations(double kbT)
  {
    orientationData.resize(Sim->particles.size());
  
    //std::sqrt(10.0/ (diameter * diameter))
  
    dout << "Initialising the line orientations" << std::endl;

    std::normal_distribution<> norm_dist;
    
    for (size_t i = 0; i < Sim->particles.size(); ++i)
      {
	//Assign the new velocities
	orientationData[i].orientation = Quaternion::identity();
      
	Vector angVelCrossing;
	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  angVelCrossing[iDim] = norm_dist(Sim->ranGenerator);
	
	//Ensure the initial angular velocity is perpendicular to the
	//director
	double I = Sim->species(Sim->particles[i])->getScalarMomentOfInertia(i);
	
	if (std::isinf(I))
	  orientationData[i].angularVelocity = Vector{0,0,0};
	else
	  {
	    orientationData[i].angularVelocity = Quaternion::initialDirector() ^ angVelCrossing;
	    orientationData[i].angularVelocity *= 0.5 * std::sqrt(kbT/I) * norm_dist(Sim->ranGenerator) / orientationData[i].angularVelocity.nrm();
	  }
      }
  }

  std::pair<double, Dynamics::TriangleIntersectingPart> 
  Dynamics::getSphereTriangleEvent(const Particle& part, 
				      const Vector & A, 
				      const Vector & B, 
				      const Vector & C,
				      const double d
				      ) const
  {
    M_throw() << "Not implemented";
  }

  std::pair<Vector, Vector> 
  Dynamics::getCOMPosVel(const IDRange& particles) const
  {
    if (particles.empty())
      M_throw() << "Cannot calculate the COM position and velocity from an empty IDRange";
    
    Vector pos = Vector{0,0,0}, 
      vel = Vector{0,0,0};

    Vector pos0 = Sim->particles[*(particles.begin())].getPosition(), 
      vel0 = Sim->particles[*(particles.begin())].getVelocity();

    double totMass = 0;

    for (size_t ID : particles)
      {
	const Particle& part = Sim->particles[ID];
	double mass = Sim->species(part)->getMass(ID);

	//Take everything relative to the first particle's position to
	//minimise issues with PBC wrapping the particles.
	Vector r12 = part.getPosition() - pos0;
	Vector v12 = part.getVelocity() - vel0;
	Sim->BCs->applyBC(r12, v12);

	pos += mass * r12;
	vel += mass * v12;

	totMass += mass;
      }

    pos /= totMass;
    vel /= totMass;
    
    return std::make_pair(pos + pos0, vel + vel0);
  }

}
