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
#include <dynamo/2particleEventData.hpp>
#include <dynamo/units/units.hpp>
#include <dynamo/BC/LEBC.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>

namespace dynamo {
  magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream& XML, const Dynamics& g)
  {
    g.outputXML(XML);
    return XML;
  }

  shared_ptr<Dynamics>
  Dynamics::getClass(const magnet::xml::Node& XML, dynamo::Simulation* tmp)
  {
    if (!strcmp(XML.getAttribute("Type"),"Newtonian"))
      return shared_ptr<Dynamics>(new DynNewtonian(tmp));
    if (!strcmp(XML.getAttribute("Type"),"NewtonianGravity"))
      return shared_ptr<Dynamics>(new DynGravity(tmp, XML));
    else if (!strcmp(XML.getAttribute("Type"),"NewtonianMC"))
      return shared_ptr<Dynamics>(new DynNewtonianMC(tmp, XML));
    else 
      M_throw() << XML.getAttribute("Type")
		<< ", Unknown type of Dynamics encountered";
  }


  void 
  Dynamics::initialise()
  {
    streamFreq = 10 * Sim->N;
    
    if (hasOrientationData())
      {
	double sumEnergy(0.0);
	BOOST_FOREACH(const Particle& part, Sim->particleList)  
	  sumEnergy += Sim->species[part].getScalarMomentOfInertia(part.getID())
	  * orientationData[part.getID()].angularVelocity.nrm2();
      
	//Check if any of the species are overridden
	bool hasInertia(false);
	BOOST_FOREACH(const shared_ptr<Species>& spec, Sim->species)
	  if (std::tr1::dynamic_pointer_cast<SpInertia>(spec))
	    hasInertia = true;

	if (!hasInertia)
	  M_throw() << "No species have inertia, yet the particles have orientational degrees of freedom set!";

	sumEnergy *= 0.5 / Sim->units.unitEnergy();
  
	dout << "System Rotational Energy " << sumEnergy
	     << "\nRotational kT " << sumEnergy / Sim->N << std::endl;
      }
  }

  PairEventData 
  Dynamics::runLineLineCollision(const IntEvent&,
				    const double&, const double&) const
  { M_throw() << "Not implemented for this Dynamics."; }

  std::pair<bool, double>
  Dynamics::getLineLineCollision(const double, 
				    const Particle&, const Particle&,
				    double) const
  { M_throw() << "Not implemented for this Dynamics."; }

  PairEventData 
  Dynamics::runOffCenterSphereOffCenterSphereCollision(const IntEvent&,
							  const double&, const double&,const double&) const
  { M_throw() << "Not implemented for this Dynamics."; }

  bool
  Dynamics::getOffCenterSphereOffCenterSphereCollision(const double, const double,  
							  const Particle&, const Particle&,
							  const double) const
  { M_throw() << "Not implemented for this Dynamics."; }

  double 
  Dynamics::getPBCSentinelTime(const Particle&, const double&) const
  { M_throw() << "Not implemented for this Dynamics."; }

  void 
  Dynamics::loadParticleXMLData(const magnet::xml::Node& XML)
  {
    dout << "Loading Particle Data" << std::endl;

    bool outofsequence = false;  
  
    for (magnet::xml::Node node = XML.getNode("ParticleData").fastGetNode("Pt"); 
	 node.valid(); ++node)
      {
	if (!node.hasAttribute("ID")
	    || node.getAttribute("ID").as<size_t>() != Sim->particleList.size())
	  outofsequence = true;
      
	Particle part(node, Sim->particleList.size());
	part.getVelocity() *= Sim->units.unitVelocity();
	part.getPosition() *= Sim->units.unitLength();
	Sim->particleList.push_back(part);
      }

    if (outofsequence)
      dout << "Particle ID's out of sequence!\n"
	   << "This can result in incorrect capture map loads etc.\n"
	   << "Erase any capture maps in the configuration file so they are regenerated." << std::endl;

    Sim->N = Sim->particleList.size();

    dout << "Particle count " << Sim->N << std::endl;

    if (XML.getNode("ParticleData").hasAttribute("OrientationData"))
      {
	orientationData.resize(Sim->N);
	size_t i(0);
	for (magnet::xml::Node node = XML.getNode("ParticleData").fastGetNode("Pt"); 
	     node.valid(); ++node, ++i)
	  {
	    orientationData[i].orientation << node.getNode("U");
	    orientationData[i].angularVelocity << node.getNode("O");
      
	    double oL = orientationData[i].orientation.nrm();
      
	    if (!(oL > 0.0))
	      M_throw() << "Particle ID " << i 
			<< " orientation vector is zero!";
      
	    //Makes the vector a unit vector
	    orientationData[i].orientation /= oL;
	  }
      }
  }

  void 
  Dynamics::outputParticleXMLData(magnet::xml::XmlStream& XML, bool applyBC) const
  {
    XML << magnet::xml::tag("ParticleData");
  
    if (hasOrientationData())
      XML << magnet::xml::attr("OrientationData") << "Y";

    for (size_t i = 0; i < Sim->N; ++i)
      {
	Particle tmp(Sim->particleList[i]);
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
    double energy(0);
    if (std::tr1::dynamic_pointer_cast<BCLeesEdwards>(Sim->BCs))
      {
	const BCLeesEdwards& bc = static_cast<const BCLeesEdwards&>(*Sim->BCs);

	energy += bc.getPeculiarVelocity(part).nrm2()
	  * Sim->species[part].getMass(part.getID());
      }
    else
      energy += part.getVelocity().nrm2()
	* Sim->species[part].getMass(part.getID());
  
    if (hasOrientationData())
      energy += orientationData[part.getID()].angularVelocity.nrm2()
	* Sim->species[part].getScalarMomentOfInertia(part.getID());

    return 0.5 * energy;
  }

  double 
  Dynamics::getSystemKineticEnergy() const
  {
    double sumEnergy(0);

    BOOST_FOREACH(const Particle& part, Sim->particleList)
      sumEnergy += getParticleKineticEnergy(part);

    return sumEnergy;
  }

  void
  Dynamics::rescaleSystemKineticEnergy(const double& scale)
  {
    double scalefactor(sqrt(scale));

    if (std::tr1::dynamic_pointer_cast<BCLeesEdwards>(Sim->BCs))
      {
	const BCLeesEdwards& bc = static_cast<const BCLeesEdwards&>(*Sim->BCs);

	BOOST_FOREACH(Particle& part, Sim->particleList)
	  part.getVelocity() = Vector(bc.getPeculiarVelocity(part) * scalefactor
				      + bc.getStreamVelocity(part));
      }
    else
      {
	double scalefactor(sqrt(scale));
      
	BOOST_FOREACH(Particle& part, Sim->particleList)
	  part.getVelocity() *= scalefactor;
      }

    BOOST_FOREACH(rotData& rdat, orientationData)
      rdat.angularVelocity *= scalefactor;      
  }

  PairEventData 
  Dynamics::parallelCubeColl(const IntEvent& event, 
				const double& e, 
				const double& d,
				const EEventType& eType) const
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

  ParticleEventData 
  Dynamics::runSphereWallCollision(Particle&, 
				      const Vector &,
				      const double&
				      ) const
  {
    M_throw() << "Not Implemented";
  }

  PairEventData 
  Dynamics::RoughSpheresColl(const IntEvent& event, 
				const double& e, 
				const double& et, 
				const double& d2, 
				const EEventType& eType
				) const
  {
    M_throw() << "Not Implemented, you need rotational dynamics";
  }

  ParticleEventData 
  Dynamics::runRoughWallCollision(Particle& part, 
				     const Vector & vNorm,
				     const double& e,
				     const double& et,
				     const double& r
				     ) const
  {
    M_throw() << "Not Implemented, you need rotational dynamics";
  }

  void 
  Dynamics::initOrientations(double ToI)
  {
    orientationData.resize(Sim->particleList.size());
  
    dout << "Initialising the line orientations" << std::endl;

    double factor = ToI * 0.5;

    Vector angVelCrossing;

    for (size_t i = 0; i < Sim->particleList.size(); ++i)
      {
	//Assign the new velocities
	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  orientationData[i].orientation[iDim] = Sim->normal_sampler();
      
	orientationData[i].orientation /= orientationData[i].orientation.nrm();
      
	for (size_t iDim = 0; iDim < NDIM; ++iDim)
	  angVelCrossing[iDim] = Sim->normal_sampler();
      
	orientationData[i].angularVelocity
	  = orientationData[i].orientation ^ angVelCrossing;
      
	orientationData[i].angularVelocity *= Sim->normal_sampler() * factor 
	  / orientationData[i].angularVelocity.nrm();
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
  Dynamics::getCOMPosVel(const Range& particles) const
  {
    if (particles.empty())
      M_throw() << "Cannot calculate the COM position and velocity from an empty Range";
    
    Vector pos = Vector(0,0,0), 
      vel = Vector(0,0,0);

    Vector pos0 = Sim->particleList[*(particles.begin())].getPosition(), 
      vel0 = Sim->particleList[*(particles.begin())].getVelocity();

    double totMass = 0;

    BOOST_FOREACH(size_t ID, particles)
      {
	const Particle& part = Sim->particleList[ID];
	double mass = Sim->species[part].getMass(ID);

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
    
    return std::make_pair(pos, vel);
  }

}
