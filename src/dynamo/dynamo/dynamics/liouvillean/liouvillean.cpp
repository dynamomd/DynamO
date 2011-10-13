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

#include <dynamo/dynamics/liouvillean/include.hpp>
#include <dynamo/dynamics/species/inertia.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/dynamics/2particleEventData.hpp>
#include <dynamo/dynamics/units/units.hpp>
#include <dynamo/dynamics/BC/LEBC.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>

namespace dynamo {
  magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream& XML, const Liouvillean& g)
  {
    g.outputXML(XML);
    return XML;
  }

  Liouvillean* 
  Liouvillean::loadClass(const magnet::xml::Node& XML, dynamo::SimData* tmp)
  {
    if (!strcmp(XML.getAttribute("Type"),"Newtonian"))
      return new LNewtonian(tmp);
    if (!strcmp(XML.getAttribute("Type"),"NewtonianGravity"))
      return new LNewtonianGravity(tmp, XML);
    else if (!strcmp(XML.getAttribute("Type"),"SLLOD"))
      return new LSLLOD(tmp);
    else if (!strcmp(XML.getAttribute("Type"),"NewtonianMC"))
      return new LNewtonianMC(tmp, XML);
    else 
      M_throw() << XML.getAttribute("Type")
		<< ", Unknown type of Liouvillean encountered";
  }


  void 
  Liouvillean::initialise()
  {
    streamFreq = 10 * Sim->N;
    
    if (hasOrientationData())
      {
	double sumEnergy(0.0);
	BOOST_FOREACH(const Particle& part, Sim->particleList)  
	  sumEnergy += Sim->dynamics.getSpecies(part).getScalarMomentOfInertia(part.getID())
	  * orientationData[part.getID()].angularVelocity.nrm2();
      
	//Check if any of the species are overridden
	bool hasInertia(false);
	BOOST_FOREACH(const std::tr1::shared_ptr<Species>& spec, Sim->dynamics.getSpecies())
	  if (std::tr1::dynamic_pointer_cast<SpInertia>(spec))
	    hasInertia = true;

	if (!hasInertia)
	  M_throw() << "No species have inertia, yet the particles have orientational degrees of freedom set!";

	sumEnergy *= 0.5 / Sim->dynamics.units().unitEnergy();
  
	dout << "System Rotational Energy " << sumEnergy
	     << "\nRotational kT " << sumEnergy / Sim->N << std::endl;
      }
  }

  PairEventData 
  Liouvillean::runLineLineCollision(const IntEvent&,
				    const double&, const double&) const
  { M_throw() << "Not implemented for this Liouvillean."; }

  bool
  Liouvillean::getLineLineCollision(CPDData&, const double&, 
				    const Particle&, const Particle&,
				    double) const
  { M_throw() << "Not implemented for this Liouvillean."; }

  PairEventData 
  Liouvillean::runOffCenterSphereOffCenterSphereCollision(const IntEvent&,
							  const double&, const double&,const double&) const
  { M_throw() << "Not implemented for this Liouvillean."; }

  bool
  Liouvillean::getOffCenterSphereOffCenterSphereCollision(CPDData&, const double&, const double&,  
							  const Particle&, const Particle&
							  ) const
  { M_throw() << "Not implemented for this Liouvillean."; }
  double 
  Liouvillean::getPBCSentinelTime(const Particle&, const double&) const
  { M_throw() << "Not implemented for this Liouvillean."; }

  void 
  Liouvillean::loadParticleXMLData(const magnet::xml::Node& XML)
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
	part.getVelocity() *= Sim->dynamics.units().unitVelocity();
	part.getPosition() *= Sim->dynamics.units().unitLength();
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
  Liouvillean::outputParticleXMLData(magnet::xml::XmlStream& XML, bool applyBC) const
  {
    XML << magnet::xml::tag("ParticleData");
  
    if (hasOrientationData())
      XML << magnet::xml::attr("OrientationData") << "Y";

    for (size_t i = 0; i < Sim->N; ++i)
      {
	Particle tmp(Sim->particleList[i]);
	if (applyBC) 
	  Sim->dynamics.BCs().applyBC(tmp.getPosition(), tmp.getVelocity());
      
	tmp.getVelocity() *= (1.0 / Sim->dynamics.units().unitVelocity());
	tmp.getPosition() *= (1.0 / Sim->dynamics.units().unitLength());
      
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
  Liouvillean::getParticleKineticEnergy(const Particle& part) const
  {
    double energy(0);
    if (Sim->dynamics.BCTypeTest<BCLeesEdwards>())
      {
	const BCLeesEdwards& bc = static_cast<const BCLeesEdwards&>(Sim->dynamics.BCs());

	energy += bc.getPeculiarVelocity(part).nrm2()
	  * Sim->dynamics.getSpecies(part).getMass(part.getID());
      }
    else
      energy += part.getVelocity().nrm2()
	* Sim->dynamics.getSpecies(part).getMass(part.getID());
  
    if (hasOrientationData())
      energy += orientationData[part.getID()].angularVelocity.nrm2()
	* Sim->dynamics.getSpecies(part).getScalarMomentOfInertia(part.getID());

    return 0.5 * energy;
  }

  double 
  Liouvillean::getSystemKineticEnergy() const
  {
    double sumEnergy(0);

    BOOST_FOREACH(const Particle& part, Sim->particleList)
      sumEnergy += getParticleKineticEnergy(part);

    return sumEnergy;
  }

  void
  Liouvillean::rescaleSystemKineticEnergy(const double& scale)
  {
    double scalefactor(sqrt(scale));

    if (Sim->dynamics.BCTypeTest<BCLeesEdwards>())
      {
	const BCLeesEdwards& bc = static_cast<const BCLeesEdwards&>(Sim->dynamics.BCs());

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
  Liouvillean::parallelCubeColl(const IntEvent& event, 
				const double& e, 
				const double& d,
				const Matrix& Rot,
				const EEventType& eType) const
  { M_throw() << "Not Implemented"; }

  std::pair<bool,double>
  Liouvillean::getPointPlateCollision(const Particle& np1, const Vector& nrw0,
				      const Vector& nhat, const double& Delta,
				      const double& Omega, const double& Sigma,
				      const double& t, bool) const
  {
    M_throw() << "Not Implemented";
  }

  ParticleEventData 
  Liouvillean::runOscilatingPlate
  (const Particle& part, const Vector& rw0, const Vector& nhat, double& delta, 
   const double& omega0, const double& sigma, const double& mass, const double& e,
   double& t, bool strongPlate) const
  {
    M_throw() << "Not Implemented";
  }


  double 
  Liouvillean::getCylinderWallCollision(const Particle& part, 
					const Vector& origin, 
					const Vector& norm,
					const double&
					) const
  {
    M_throw() << "Not Implemented";
  }

  ParticleEventData 
  Liouvillean::runCylinderWallCollision(const Particle&, 
					const Vector &,
					const Vector &,
					const double&
					) const
  {
    M_throw() << "Not Implemented";
  }

  ParticleEventData 
  Liouvillean::runSphereWallCollision(const Particle&, 
				      const Vector &,
				      const double&
				      ) const
  {
    M_throw() << "Not Implemented";
  }

  PairEventData 
  Liouvillean::RoughSpheresColl(const IntEvent& event, 
				const double& e, 
				const double& et, 
				const double& d2, 
				const EEventType& eType
				) const
  {
    M_throw() << "Not Implemented, you need rotational dynamics";
  }

  ParticleEventData 
  Liouvillean::runRoughWallCollision(const Particle& part, 
				     const Vector & vNorm,
				     const double& e,
				     const double& et,
				     const double& r
				     ) const
  {
    M_throw() << "Not Implemented, you need rotational dynamics";
  }

  void 
  Liouvillean::initOrientations(double ToI)
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

  std::pair<double, Liouvillean::TriangleIntersectingPart> 
  Liouvillean::getSphereTriangleEvent(const Particle& part, 
				      const Vector & A, 
				      const Vector & B, 
				      const Vector & C,
				      const double d
				      ) const
  {
    M_throw() << "Not implemented";
  }
}
