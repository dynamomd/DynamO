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

#pragma once

#include <dynamo/species/species.hpp>
#include <dynamo/BC/LEBC.hpp>
#include <magnet/xmlwriter.hpp>
#include <memory>

namespace Gtk {
  class VBox;
  class RadioButton;
}

namespace dynamo {
  class SpPoint: public Species
  {
  public:
    template<class T1>
    SpPoint(dynamo::Simulation* sim, IDRange* r, T1 nmass, std::string nName, unsigned int ID):
      Species(sim, "SpPoint", r, nmass, nName, ID)
    {}
  
    SpPoint(const magnet::xml::Node& XML, dynamo::Simulation* nSim, unsigned int nID):
      Species(nSim, "", NULL, 0, "", nID)
    { operator<<(XML); }

    virtual void initialise() {}

    virtual void operator<<(const magnet::xml::Node& XML);

    virtual double getScalarMomentOfInertia(size_t ID) const 
    { M_throw() << "Species has no intertia"; }

    virtual double getParticleKineticEnergy(size_t ID)  const {
      const double mass = getMass(ID);
      const Particle& part = Sim->particles[ID];

      if (std::isinf(mass))
	return 0;      
      
      if (std::dynamic_pointer_cast<BCLeesEdwards>(Sim->BCs))
	return 0.5 * static_cast<const BCLeesEdwards&>(*Sim->BCs).getPeculiarVelocity(part).nrm2() * mass;
      else
	return 0.5 * part.getVelocity().nrm2() * mass;
    }
    
  protected:
    virtual void outputXML(magnet::xml::XmlStream& XML) const;
  };
}
