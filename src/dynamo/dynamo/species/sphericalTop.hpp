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

#include <dynamo/species/inertia.hpp>

namespace dynamo {
  class SpSphericalTop: public SpInertia
  {
  public:
    SpSphericalTop(dynamo::Simulation*, IDRange*, double nMass, std::string nName, unsigned int ID, double iC);
  
    SpSphericalTop(const magnet::xml::Node&, dynamo::Simulation*, unsigned int ID);

    virtual double getScalarMomentOfInertia(size_t ID) const 
    { return inertiaConstant * getMass(ID); }

    virtual void operator<<(const magnet::xml::Node&);

    virtual double getParticleKineticEnergy(size_t ID) const;

    virtual double getDOF() const { return NDIM + 2; }

  protected:

    virtual void outputXML(magnet::xml::XmlStream& XML) const { outputXML(XML, "SphericalTop"); }

    void outputXML(magnet::xml::XmlStream& XML, std::string type) const;
  
    double inertiaConstant;
  };
}
