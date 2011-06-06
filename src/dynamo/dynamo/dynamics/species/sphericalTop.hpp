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

#pragma once

#include "inertia.hpp"

class SpSphericalTop: public SpInertia
{
public:
  SpSphericalTop(dynamo::SimData*, CRange*, double nMass, std::string nName, 
		 unsigned int ID, double iC, std::string nIName="Bulk");
  
  SpSphericalTop(const magnet::xml::Node&, dynamo::SimData*, unsigned int ID);

  virtual Species* Clone() const { return new SpSphericalTop(*this); }

  virtual double getScalarMomentOfInertia(size_t ID) const 
  { return inertiaConstant * getMass(ID); }

  virtual void operator<<(const magnet::xml::Node&);

protected:

  virtual void outputXML(xml::XmlStream& XML) const { outputXML(XML, "SphericalTop"); }

  void outputXML(xml::XmlStream& XML, std::string type) const;
  
  double inertiaConstant;
};
