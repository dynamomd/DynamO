/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef CSSphericalTop_H
#define CSSphericalTop_H

#include "species.hpp"

class CSSphericalTop: public CSpecInertia
{
public:
  CSSphericalTop(DYNAMO::SimData*, CRange*, Iflt nMass, std::string nName, 
		 unsigned int ID, Iflt iC, std::string nIName="Bulk");
  
  CSSphericalTop(const XMLNode&, DYNAMO::SimData*, unsigned int ID);


  virtual CSpecies* Clone() const { return new CSSphericalTop(*this); }

  virtual Iflt getScalarMomentOfInertia() const { return inertiaConstant * mass; }

  virtual void operator<<(const XMLNode&);

protected:

  virtual void outputXML(xmlw::XmlStream&) const;
  
  Iflt inertiaConstant;
};

#endif
