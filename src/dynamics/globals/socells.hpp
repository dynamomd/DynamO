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

#ifndef CGSOCells_HPP
#define CGSOCells_HPP

#include "global.hpp"
#include "../../extcode/mathtemplates.hpp"
#include "../../datatypes/vector.hpp"
#include "../../simulation/particle.hpp"
#include <vector>

class CGSOCells: public Global
{
public:
  CGSOCells(const XMLNode&, DYNAMO::SimData*);

  CGSOCells(DYNAMO::SimData*, const std::string&);

  virtual ~CGSOCells() {}

  virtual Global* Clone() const 
  { 
    return new CGSOCells(*this); 
  }

  virtual GlobalEvent getEvent(const Particle &) const;

  virtual void runEvent(const Particle&, const double) const;

  virtual void initialise(size_t);

  virtual void operator<<(const XMLNode&);

  virtual void outputXML(xml::XmlStream& XML) const;

protected:

  CVector<int> cellCount;
  Vector  cellDimension;
  size_t cuberootN;
};

#endif
