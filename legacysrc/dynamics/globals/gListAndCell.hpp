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

#ifndef CGListAndCell_HPP
#define CGListAndCell_HPP

#include "gcells.hpp"
#include "../ranges/1range.hpp"


class CGListAndCell: public CGCells
{
public:
  CGListAndCell(const XMLNode&, DYNAMO::SimData*);

  CGListAndCell(DYNAMO::SimData*, const std::string&);

  virtual ~CGListAndCell() {}

  virtual Global* Clone() const 
  { return new CGListAndCell(*this); }

  virtual void initialise(size_t);

  virtual void getParticleNeighbourhood(const Particle&, 
					const nbHoodFunc&) const;

  virtual void operator<<(const XMLNode&);

  virtual Iflt getMaxInteractionLength() const;

protected:
  virtual void outputXML(xml::XmlStream&) const;

  ClonePtr<CRange> largestParticles;
};

#endif
