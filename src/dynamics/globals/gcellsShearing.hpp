/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifndef CGCellsShearing_HPP
#define CGCellsShearing_HPP

#include "gcells.hpp"
#include "../ranges/1range.hpp"


class CGCellsShearing: public CGCells
{
public:
  CGCellsShearing(const XMLNode&, DYNAMO::SimData*);
  
  CGCellsShearing(DYNAMO::SimData*, const std::string&);
  
  virtual ~CGCellsShearing() {}

  virtual void initialise(size_t);
  
  virtual CGlobal* Clone() const 
  { return new CGCellsShearing(*this); }
  
  virtual void runEvent(const Particle&) const;

  virtual void getParticleNeighbourhood(const Particle&, 
					const nbHoodFunc&) const;

  void getExtraLEParticleNeighbourhood(const Particle& part,
				       const nbHoodFunc& func) const;

protected:
  virtual void outputXML(xmlw::XmlStream&) const;
};

#endif
