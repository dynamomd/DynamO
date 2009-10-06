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

#ifndef CSNBListCompressionFix_HPP
#define CSNBListCompressionFix_HPP

#include "system.hpp"
#include "../../datatypes/vector.hpp"

class CSNBListCompressionFix: public CSystem
{
public:
  CSNBListCompressionFix(DYNAMO::SimData*, Iflt, size_t);
  
  virtual CSystem* Clone() const { return new CSNBListCompressionFix(*this); }

  virtual void runEvent() const;

  virtual void initialise(size_t);
  
  virtual void operator<<(const XMLNode&) {}
  
protected:
  virtual void outputXML(xmlw::XmlStream&) const {}
  
  Iflt growthRate;
  size_t cellID;
};

#endif
