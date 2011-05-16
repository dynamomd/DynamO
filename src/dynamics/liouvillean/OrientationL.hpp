/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>
    Copyright (C) 2011  Sebastian Gonzalez <tsuresuregusa@gmail.com>

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
#include "NewtonL.hpp"
#include "../../datatypes/vector.hpp"
#include <vector>

class CLinesFunc;
class CDumbbellsFunc;
class CShape;

class LNOrientation: public LNewtonian
{
public:  
  LNOrientation(dynamo::SimData* Sim, const magnet::xml::Node& XML):
    LNewtonian(Sim) {}

  LNOrientation(dynamo::SimData* Sim): LNewtonian(Sim) {}

  virtual void initialise();

  virtual Liouvillean* Clone() const { return new LNOrientation(*this); }

  virtual void loadParticleXMLData(const magnet::xml::Node&);
    
protected:

  virtual void outputXML(xml::XmlStream&) const;

};
