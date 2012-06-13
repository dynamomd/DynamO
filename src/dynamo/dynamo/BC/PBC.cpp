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

#include <dynamo/BC/PBC.hpp>
#include <dynamo/simulation.hpp>
#include <magnet/xmlwriter.hpp>

namespace dynamo {
  BCPeriodic::BCPeriodic(const dynamo::Simulation* tmp):
    BoundaryCondition(tmp, "RPBC")
  {
    Sim = tmp;
  }

  void 
  BCPeriodic::applyBC(Vector & pos) const
  { 
    for (size_t n = 0; n < NDIM; ++n)
      pos[n] -= Sim->primaryCellSize[n] *
	lrint(pos[n]/Sim->primaryCellSize[n]);    
  }

  void 
  BCPeriodic::applyBC(Vector & pos, Vector&) const
  {
    for (size_t n = 0; n < NDIM; ++n)
      pos[n] -= Sim->primaryCellSize[n] *
	lrint(pos[n] / Sim->primaryCellSize[n]);    
  }

  void 
  BCPeriodic::applyBC(Vector  &pos, const double&) const 
  {
    for (size_t n = 0; n < NDIM; ++n)
      pos[n] -= Sim->primaryCellSize[n] *
	lrint(pos[n] / Sim->primaryCellSize[n]);    
  }

  void 
  BCPeriodic::outputXML(magnet::xml::XmlStream &XML) const
  {
    XML << magnet::xml::attr("Type") << "PBC";
  }

  void 
  BCPeriodic::operator<<(const magnet::xml::Node&) 
  {}

  BCPeriodicExceptX::BCPeriodicExceptX(const dynamo::Simulation* tmp):
    BoundaryCondition(tmp, "NoXPBC")
  { Sim = tmp; }

  void 
  BCPeriodicExceptX::outputXML(magnet::xml::XmlStream &XML) const
  {
    XML << magnet::xml::attr("Type") << "NoXPBC";
  }

  void 
  BCPeriodicExceptX::operator<<(const magnet::xml::Node&) 
  {}

  void 
  BCPeriodicExceptX::applyBC(Vector & pos) const
  { 
    double x = pos[0];

    for (size_t n = 0; n < NDIM; ++n)
      pos[n] -= Sim->primaryCellSize[n] *
	lrint(pos[n] / Sim->primaryCellSize[n]);    

    pos[0] = x;
  }
  
  void 
  BCPeriodicExceptX::applyBC(Vector & pos, Vector&) const
  { 
    double x = pos[0];

    for (size_t n = 0; n < NDIM; ++n)
      pos[n] -= Sim->primaryCellSize[n] *
	lrint(pos[n]/Sim->primaryCellSize[n]);    

    pos[0] = x;
  }

  void 
  BCPeriodicExceptX::applyBC(Vector  &pos, const double&) const 
  { 
    double x = pos[0];

    for (size_t n = 0; n < NDIM; ++n)
      pos[n] -= Sim->primaryCellSize[n] *
	lrint(pos[n] / Sim->primaryCellSize[n]);    
  
    applyBC(pos); 

    pos[0] = x;
  }

  BCPeriodicXOnly::BCPeriodicXOnly(const dynamo::Simulation* tmp):
    BoundaryCondition(tmp, "NoXPBC")
  {
    Sim = tmp;
  }

  void 
  BCPeriodicXOnly::outputXML(magnet::xml::XmlStream &XML) const
  {
    XML << magnet::xml::attr("Type") << "OnlyXPBC";
  }

  void 
  BCPeriodicXOnly::operator<<(const magnet::xml::Node&) 
  {}

  void 
  BCPeriodicXOnly::applyBC(Vector & pos) const
  { 
    pos[0] -= Sim->primaryCellSize[0] 
      * lrint(pos[0] / Sim->primaryCellSize[0]);
  }
  
  void 
  BCPeriodicXOnly::applyBC(Vector & pos, Vector&) const
  { 
    applyBC(pos);
  }

  void 
  BCPeriodicXOnly::applyBC(Vector  &pos, const double&) const 
  { 
    applyBC(pos); 
  }
}
