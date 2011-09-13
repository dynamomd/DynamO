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

#include <dynamo/dynamics/BC/LEBC.hpp>
#include <dynamo/dynamics/interactions/intEvent.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <cmath>

namespace dynamo {
  BCLeesEdwards::BCLeesEdwards(const dynamo::SimData* tmp):
    BoundaryCondition(tmp, "LEBC"),
    _dxd(0) ,
    _shearRate(1)
  {
    Sim = tmp;
    dout << " Lee's Edwards BC loaded" << std::endl; 
  }

  BCLeesEdwards::BCLeesEdwards(const magnet::xml::Node& XML, 
			       const dynamo::SimData* tmp):
    BoundaryCondition(tmp, "LEBC"),
    _dxd(0) 
  {
    Sim = tmp;
    operator<<(XML);
    dout << "Lee's Edwards BC loaded\n"
	 << "DXD = " << _dxd
	 << "Shear Rate = " 
	 << _shearRate * Sim->dynamics.units().unitTime() << std::endl;
  }

  void 
  BCLeesEdwards::outputXML(magnet::xml::XmlStream &XML) const
  {
    XML << magnet::xml::attr("Type") << "LE"
	<< magnet::xml::attr("DXD") << _dxd / Sim->dynamics.units().unitLength()
	<< magnet::xml::attr("Rate") << _shearRate * Sim->dynamics.units().unitTime();
  }

  void 
  BCLeesEdwards::operator<<(const magnet::xml::Node& XML)
  { 
    try 
      {
	if (XML.hasAttribute("DXD"))
	  _dxd = XML.getAttribute("DXD").as<double>();
	_dxd *= Sim->dynamics.units().unitLength();
      
	if (XML.hasAttribute("Rate"))
	  _shearRate = XML.getAttribute("Rate").as<double>();      
	_shearRate /= Sim->dynamics.units().unitTime();
      }
    catch (boost::bad_lexical_cast &)
      { M_throw() << "Failed a lexical cast in LEBC"; }
  }

  void 
  BCLeesEdwards::applyBC(Vector& pos) const
  {
    //Shift the x distance due to the Lee's Edwards conditions
    pos[0] -= rint(pos[1] / Sim->primaryCellSize[1]) * _dxd;
  
    for (size_t n = 0; n < NDIM; ++n)
      pos[n] -= Sim->primaryCellSize[n] *
	lrint(pos[n] / Sim->primaryCellSize[n]);    
  }

  void 
  BCLeesEdwards::applyBC(Vector  &pos, Vector &vel) const 
  {
    //Shift the x distance due to the Lee's Edwards conditions
    pos[0] -= rint(pos[1] / Sim->primaryCellSize[1]) * _dxd;
  
    //Adjust the velocity due to the box shift
    vel[0] -= rint(pos[1] / Sim->primaryCellSize[1]) 
      * _shearRate * Sim->primaryCellSize[1];
  
    for (size_t n = 0; n < NDIM; ++n)
      pos[n] -= Sim->primaryCellSize[n] *
	lrint(pos[n] / Sim->primaryCellSize[n]);    
  }

  void 
  BCLeesEdwards::applyBC(Vector& posVec, const double& dt) const 
  { 
    double localdxd = _dxd + dt * _shearRate * Sim->primaryCellSize[1];
  
    //Shift the x distance due to the Lee's Edwards conditions
    posVec[0] -= rint(posVec[1] / Sim->primaryCellSize[1]) * localdxd;
  
    for (size_t n = 0; n < NDIM; ++n)
      posVec[n] -= Sim->primaryCellSize[n] *
	lrint(posVec[n] / Sim->primaryCellSize[n]);    
  }

  void 
  BCLeesEdwards::update(const double& dt) 
  {
    //Shift the boundary of the system v_box = \gamma*L
    _dxd += dt * _shearRate * Sim->primaryCellSize[1];
  
    //PBC for the shift to keep accuracy?
    _dxd -= floor(_dxd/Sim->primaryCellSize[0])*Sim->primaryCellSize[0];
  }

  Vector
  BCLeesEdwards::getStreamVelocity(const Particle& part) const
  { return Vector(part.getPosition()[1] * _shearRate, 0, 0); }

  Vector
  BCLeesEdwards::getPeculiarVelocity(const Particle& part) const
  { return part.getVelocity() - getStreamVelocity(part); }
}
