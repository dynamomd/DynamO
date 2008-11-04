/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "OrientationL.hpp"
#include "../2particleEventData.hpp"

void
CLNOrientation::operator<< (const XMLNode& XML)
{}

void
CLNOrientation::outputXML(xmlw::XmlStream& XML) const
{}

Iflt 
CLNOrientation::getLineLineCollision() const
{
}

C2ParticleData 
CLNOrientation::runLineLineCollision() const
{

}

void 
CLNOrientation::streamParticle(CParticle& part, const Iflt& dt) const
{
 //First hand over to the newtonian dynamics
 CLNewton::streamParticle(part, dt);

 //Now stream the orientation dynamics
}
