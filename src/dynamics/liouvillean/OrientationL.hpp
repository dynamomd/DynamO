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

#ifndef CLNOrientation_H
#define CLNOrientation_H

#include "NewtonL.hpp"
#include <vector>
#include "../../datatypes/vector.hpp"

class CLNOrientation: public CLNewton
{
public:
  CLNOrientation(DYNAMO::SimData* Sim, const XMLNode& XML):
    CLNewton(Sim)
  {
    operator<<(XML);
  }

  virtual CLiouvillean* Clone() const { return new CLNOrientation(*this); }

  //Remember to update the signatures in the Liouvillean base clas when you change these
  virtual Iflt getLineLineCollision() const;
  
  //Remember to update the signatures in the Liouvillean base clas when you change these
  virtual C2ParticleData runLineLineCollision() const;

  virtual void operator<<(const XMLNode&);

protected:
  virtual void outputXML(xmlw::XmlStream&) const;

  virtual void streamParticle(CParticle&, const Iflt&) const;

  struct rotData
  {
    CVector<> orientation;
    CVector<> angularMomentum;
  };
  
  std::vector<rotData> orientationData;
};
#endif
