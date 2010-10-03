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

#ifndef LNewtonianMC_H
#define LNewtonianMC_H

#include "NewtonL.hpp"
#include <boost/unordered_map.hpp>

class LNewtonianMC: public LNewtonian
{
public:
  LNewtonianMC(DYNAMO::SimData* tmp, const XMLNode&);

  //Pair particle dynamics
  virtual PairEventData SphereWellEvent(const IntEvent&, const double&, 
					 const double&) const;
  virtual NEventData multibdyWellEvent(const CRange&, const CRange&, 
					   const double&, const double&, 
					   EEventType&) const;

  //Cloning
  virtual Liouvillean* Clone() const { return new LNewtonianMC(*this); }

  virtual void initialise();

protected:
  virtual void outputXML(xml::XmlStream& ) const;

  boost::unordered_map<int, double> _MCEnergyPotential; 

  double EnergyPotentialStep;
  
};
#endif
