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

#pragma once
#include "../../../simulation/particle.hpp"
#include "../../BC/BC.hpp"
#include "../../dynamics.hpp"
#include "../../../base/is_simdata.hpp"
#include "shape.hpp"
#include <iomanip>

class CDumbbellsFunc : public CShape {
public:
  CDumbbellsFunc(const Vector& nr12, const Vector& nv12,
		const Vector& nw1, const Vector& nw2, 
		 const Vector& nu1, const Vector& nu2,const double nl,const double nd):
    w1(nw1), w2(nw2), u1(nu1), u2(nu2),
    w12(nw1-nw2), r12(nr12), v12(nv12),L(nl),diameter(nd)
  {}

  void stream(const double& dt)
  {
    u1 = Rodrigues(w1 * dt) * Vector(u1);
    u2 = Rodrigues(w2 * dt) * Vector(u2);
    r12 += v12 * dt;
  }
  
  Vector getCollisionPoints() const
  { // Marcus, can you check the math?
    M_throw() << "We return " << (r12 + u1 - u2).nrm();
    return 0.5*(r12 + u1 - u2); 
  }
  
  //Distance between 2 particles
   double F_zeroDeriv() const
  // For the moment we will assume only one sided dumbbell
  // so the equation get simpler. 
  { return   ((r12 + u1*L/2.0 - u2*L/2.0)|(r12 + u1*L/2.0 - u2*L/2.0)) - (diameter*diameter);}
			      

  double F_firstDeriv() const
  {    
    // Simply chain rule
    return 2.0* (r12 + u1*L/2.0 - u2*L/2.0)|(v12 + w1^u1*L/2.0 - w2^u2*L/2.0);
  }

  double F_firstDeriv_max() const
  { 
    return  2.0*(3.0*L + diameter)*(v12.nrm() + w1.nrm()*L/2.0 + w2.nrm()*L/2.0);
  }
  
  double F_secondDeriv() const
  {
    return 2.0* (((r12 + u1*L/2.0 - u2*L/2.0)|(- w1.nrm()*w1.nrm()*u1*L/2.0 + w2.nrm()*w2.nrm()*u2*L/2.0 )) 
    		 + ((v12 + w1^u1*L/2.0 - w2^u2*L/2.0)|(v12 + w1^u1*L/2.0 - w2^u2*L/2.0))) ;
   
  }
  
  double F_secondDeriv_max() const
  {
    return 2.0* ((3.0*L + diameter)*(w1.nrm()*w1.nrm()*L/2.0 + w2.nrm()*w2.nrm()*L/2.0 ) + 
   		 (v12.nrm() + w1.nrm()*L/2.0 + w2.nrm()*L/2.0)*(v12.nrm() + w1.nrm()*L/2.0 + w2.nrm()*L/2.0));
  }
  
    
  
  const Vector& getu1() const { return u1; }
  const Vector& getu2() const { return u2; }
  const Vector& getw1() const { return w1; }
  const Vector& getw2() const { return w2; }
  const Vector& getw12() const { return w12; }
  const Vector& getr12() const { return r12; }
  const Vector& getv12() const { return v12; }
  
  virtual CShape* Clone() const { return new CDumbbellsFunc(*this); };
  
  virtual bool test_root() const
  {
    Vector cp = getCollisionPoints();
    M_throw() << "We return " <<(cp.nrm() - diameter < 1e-16);
    return  (cp.nrm() - diameter) < 1e-16 ;
  }
  
private:
  const Vector& w1;
  const Vector& w2;
  Vector u1;
  Vector u2;
  Vector w12;
  Vector r12;
  Vector v12;
  const double L;
  const double diameter;
};
