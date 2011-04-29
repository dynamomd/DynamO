/*  dynamo:- Event driven molecular dynamics simulator 
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

class CCapsulesFunc : public CShape {
public:
  CCapsulesFunc(const Vector& nr12, const Vector& nv12,
		const Vector& nw1, const Vector& nw2, 
		const Vector& nu1, const Vector& nu2):
    w1(nw1), w2(nw2), u1(nu1), u2(nu2),
    w12(nw1-nw2), r12(nr12), v12(nv12)
  {}

  void stream(const double& dt)
  {
    u1 = Rodrigues(w1 * dt) * Vector(u1);
    u2 = Rodrigues(w2 * dt) * Vector(u2);
    r12 += v12 * dt;
  }
  
  std::pair<double, double> getCollisionPoints() const
  {
    M_throw() << "Don't use!";
  }
  
  //Distance between 2 particles
   double F_zeroDeriv() const
  // For the moment we will assume only one sided dumbbell
  // so the equation get simpler.
  { return (r12 + (u1 + u2) * L * 0.5).nrm2()  - Diameter * Diameter;}
			      

  double F_firstDeriv() const
  {    
    // Simply chain rule
    return 2.0* (r12 + (u1 + u2) * L * 0.5 ) | (v12 + ((w1 ^ u1) + (w2 ^ u2)) * L * 0.5);
  }

  double F_firstDeriv_max(const double& length) const
  { 
    return  2 * (3 * L + Diameter) * (v12.nrm() + (w1.nrm() + w2.nrm()) * L / 2);
  }

  double F_secondDeriv() const
  {
    return 2.0* ((r12 + u1*L/2.0 + u2*L/2.0)|(- w1.nrm()*w1.nrm()*u1*L/2.0 - w2.nrm()*w2.nrm()*u2*L/2.0 ) 
		 + (v12 + w1^u1*L/2.0 + w2^u2*L/2.0)|(v12 + w1^u1*L/2.0 + w2^u2*L/2.0)) ;
  }

    double F_secondDeriv_max(const double& length) const
    {
      return 2.0* ((2*L)*(+ w1.nrm()*w1.nrm()*L/2.0 + w2.nrm()*w2.nrm()*L/2.0 ) + 
		   (v12.nrm() + w1.nrm()*L/2.0 + w2.nrm()*L/2.0)*(v12.nrm() + w1.nrm()*L/2.0 + w2.nrm()*L/2.0));
    }

    std::pair<double, double> discIntersectionWindow(const double& length) const
    {
      //I think this remains the same exept the length goes to lenght + diameter 
      Vector  Ahat = w1 / w1.nrm();
      double dotproduct = (w1 | w2) / (w2.nrm() * w1.nrm());
      double signChangeTerm = (length / 2.0 + r) * sqrt(1.0 - pow(dotproduct, 2.0));
    
      std::pair<double,double> 
	retVal(((-1.0 * (r12 | Ahat)) - signChangeTerm) / (v12 | Ahat),
	       ((-1.0 * (r12 | Ahat)) + signChangeTerm) / (v12 | Ahat));
  
      if(retVal.second < retVal.first) std::swap(retVal.first, retVal.second);

      return retVal;
    }

    const Vector& getu1() const { return u1; }
    const Vector& getu2() const { return u2; }
    const Vector& getw1() const { return w1; }
    const Vector& getw2() const { return w2; }
    const Vector& getw12() const { return w12; }
    const Vector& getr12() const { return r12; }
    const Vector& getv12() const { return v12; }

    virtual CShape* Clone() const { return new CCapsulesFunc(*this); };

    virtual bool test_root(const double& length) const
    {
      double cp = getCollisionPoints();
    
      return  fabs(cp) < 1e-16 ;
    }
  
  private:
    const Vector& w1;
    const Vector& w2;
    Vector u1;
    Vector u2;
    Vector w12;
    Vector r12;
    Vector v12;
  };
