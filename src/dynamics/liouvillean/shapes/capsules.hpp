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
    double rijdotui = (r12 | u1);
    double rijdotuj = (r12 | u2);
    double uidotuj = (u1 | u2);

    return std::make_pair(- (rijdotui - (rijdotuj * uidotuj)) / (1.0 - uidotuj*uidotuj),
			  (rijdotuj - (rijdotui * uidotuj)) / (1.0 - uidotuj*uidotuj));
  }
  //Distance between the particles
  //for the case of the capsule it is the minimum
  //between the line and the 2 spheres.
  //Maybe a better option would be to take the time for each 
  //component and then take the minimum
  double F_zeroDeriv() const
  // For the moment we will assume only one sided dumbbell
  // so the equation get simpler. 
  { return   (r12 + u1*L/2.0 + u2*L/2.0)| (r12 + u1*L/2.0 + u2*L/2.0)  - Diameter*Diameter;}
			      

  //for the following functions we need to save the value of F
  //to decide which derivative to use. 
  double F_firstDeriv() const
  {    
    // Simply chain rule
    return 2.0* (r12 + u1*L/2.0 + u2*L/2.0)|(v12 + w1^u1*L/2.0 + + w2^u2*L/2.0);
  }

  double F_firstDeriv_max(const double& length) const
  { return  2.0* (2.0*L)*(v12.nrm() + w1.nrm()*L/2.0 + + w2.nrm()*L/2.0);

    double F_secondDeriv() const
    {
      return 2.0* ((r12 + u1*L/2.0 + u2*L/2.0)|(- w1.nrm()*w1.nrm()*u1*L/2.0 - w2.nrm()*w2.nrm()*u2*L/2.0 ) + (v12 + w1^u1*L/2.0 + + w2^u2*L/2.0)|(v12 + w1^u1*L/2.0 + + w2^u2*L/2.0)) ;
    }

    double F_secondDeriv_max(const double& length) const
    {
      return 2.0* ((2*L)*(+ w1.nrm()*w1.nrm()*L/2.0 + w2.nrm()*w2.nrm()*L/2.0 ) + 
		   (v12.nrm() + w1.nrm()*L/2.0 + w2.nrm()*L/2.0)*(v12.nrm() + w1.nrm()*L/2.0 + w2.nrm()*L/2.0));
    }

    std::pair<double, double> discIntersectionWindow(const double& length) const
    {
      Vector  Ahat = w1 / w1.nrm();
      double dotproduct = (w1 | w2) / (w2.nrm() * w1.nrm());
      double signChangeTerm = (length / 2.0) * sqrt(1.0 - pow(dotproduct, 2.0));
    
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
      std::pair<double,double> cp = getCollisionPoints();
    
      return (fabs(cp.first) < length / 2.0 && fabs(cp.second) < length / 2.0);
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
