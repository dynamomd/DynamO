/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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
#include <dynamo/simulation/particle.hpp>
#include <dynamo/dynamics/BC/BC.hpp>

#include <dynamo/base/is_simdata.hpp>
#include <dynamo/dynamics/liouvillean/shapes/shape.hpp>
#include <iomanip>

namespace dynamo {
  class SFDumbbells : public ShapeFunc {
  public:
    SFDumbbells(const Vector& nr12, const Vector& nv12,
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
    {
      M_throw() << "Don't need this function!";
    }
  
    //Distance between 2 particles
    double F_zeroDeriv() const
    // For the moment we will assume only one sided dumbbell
    // so the equation get simpler.
    { return (r12 + (u1 + u2) * L * 0.5).nrm2()  - diameter * diameter;}
  
    double F_firstDeriv() const
    {    
      // Simply chain rule
      return 2 * (r12 + (u1 + u2) * L / 2) | (v12 + ((w1 ^ u1) + (w2 ^ u2)) * L / 2);
    }

    double F_firstDeriv_max() const
    { 
      return 2 * (3 * L + diameter) * (v12.nrm() + (w1.nrm() + w2.nrm()) * L / 2);
    }


    double F_secondDeriv() const
    {
      return 2 * (((r12 + (u1 + u2) * L / 2) | (-w1.nrm2() * u1 - w2.nrm2() * u2)) * L / 2
		  + (v12 + ((w1 ^ u1) + (w2 ^ u2)) * L / 2).nrm2());
    }
  
    double F_secondDeriv_max() const
    {
      double term1 = v12.nrm() + (w1.nrm() + w2.nrm()) * L / 2;
      return 2 * ((L / 2) * (3 * L + diameter) * (w1.nrm2() + w2.nrm2()) 
		  + term1 * term1);
    }
  
    const Vector& getu1() const { return u1; }
    const Vector& getu2() const { return u2; }
    const Vector& getw1() const { return w1; }
    const Vector& getw2() const { return w2; }
    const Vector& getw12() const { return w12; }
    const Vector& getr12() const { return r12; }
    const Vector& getv12() const { return v12; }
  
    virtual bool test_root() const
    {
      return true;
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
}

