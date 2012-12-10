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

#pragma once
#include <dynamo/particle.hpp>
#include <dynamo/BC/BC.hpp>
#include <dynamo/simulation.hpp>
#include <magnet/math/matrix.hpp>

namespace dynamo {
  class SFLines 
  {
  public:
    SFLines(const Vector& nr12, const Vector& nv12,
	    const Vector& nw1, const Vector& nw2, 
	    const Vector& nu1, const Vector& nu2,
	    const double& length):
      w1(nw1), w2(nw2), u1(nu1), u2(nu2),
      w12(nw1 - nw2), r12(nr12), v12(nv12),
      _length(length)
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
    
    template<size_t deriv> 
    double eval() const
    {
      switch (deriv)
	{
	case 0:
	  return ((u1 ^ u2) | r12);
	case 1:
	  return ((u1 | r12) * (w12 | u2)) 
	    + ((u2 | r12) * (w12 | u1)) 
	    - ((w12 | r12) * (u1 | u2)) 
	    + (((u1 ^ u2) | v12));
	case 2:
	  return 2.0 
	    * (((u1 | v12) * (w12 | u2)) 
	       + ((u2 | v12) * (w12 | u1))
	       - ((u1 | u2) * (w12 | v12)))
	    - ((w12 | r12) * (w12 | (u1 ^ u2))) 
	    + ((u1 | r12) * (u2 | (w1 ^ w2))) 
	    + ((u2 | r12) * (u1 | (w1 ^ w2)))
	    + ((w12 | u1) * (r12 | (w2 ^ u2)))
	    + ((w12 | u2) * (r12 | (w1 ^ u1))); 
	default:
	  M_throw() << "Invalid access";
	}
    }
    
    template<size_t deriv> 
    double max() const
    {
      switch (deriv)
	{
	case 1:
	  return _length * w12.nrm() + v12.nrm();
	case 2:
	  return w12.nrm() 
	    * ((2 * v12.nrm()) + (_length * (w1.nrm() + w2.nrm())));
	default:
	  M_throw() << "Invalid access";
	}
    }

    std::pair<double, double> discIntersectionWindow() const
    {
      Vector  Ahat = w1 / w1.nrm();
      double dotproduct = (w1 | w2) / (w2.nrm() * w1.nrm());
      double signChangeTerm = (_length / 2.0) * sqrt(1.0 - pow(dotproduct, 2.0));
    
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

    virtual bool test_root() const
    {
      std::pair<double,double> cp = getCollisionPoints();
    
      return (fabs(cp.first) < _length / 2.0 && fabs(cp.second) < _length / 2.0);
    }
  
  private:
    const Vector& w1;
    const Vector& w2;
    Vector u1;
    Vector u2;
    Vector w12;
    Vector r12;
    Vector v12;

    const double _length;
  };
}
