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
#include <dynamo/dynamics/shapes/shape.hpp>
#include <limits>

namespace dynamo {
  class SFOscillatingPlate : public ShapeFunc {
  public:
    SFOscillatingPlate(const Vector& nvp, const Vector& nnhat,
			  const Vector& nrp, 
			  const double& nt, const double& nDelta,
			  const double& nOmega, const double& nSigma):
      vp(nvp), nhat(nnhat), rp(nrp),
      t(nt), Delta(nDelta), Omega(nOmega), Sigma(nSigma)
    {    
    }

    void stream(const double& dt)
    {
      t += dt;
      rp += vp * dt;
    }
  
    double velnHatWall() const
    {
      double retval = - Delta * Omega * std::sin(Omega * t);
      return  retval;
    }

    double maxWallVel() const
    {    
      return  Delta * Omega;
    }

    Vector wallPosition() const
    {
      return  nhat * wallnHatPosition();
    }

    double wallnHatPosition() const
    {
      return  (Delta * std::cos(Omega * t));
    }

    Vector wallVelocity() const
    {
      return  nhat * velnHatWall();
    }

    double F_zeroDeriv() const
    {
      return (rp | nhat) - ( Sigma + wallnHatPosition());
    }

    void fixFZeroSign(bool sign)
    {
      rp -= (rp | nhat) * nhat;
      rp += nhat * (wallnHatPosition() + Sigma);
      double Fval = F_zeroDeriv();
      size_t loop(1);
      while (sign ? (Fval < 0) : (Fval > 0))
	{
	  rp -= nhat * ((loop++) * std::numeric_limits<double>::epsilon())
	    * Sigma;
	  Fval = F_zeroDeriv();
	}
    }

    double F_zeroDerivFlip() const
    { 
      return ((rp - wallPosition()) | nhat) + Sigma;
    }

    double F_firstDeriv() const
    {    
      return (vp | nhat) - velnHatWall();
    }

    double F_firstDeriv_max() const
    {
      return std::fabs(vp | nhat) + Delta * Omega; 
    }

    double F_secondDeriv() const
    {
      return Delta * Omega * Omega * std::cos(Omega * t);
    }

    double F_secondDeriv_max() const
    {
      return Delta * Omega * Omega;
    }

    virtual bool test_root() const
    {
      return (((vp | nhat) - velnHatWall()) 
	      * ((rp | nhat) - wallnHatPosition())) > 0;
    }

    void flipSigma() { Sigma = -Sigma; }
  
  private:
    const Vector& vp;
    const Vector& nhat;
    Vector rp;
    long double t;
    double Delta;
    double Omega;
    double Sigma;  
  };
}
