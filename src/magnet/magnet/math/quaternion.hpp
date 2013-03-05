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
#include <magnet/math/vector.hpp>

namespace magnet {
  namespace math {
    
    /*! \brief A quaternion class.
     */
    class Quaternion
    {
      Vector _imaginary;
      double _real;

    public:
      /*! \brief Create a quaternion from a rotation vector.
       */
      Quaternion(const Vector& vec)
      {

	//Calculate the parameters of the rotation
	Vector axis(0,0,1);
	double vecnrm = vec.nrm();
	double cosangle = (vec | axis) / vecnrm;

	if ((vecnrm == 0) || (cosangle == 1))
	  {
	    //Special case of no rotation or zero length vector
	    _imaginary = Vector(0.0, 0.0, 0.0);
	    _real = 1.0;
	    return;
	  }

	if (cosangle == -1)
	  {
	    //Special case where vec and axis are opposites
	    _imaginary = Vector(1.0, 0.0, 0.0);
	    _real = 0.0;
	    return;
	  }
	
	//Calculate the rotation axis and store in the quaternion
	_imaginary = (vec ^ axis) / vecnrm;
	_real = cosangle;
	normalise();

	//Perform a half angle converison
	_imaginary += 1;
	normalise()
      }

      Vector imaginary() const { return _imaginary; }
      
      double real() const { return _real; }

      double nrm2() const {
	return _imaginary.nrm2() + _real * _real;
      }

      double nrm() const {
	return std::sqrt(nrm2());
      }
      
      void normalise() {
	double nrm = nrm();
	double inv_norm = 1.0 / (nrm + (nrm == 0));
	_imaginary *= inv_norm; 
	_real *= inv_norm;
      }

      /*! \brief Rotation of a vector */
      Vector operator*(const Vector& vec) const {
	return vec + 2.0 * (((vec ^ _imaginary) + _real * vec) ^ _imaginary);
      }

      Quaternion operator*(const Quaternion& q2) const {
	return *this;
      }
    };
  }
}
