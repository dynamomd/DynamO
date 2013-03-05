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
      double _real;
      Vector _imaginary;

    public:
      Quaternion(double r, double i, double j, double k):
	_real(r),
	_imaginary(i,j,k)
      {}

      /*! \brief Create a quaternion from a rotation vector.
	
	The quaternion is a rotation that takes the vector (0,0,1)
	into the vector passed as an argument.
       */
      Quaternion(const Vector& vec)
      {
	Vector unrotated_director(0,0,1);
	//Calculate the parameters of the rotation
	double vecnrm = vec.nrm();
	//get the cosine of the angle between the unrotated vector and
	//this vector
	double cosangle = (vec | unrotated_director) / vecnrm;

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
	
	//Calculate the rotation axis and store in the imaginary part
	//of the quaternion
	_imaginary = (vec ^ unrotated_director) / vecnrm;
	//Store the angle in the real part
	_real = cosangle;
	normalise();
	
	//The current quaternion represents a rotation which is twice
	//the required angle.
	//Perform a half angle converison
	halfRotation();
      }

      /*! \brief Half the rotation of the current quaternion (assuming
          it is already normalised).
       */
      void halfRotation() {
	_real += 1;
	normalise();
      }
      
      double operator[](size_t i) const {
	if (i == 0) return _real;
	return _imaginary[i - 1];
      }

      double& operator[](size_t i) {
	if (i == 0) return _real;
	return _imaginary[i];
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
	double norm = nrm();
	double inv_norm = 1.0 / (norm + (norm == 0));
	_imaginary *= inv_norm; 
	_real *= inv_norm;
      }

      /*! \brief Rotation of a vector */
      Vector operator*(const Vector& vec) const {
	return vec + 2.0 * (((vec ^ _imaginary) + _real * vec) ^ _imaginary);
      }

      Quaternion operator*(const Quaternion& r) const {
	const Quaternion& q(*this);
	return Quaternion(r[0] * q[0] - r[1] * q[1] - r[2] * q[2] - r[3] * q[3],
			  r[0] * q[1] + r[1] * q[0] - r[2] * q[3] + r[3] * q[2],
			  r[0] * q[2] + r[1] * q[3] + r[2] * q[0] - r[3] * q[1],
			  r[0] * q[3] - r[1] * q[2] + r[2] * q[1] + r[3] * q[0]);
      }

      Quaternion conjugate() const {
	return Quaternion(_real, - _imaginary[0], - _imaginary[1], - _imaginary[2]);
      }

      Quaternion inverse() const {
	const double inv_nrm2 = 1.0 / nrm2();
	return Quaternion(_real * inv_nrm2, - _imaginary[0] * inv_nrm2, - _imaginary[1] * inv_nrm2, - _imaginary[2] * inv_nrm2);
      }      
    };
  }
}
