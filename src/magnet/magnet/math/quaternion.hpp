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
      Quaternion():
	_real(1), _imaginary(0,0,0)
      {}

      Quaternion(double r, double i, double j, double k):
	_real(r),
	_imaginary(i,j,k)
      {}

      Quaternion(double real, Vector imaginary):
	_real(real),
	_imaginary(imaginary)
      {}

      /*! \brief Returns the unrotated director.
	
	When quaternions are used to store an orientation, they
	actually encode a rotation from a reference director to the
	encoded direction. This function returns the unit reference vector
	used in the fromOrientation() function.
       */
      static Vector initialDirector() { return Vector(0,0,1); }

      /*! \brief Create a quaternion from a rotation vector.
	
	The quaternion is a rotation that takes the initialDirector()
	vector into the vector passed as an argument. This cannot take
	into account the additional rotation about the axis of the
	vector.
       */
      static Quaternion fromOrientation(const Vector& vec)
      {
	//Calculate the parameters of the rotation
	double vecnrm = vec.nrm();
	//get the cosine of the angle between the unrotated vector and
	//this vector
	double cosangle = (vec | initialDirector()) / vecnrm;

	if ((vecnrm == 0) || (cosangle == 1))
	  //Special case of no rotation
	  return Quaternion(1,0,0,0);

	if (cosangle == -1)
	  //Special case where vec and axis are opposites
	  return Quaternion(0,1,0,0);
	
	//Calculate the rotation axis and store in the imaginary part
	//of the quaternion
	Quaternion retval(cosangle, (vec ^ initialDirector()) / vecnrm);
	retval.normalise();
	
	//The current quaternion represents a rotation which is twice
	//the required angle.
	//Perform a half angle converison
	retval.halfRotation();
	return retval;
      }

      /*! \brief This calculates a quaternion from a rotation axis,
          where the magnitude of the vector is the angle of rotation.
	  
	  This definition does not quite align with other examples on
	  the web, this is because those definitions have an opposite
	  handedness. This form is coherent with the Matrix Rodrigues
	  function.
       */
      static Quaternion fromRotationAxis(const Vector& vec)
      {
	double angle = vec.nrm();

	if (angle == 0)
	  //Special case of no rotation
	  return Quaternion(1,0,0,0);

	double half_angle = angle * 0.5;
	double sin_half_angle = std::sin(half_angle);
	double cos_half_angle = std::cos(half_angle);

	return Quaternion(cos_half_angle, -vec * (sin_half_angle / angle));
      }

      /*! \brief Returns an identity quaternion.
       */
      static Quaternion identity() {
	return Quaternion(1, Vector(0,0,0));
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
	return _imaginary[i - 1];
      }

      Vector imaginary() const { return _imaginary; }      
      Vector& imaginary() { return _imaginary; }      
      double real() const { return _real; }
      double& real() { return _real; }

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
	Vector uv = _imaginary ^ vec;
	Vector uuv = _imaginary ^ uv;
	uv *= 2.0 * _real;
	uuv *= 2.0;
	
	//return vec + uv + uuv;
	return vec + 2.0 * (((vec ^ _imaginary) + _real * vec) ^ _imaginary);
      }

      Quaternion operator*(const Quaternion& r) const {
	const Quaternion& q(*this);
	return Quaternion(r[0] * q[0] - r[1] * q[1] - r[2] * q[2] - r[3] * q[3],
			  r[0] * q[1] + r[1] * q[0] + r[2] * q[3] - r[3] * q[2],
			  r[0] * q[2] + r[2] * q[0] + r[3] * q[1] - r[1] * q[3],
			  r[0] * q[3] + r[3] * q[0] + r[1] * q[2] - r[2] * q[1]);
      }

      Quaternion conjugate() const {
	return Quaternion(_real, -_imaginary);
      }

      Quaternion inverse() const {
	const double inv_nrm2 = 1.0 / nrm2();
	return Quaternion(_real * inv_nrm2, - _imaginary * inv_nrm2);
      }

      std::string toString() const {
	std::ostringstream os;
	os << "[" << _real << "," << _imaginary.toString() << "]";
	return os.str();
      }
    };

    inline magnet::xml::XmlStream& operator<<(magnet::xml::XmlStream& XML, const Quaternion & q)
    {
      XML << q.imaginary() << magnet::xml::attr("w") << q.real();
      return XML;
    }

    inline Quaternion& operator<<(Quaternion& q, const magnet::xml::Node& XML)
    {
      q.imaginary() << XML;
      q.real() = XML.getAttribute("w").as<double>();
      return q;
    }

  }
}
