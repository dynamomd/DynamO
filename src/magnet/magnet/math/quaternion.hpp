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
#include <magnet/math/matrix.hpp>

namespace magnet {
  namespace math {
    
    /*! \brief A quaternion class.
     */
    class Quaternion
    {
      Vector _imaginary;
      double _real;

    public:
      Quaternion(): _imaginary(0,0,0), _real(1) {}

      Quaternion(double r, double i, double j, double k):
	_imaginary(i,j,k), _real(r) {}

      Quaternion(double real, Vector imaginary):
	_imaginary(imaginary), _real(real)
      {}

      /*! \brief Returns the default unrotated director.
	
	When quaternions are used to store an orientation, they
	actually encode a rotation from a reference director to the
	encoded direction. This function returns a default unit
	reference vector.
       */
      static Vector initialDirector() { return Vector(0,0,1); }

      /*! \brief Create a quaternion from the cosine of the rotation
          angle and a rotation axis.
       */
      static Quaternion fromCosAngleAxis(double cosangle, Vector axis)
      {
	Quaternion retval(cosangle, axis);
	retval.normalise();
	
	//The current quaternion represents a rotation which is twice
	//the required angle.
	//Perform a half angle converison
	retval.halfRotation();
	return retval;
      }
      
      /*! \brief Create a quaternion from a rotation angle and a
          rotation axis.
       */
      static Quaternion fromAngleAxis(double angle, Vector axis)
      {
	if (angle == 0) return identity();
        double half_angle = 0.5 * angle;
	return Quaternion(std::cos(half_angle), std::sin(half_angle) * axis);
      }
      

      /*! \brief Create a quaternion from the shortest rotation
          between two unit vectors.
	
	The two vectors passed as arguments must be normalised!

	The quaternion is a rotation that takes the from vector into
	the to vector using the shortest arc. This cannot calculate
	the additional possible rotation about the axis of the vector.
       */
      static Quaternion fromToVector(Vector to, Vector from = initialDirector())
      {
	//get the cosine of the angle between the unrotated vector and
	//this vector
	double cosangle = (from | to);

	//Special case of no rotation
	if (cosangle >= 1) return identity();

	//Special case where vec and axis are opposites
	if (cosangle <= -1) return Quaternion(0,1,0,0);
	
	return fromCosAngleAxis(cosangle, from ^ to);
      }

      /*! \brief This calculates a quaternion from a rotation axis,
          where the magnitude of the vector is the angle of rotation.
      */
      static Quaternion fromRotationAxis(Vector axis)
      {
	double angle = axis.nrm();
	axis /= angle + (angle==0);
	return fromAngleAxis(angle, axis);
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

      const Vector& imaginary() const { return _imaginary; }      
      Vector& imaginary() { return _imaginary; }      
      const double& real() const { return _real; }
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

      /*! \brief Rotation of a vector (assuming the vector is already normalised) */
      Vector operator*(const Vector& vec) const {
	Vector img = imaginary();
	return vec + 2.0 * (img ^ ((img ^ vec) + real() * vec));
	//Equivalent but slower
	//return (((*this) * Quaternion(0,vec)) * ((*this).conjugate()))._imaginary;
      }

      Quaternion operator*(const Quaternion& q) const {
	const Quaternion& r(*this);
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

      Matrix toMatrix() const {
	double xx = _imaginary[0] * _imaginary[0];
	double xy = _imaginary[0] * _imaginary[1];
	double xz = _imaginary[0] * _imaginary[2];
	double xw = _imaginary[0] * _real;
	double yy = _imaginary[1] * _imaginary[1];
	double yz = _imaginary[1] * _imaginary[2];
	double yw = _imaginary[1] * _real;
	double zz = _imaginary[2] * _imaginary[2];
	double zw = _imaginary[2] * _real;
	
	return Matrix(1 - 2 * ( yy + zz ), 2 * ( xy - zw ), 2 * ( xz + yw ),
		      2 * ( xy + zw ), 1 - 2 * ( xx + zz ), 2 * ( yz - xw ),
		      2 * ( xz - yw ), 2 * ( yz + xw ), 1 - 2 * ( xx + yy ));
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

namespace dynamo { typedef ::magnet::math::Quaternion Quaternion; }
