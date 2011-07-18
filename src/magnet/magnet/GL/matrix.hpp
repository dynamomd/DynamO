/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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
#include <tr1/array>

namespace magnet {
  namespace GL {
    /*! \brief A 4x4 matrix class for projection/model view matrix
     * math.
     */
    class Matrix : public std::tr1::array<GLfloat, 4 * 4>
    {
      typedef std::tr1::array<GLfloat, 4 * 4> Base;
    public:
      Matrix() {}

      Matrix(const Base& o):
	Base(o)
      {}

      Matrix(const Matrix& o):
	Base(o)
      {}

      /*! \brief Return an identity matrix.
       */
      inline static Matrix identity()
      {
	Base retval = {{1,0,0,0, 
			0,1,0,0,
			0,0,1,0,
			0,0,0,1}};
	return retval;
      }
      
      /*! \brief Return a matrix corresponding to a translation.
       */
      inline static Matrix translate(const Vector& vec)
      {
	Base retval = {{1,0,0,0, 
			0,1,0,0,
			0,0,1,0,
			vec[0],vec[1],vec[2],1}};
	return retval;
      }

      /*! \brief Return a matrix corresponding to a rotation.
       * \param angle The angle of rotation (in degrees).
       * \param axis The axis of rotation.
       */
      inline static Matrix rotate(const GLfloat& angle, const Vector& axis)
      {
	::Matrix m = Rodrigues((angle * M_PI / 180.0f) * axis);

	Base retval =
	  {{m(0,0),m(1,0),m(2,0),0,
	    m(0,1),m(1,1),m(2,1),0,
	    m(0,2),m(1,2),m(2,2),0,
	    0,0,0,1}};
	
	return retval;	
      }

      inline Matrix operator*=(const Matrix& om)
      {
	(*this) = (*this) * om; 
	return *this;
      }

      inline Matrix operator*(const Matrix& om) const
      {
	Matrix retval;
	for (size_t i(0); i < 4; ++i)
	  for (size_t j(0); j < 4; ++j)
	    {
	      GLfloat sum(0);
	      for (size_t k(0); k < 4; ++k)
		sum += operator[](k * 4 + j) * om[i * 4 + k];
	      retval[i * 4 + j] = sum;
	    }
	return retval;
      }
      
      inline void transpose()
      {
	for (size_t i(0); i < 4; ++i)
	  for (size_t j(i+1); j < 4; ++j)
	    std::swap((*this)[4*i+j], (*this)[4*j+i]);
      }


      Matrix inverse() const
      {
	Matrix result;

	GLfloat tmp[12];											//temporary pair storage

	//calculate pairs for first 8 elements (cofactors)
	tmp[0] = (*this)[10] * (*this)[15];
	tmp[1] = (*this)[11] * (*this)[14];
	tmp[2] = (*this)[9] * (*this)[15];
	tmp[3] = (*this)[11] * (*this)[13];
	tmp[4] = (*this)[9] * (*this)[14];
	tmp[5] = (*this)[10] * (*this)[13];
	tmp[6] = (*this)[8] * (*this)[15];
	tmp[7] = (*this)[11] * (*this)[12];
	tmp[8] = (*this)[8] * (*this)[14];
	tmp[9] = (*this)[10] * (*this)[12];
	tmp[10] = (*this)[8] * (*this)[13];
	tmp[11] = (*this)[9] * (*this)[12];

	//calculate first 8 elements (cofactors)
	result[0] = tmp[0]*(*this)[5] + tmp[3]*(*this)[6] + tmp[4]*(*this)[7]-tmp[1]*(*this)[5] - tmp[2]*(*this)[6] - tmp[5]*(*this)[7];
	result[1] = tmp[1]*(*this)[4] + tmp[6]*(*this)[6] + tmp[9]*(*this)[7]-tmp[0]*(*this)[4] - tmp[7]*(*this)[6] - tmp[8]*(*this)[7];
	result[2] = tmp[2]*(*this)[4] + tmp[7]*(*this)[5] + tmp[10]*(*this)[7]-tmp[3]*(*this)[4] - tmp[6]*(*this)[5] - tmp[11]*(*this)[7];
	result[3] = tmp[5]*(*this)[4] + tmp[8]*(*this)[5] + tmp[11]*(*this)[6]-tmp[4]*(*this)[4] - tmp[9]*(*this)[5] - tmp[10]*(*this)[6];
	result[4] = tmp[1]*(*this)[1] + tmp[2]*(*this)[2] + tmp[5]*(*this)[3]-tmp[0]*(*this)[1] - tmp[3]*(*this)[2] - tmp[4]*(*this)[3];
	result[5] = tmp[0]*(*this)[0] + tmp[7]*(*this)[2] + tmp[8]*(*this)[3]-tmp[1]*(*this)[0] - tmp[6]*(*this)[2] - tmp[9]*(*this)[3];
	result[6] = tmp[3]*(*this)[0] + tmp[6]*(*this)[1] + tmp[11]*(*this)[3]-tmp[2]*(*this)[0] - tmp[7]*(*this)[1] - tmp[10]*(*this)[3];
	result[7] = tmp[4]*(*this)[0] + tmp[9]*(*this)[1] + tmp[10]*(*this)[2]-tmp[5]*(*this)[0] - tmp[8]*(*this)[1] - tmp[11]*(*this)[2];

	//calculate pairs for second 8 elements (cofactors)
	tmp[0] = (*this)[2]*(*this)[7];
	tmp[1] = (*this)[3]*(*this)[6];
	tmp[2] = (*this)[1]*(*this)[7];
	tmp[3] = (*this)[3]*(*this)[5];
	tmp[4] = (*this)[1]*(*this)[6];
	tmp[5] = (*this)[2]*(*this)[5];
	tmp[6] = (*this)[0]*(*this)[7];
	tmp[7] = (*this)[3]*(*this)[4];
	tmp[8] = (*this)[0]*(*this)[6];
	tmp[9] = (*this)[2]*(*this)[4];
	tmp[10] = (*this)[0]*(*this)[5];
	tmp[11] = (*this)[1]*(*this)[4];

	//calculate second 8 elements (cofactors)
	result[8 ] = tmp[0]*(*this)[13] + tmp[3]*(*this)[14] + tmp[4]*(*this)[15]-tmp[1]*(*this)[13] - tmp[2]*(*this)[14] - tmp[5]*(*this)[15];
	result[9 ] = tmp[1]*(*this)[12] + tmp[6]*(*this)[14] + tmp[9]*(*this)[15]-tmp[0]*(*this)[12] - tmp[7]*(*this)[14] - tmp[8]*(*this)[15];
	result[10] = tmp[2]*(*this)[12] + tmp[7]*(*this)[13] + tmp[10]*(*this)[15]-tmp[3]*(*this)[12] - tmp[6]*(*this)[13] - tmp[11]*(*this)[15];
	result[11] = tmp[5]*(*this)[12] + tmp[8]*(*this)[13] + tmp[11]*(*this)[14]-tmp[4]*(*this)[12] - tmp[9]*(*this)[13] - tmp[10]*(*this)[14];
	result[12] = tmp[2]*(*this)[10] + tmp[5]*(*this)[11] + tmp[1]*(*this)[9]-tmp[4]*(*this)[11] - tmp[0]*(*this)[9] - tmp[3]*(*this)[10];
	result[13] = tmp[8]*(*this)[11] + tmp[0]*(*this)[8] + tmp[7]*(*this)[10]-tmp[6]*(*this)[10] - tmp[9]*(*this)[11] - tmp[1]*(*this)[8];
	result[14] = tmp[6]*(*this)[9] + tmp[11]*(*this)[11] + tmp[3]*(*this)[8]-tmp[10]*(*this)[11] - tmp[2]*(*this)[8] - tmp[7]*(*this)[9];
	result[15] = tmp[10]*(*this)[10] + tmp[4]*(*this)[8] + tmp[9]*(*this)[9]-tmp[8]*(*this)[9] - tmp[11]*(*this)[10] - tmp[5]*(*this)[8];

	// calculate determinant
	GLfloat det = (*this)[0]*result[0]+(*this)[1]*result[1]+(*this)[2]*result[2]
	  +(*this)[3]*result[3];

	for (size_t i(0); i < 4*4; ++i)
	  result[i] /= det;

	if(det==0.0f)
	  for (size_t i(0); i < 4; ++i)
	    for (size_t j(0); j < 4; ++j)
	      result[4*i+j] = (i==j) ? 1 : 0;

	//Now we need to transpose the matrix as we copy out
	result.transpose();
	return result;
      }
    };
  }
}
