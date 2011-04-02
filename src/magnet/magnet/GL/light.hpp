/*    DYNAMO:- Event driven molecular dynamics simulator 
 *    http://www.marcusbannerman.co.uk/dynamo
 *    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>
 *
 *    This program is free software: you can redistribute it and/or
 *    modify it under the terms of the GNU General Public License
 *    version 3 as published by the Free Software Foundation.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include <limits>
#include <magnet/GL/viewPort.hpp>

namespace magnet {
  namespace GL {        
    class lightInfo: public viewPort
    {
    public:
      //We need a default constructor as viewPorts may be created without GL being initialized
      inline lightInfo():
	viewPort(1,1) {}

      inline lightInfo(Vector position, 
		       Vector lookAtPoint,
		       GLenum lightHandle = GL_LIGHT0,
		       GLfloat fovY = 45.0f,
		       GLfloat zNearDist = 0.01f, GLfloat zFarDist = 10.0f,
		       Vector up = Vector(0,1,0)):
	viewPort(1,1,position, lookAtPoint, fovY, zNearDist, zFarDist, up, 1.0),
	_lightHandle(lightHandle)
      {
	const GLfloat white[] = {1.0f, 1.0f, 1.0f, 1.0f};
	//Setup a bright light
	glLightfv(lightHandle, GL_DIFFUSE, white);
	glLightfv(lightHandle, GL_SPECULAR, white);
      }

      inline void glUpdateLight()
      {
	Vector cameraLocation(getEyeLocation());
	
	GLfloat light_position[4] = {cameraLocation[0], cameraLocation[1], cameraLocation[2], 1.0f};
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
      }

      inline void drawLight()
      {
	glColor3f(1.0f, 1.0f, 1.0f);
	
	GLfloat rotationAngle 
	  = (180.0 / M_PI) * std::acos(Vector(0,0,-1) | _cameraDirection);
	
	Vector RotationAxis = Vector(0,0,-1) ^ _cameraDirection;
	float norm = RotationAxis.nrm();
	RotationAxis /= norm;
	if (norm < std::numeric_limits<double>::epsilon())
	  RotationAxis = Vector(1,0,0);

	Vector cameraLocation(getEyeLocation());
	
	glPushMatrix();
	glTranslatef(cameraLocation[0], cameraLocation[1], cameraLocation[2]);
	glRotatef(rotationAngle, RotationAxis.x, RotationAxis.y, RotationAxis.z);
	
	GLfloat r = 0.05f;
	
	glTranslatef(0.0f,0.0f,0.0025f);
	
	double fovY = getFOVY();
	glutSolidCone(r * std::sin(fovY * M_PI / 360.0f), 
		      r * std::cos(fovY * M_PI / 360.0f), 
		      15, 15);
	glPopMatrix();
      }

      lightInfo& operator=(const viewPort& vp)
      {
	//The FOV and the aspect ratio of the light must be maintained
	size_t width = _width;
	size_t height = _height;
	double fovY = getFOVY();
	viewPort::operator=(vp);
	setFOVY(fovY);
	_width = width;
	_height = height;
	
	buildMatrices();
	return *this; 
      }

      inline void loadShadowTextureMatrix(const viewPort& vp)
      {
	//Build the texture matrix
	glLoadIdentity();
	glTranslatef(0.5f, 0.5f, 0.5f);
	glScalef(0.5f, 0.5f, 0.5f);
	glMultMatrixf(_projectionMatrix);
	glMultMatrixf(_viewMatrix);

	GLfloat vp_viewMatrix[4*4];
	for (size_t i(0); i < 4*4; ++i)
	  vp_viewMatrix[i] = vp.getViewMatrix()[i]; 

	invert(vp_viewMatrix);
	glMultMatrixf(vp_viewMatrix);
      }

      inline void invert(GLfloat* mat4x4)
      {
	GLfloat result[4*4];

	GLfloat tmp[12];											//temporary pair storage

	//calculate pairs for first 8 elements (cofactors)
	tmp[0] = mat4x4[10] * mat4x4[15];
	tmp[1] = mat4x4[11] * mat4x4[14];
	tmp[2] = mat4x4[9] * mat4x4[15];
	tmp[3] = mat4x4[11] * mat4x4[13];
	tmp[4] = mat4x4[9] * mat4x4[14];
	tmp[5] = mat4x4[10] * mat4x4[13];
	tmp[6] = mat4x4[8] * mat4x4[15];
	tmp[7] = mat4x4[11] * mat4x4[12];
	tmp[8] = mat4x4[8] * mat4x4[14];
	tmp[9] = mat4x4[10] * mat4x4[12];
	tmp[10] = mat4x4[8] * mat4x4[13];
	tmp[11] = mat4x4[9] * mat4x4[12];

	//calculate first 8 elements (cofactors)
	result[0] = tmp[0]*mat4x4[5] + tmp[3]*mat4x4[6] + tmp[4]*mat4x4[7]-tmp[1]*mat4x4[5] - tmp[2]*mat4x4[6] - tmp[5]*mat4x4[7];
	result[1] = tmp[1]*mat4x4[4] + tmp[6]*mat4x4[6] + tmp[9]*mat4x4[7]-tmp[0]*mat4x4[4] - tmp[7]*mat4x4[6] - tmp[8]*mat4x4[7];
	result[2] = tmp[2]*mat4x4[4] + tmp[7]*mat4x4[5] + tmp[10]*mat4x4[7]-tmp[3]*mat4x4[4] - tmp[6]*mat4x4[5] - tmp[11]*mat4x4[7];
	result[3] = tmp[5]*mat4x4[4] + tmp[8]*mat4x4[5] + tmp[11]*mat4x4[6]-tmp[4]*mat4x4[4] - tmp[9]*mat4x4[5] - tmp[10]*mat4x4[6];
	result[4] = tmp[1]*mat4x4[1] + tmp[2]*mat4x4[2] + tmp[5]*mat4x4[3]-tmp[0]*mat4x4[1] - tmp[3]*mat4x4[2] - tmp[4]*mat4x4[3];
	result[5] = tmp[0]*mat4x4[0] + tmp[7]*mat4x4[2] + tmp[8]*mat4x4[3]-tmp[1]*mat4x4[0] - tmp[6]*mat4x4[2] - tmp[9]*mat4x4[3];
	result[6] = tmp[3]*mat4x4[0] + tmp[6]*mat4x4[1] + tmp[11]*mat4x4[3]-tmp[2]*mat4x4[0] - tmp[7]*mat4x4[1] - tmp[10]*mat4x4[3];
	result[7] = tmp[4]*mat4x4[0] + tmp[9]*mat4x4[1] + tmp[10]*mat4x4[2]-tmp[5]*mat4x4[0] - tmp[8]*mat4x4[1] - tmp[11]*mat4x4[2];

	//calculate pairs for second 8 elements (cofactors)
	tmp[0] = mat4x4[2]*mat4x4[7];
	tmp[1] = mat4x4[3]*mat4x4[6];
	tmp[2] = mat4x4[1]*mat4x4[7];
	tmp[3] = mat4x4[3]*mat4x4[5];
	tmp[4] = mat4x4[1]*mat4x4[6];
	tmp[5] = mat4x4[2]*mat4x4[5];
	tmp[6] = mat4x4[0]*mat4x4[7];
	tmp[7] = mat4x4[3]*mat4x4[4];
	tmp[8] = mat4x4[0]*mat4x4[6];
	tmp[9] = mat4x4[2]*mat4x4[4];
	tmp[10] = mat4x4[0]*mat4x4[5];
	tmp[11] = mat4x4[1]*mat4x4[4];

	//calculate second 8 elements (cofactors)
	result[8 ] = tmp[0]*mat4x4[13] + tmp[3]*mat4x4[14] + tmp[4]*mat4x4[15]-tmp[1]*mat4x4[13] - tmp[2]*mat4x4[14] - tmp[5]*mat4x4[15];
	result[9 ] = tmp[1]*mat4x4[12] + tmp[6]*mat4x4[14] + tmp[9]*mat4x4[15]-tmp[0]*mat4x4[12] - tmp[7]*mat4x4[14] - tmp[8]*mat4x4[15];
	result[10] = tmp[2]*mat4x4[12] + tmp[7]*mat4x4[13] + tmp[10]*mat4x4[15]-tmp[3]*mat4x4[12] - tmp[6]*mat4x4[13] - tmp[11]*mat4x4[15];
	result[11] = tmp[5]*mat4x4[12] + tmp[8]*mat4x4[13] + tmp[11]*mat4x4[14]-tmp[4]*mat4x4[12] - tmp[9]*mat4x4[13] - tmp[10]*mat4x4[14];
	result[12] = tmp[2]*mat4x4[10] + tmp[5]*mat4x4[11] + tmp[1]*mat4x4[9]-tmp[4]*mat4x4[11] - tmp[0]*mat4x4[9] - tmp[3]*mat4x4[10];
	result[13] = tmp[8]*mat4x4[11] + tmp[0]*mat4x4[8] + tmp[7]*mat4x4[10]-tmp[6]*mat4x4[10] - tmp[9]*mat4x4[11] - tmp[1]*mat4x4[8];
	result[14] = tmp[6]*mat4x4[9] + tmp[11]*mat4x4[11] + tmp[3]*mat4x4[8]-tmp[10]*mat4x4[11] - tmp[2]*mat4x4[8] - tmp[7]*mat4x4[9];
	result[15] = tmp[10]*mat4x4[10] + tmp[4]*mat4x4[8] + tmp[9]*mat4x4[9]-tmp[8]*mat4x4[9] - tmp[11]*mat4x4[10] - tmp[5]*mat4x4[8];

	// calculate determinant
	GLfloat det = mat4x4[0]*result[0]+mat4x4[1]*result[1]+mat4x4[2]*result[2]
	  +mat4x4[3]*result[3];

	for (size_t i(0); i < 4*4; ++i)
	  result[i] /= det;

	if(det==0.0f)
	  for (size_t i(0); i < 4; ++i)
	    for (size_t j(0); j < 4; ++j)
	      result[4*i+j] = (i==j) ? 1 : 0;

	//Now we need to transpose the matrix as we copy out
	for (size_t i(0); i < 4; ++i)
	  for (size_t j(0); j < 4; ++j)
	    mat4x4[4*i+j] = result[4*j+i];
      }

      GLenum _lightHandle;
    };
  }
}
