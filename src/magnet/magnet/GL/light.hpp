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
      inline lightInfo() {}

      inline lightInfo(Vector position, 
		       Vector lookAtPoint,
		       GLenum lightHandle = GL_LIGHT0,
		       GLfloat fovY = 45.0f,
		       GLfloat zNearDist = 0.01f, GLfloat zFarDist = 10.0f,
		       Vector up = Vector(0,1,0)):
	viewPort(position, lookAtPoint, fovY, zNearDist, zFarDist, up, 1.0),
	_lightHandle(lightHandle)
      {
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
	GLdouble aspectRatio = _aspectRatio;
	double fovY = getFOVY();
	viewPort::operator=(vp);
	_aspectRatio = aspectRatio;
	setFOVY(fovY);
	
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

	MATRIX4X4 invView = vp.getViewMatrix();
	invView.Invert();
	glMultMatrixf(invView);
      }

      GLenum _lightHandle;
    };
  }
}
