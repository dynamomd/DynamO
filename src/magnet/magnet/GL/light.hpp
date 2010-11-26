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
      inline lightInfo(GLenum lightHandle = GL_LIGHT0, 
		       Vector position = Vector(0,2,0), 
		       Vector lookAtPoint = Vector(0,0,0),
		       GLfloat fovY = 45.0f,
		       GLfloat zNearDist = 0.001f, GLfloat zFarDist = 100.0f,
		       Vector up = Vector(0,1,0)):
	viewPort(position, lookAtPoint, fovY, zNearDist, zFarDist, up, 1.0),
	_lightHandle(lightHandle)
      {
	//Setup a bright light
	glLightfv(lightHandle, GL_DIFFUSE, white);
	glLightfv(lightHandle, GL_SPECULAR, white);
      }

      inline void drawLight()
      {
	glColor3f(1,1,0);
	
	GLfloat rotationAngle 
	  = (180.0 / M_PI) * std::acos(Vector(0,0,-1) | _cameraDirection);
	
	Vector RotationAxis = Vector(0,0,-1) ^ _cameraDirection;
	float norm = RotationAxis.nrm();
	RotationAxis /= norm;
	if (norm < std::numeric_limits<double>::epsilon())
	  RotationAxis = Vector(1,0,0);
	
	glPushMatrix();
	glTranslatef(_position.x, _position.y, _position.z);
	glRotatef(rotationAngle, RotationAxis.x, RotationAxis.y, RotationAxis.z);
	
	glTranslatef(0.0f,0.0f,-0.025f);
	glutSolidTorus(0.01f, 0.02f, 20, 20);
	glTranslatef(0.0f,0.0f,-0.025f);
	glutSolidCone(0.02f, 0.05f, 15, 15);
	glPopMatrix();
      }

      inline void loadShadowTextureMatrix(const viewPort& vp)
      {
	//Build the texture matrix
	glLoadIdentity();
	glTranslatef(0.5f, 0.5f, 0.5f);
	glScalef(0.5f, 0.5f, 0.5f);
	glMultMatrixf(_projectionMatrix);
	glMultMatrixf(_viewMatrix);
	MATRIX4X4 invView = vp._viewMatrix.GetInverse();
	glMultMatrixf(invView);
      }

      GLenum _lightHandle;
    };
  }
}
