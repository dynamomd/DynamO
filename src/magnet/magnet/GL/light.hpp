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

namespace magnet {
  namespace GL {        
    class lightInfo
    {
    public:
      inline lightInfo() {}
      
      inline lightInfo(Vector position, Vector lookAtPoint, GLenum lightHandle, 
		GLfloat beamAngle = 45.0f,
		GLfloat rangeMax = 3.0f, GLfloat rangeMin = 0.01):
	_position(position),
	_lookAtPoint(lookAtPoint),
	_lightHandle(lightHandle)
      {
	//Setup a bright light
	glLightfv(lightHandle, GL_DIFFUSE, white);
	glLightfv(lightHandle, GL_SPECULAR, white);

	//Build the view matrix and so on
	glPushMatrix();
	
	glLoadIdentity();
	gluPerspective(beamAngle, 1.0f, rangeMin, rangeMax);
	glGetFloatv(GL_MODELVIEW_MATRIX, _projectionMatrix);
	
	glLoadIdentity();
	Vector directionNorm = (lookAtPoint - position);
	directionNorm /= directionNorm.nrm();
	
	GLfloat rotationAngle = (180.0 / M_PI) * std::acos(Vector(0,0,-1) | directionNorm);
	
	
	Vector RotationAxis = Vector(0,0,-1) ^ directionNorm;
	float norm = RotationAxis.nrm();
	RotationAxis /= norm;
	if (norm < std::numeric_limits<double>::epsilon())
	  RotationAxis = Vector(1,0,0);
	
	glRotatef(-rotationAngle, RotationAxis.x,RotationAxis.y,RotationAxis.z);
	glTranslatef(-position.x,-position.y,-position.z);
	glGetFloatv(GL_MODELVIEW_MATRIX, _viewMatrix);
	
	glPopMatrix();
      }

      inline void drawLight()
      {
	glColor3f(1,1,0);
	
	Vector directionNorm = (_lookAtPoint - _position);
	directionNorm /= directionNorm.nrm();
    
	GLfloat rotationAngle = (180.0 / M_PI) * std::acos(Vector(0,0,-1) | directionNorm);
	
	
	Vector RotationAxis = Vector(0,0,-1) ^ directionNorm;
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

      inline void buildShadowTextureMatrix()
      {
	//Build the texture matrix
	glLoadIdentity();
	glTranslatef(0.5f, 0.5f, 0.5f);
	glScalef(0.5f, 0.5f, 0.5f);
	glMultMatrixf(_projectionMatrix);
	glMultMatrixf(_viewMatrix);
      }

      Vector _position;
      Vector _lookAtPoint;
      
      MATRIX4X4 _projectionMatrix;
      MATRIX4X4 _viewMatrix;
      
      GLenum _lightHandle;
    };
  }
}
