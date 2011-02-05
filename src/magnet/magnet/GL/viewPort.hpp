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

#define GL_GLEXT_PROTOTYPES
#include <GL/glew.h>

#include <cmath>
#include <magnet/math/matrix.hpp>
#include <coil/Maths/VECTOR4D.h>
#include <coil/Maths/MATRIX4X4.h>

namespace magnet {
  namespace GL {    
    struct viewPort
    {
      inline viewPort(Vector position = Vector(1,1,1), 
		      Vector lookAtPoint = Vector(0,0,0),
		      GLfloat fovY = 45.0f,
		      GLfloat zNearDist = 0.01f, GLfloat zFarDist = 10.0f,
		      Vector up = Vector(0,1,0),
		      GLfloat aspectRatio = 1
		      ):
	_panrotation(180),
	_tiltrotation(0),
	_position(position),
	_fovY(fovY),
	_aspectRatio(aspectRatio),
	_zNearDist(zNearDist),
	_zFarDist(zFarDist)
      {
	if (_zNearDist > _zFarDist) 
	  M_throw() << "zNearDist > _zFarDist!";

	up /= up.nrm();
	
	//Now rotate about the up vector, we do tilt seperately
	Vector directionNorm = (lookAtPoint - position);
	directionNorm /= directionNorm.nrm();
	double upprojection = (directionNorm | up);
	Vector directionInXZplane = directionNorm - upprojection * up;
	directionInXZplane /= directionInXZplane.nrm();
	_panrotation = -(180.0 / M_PI) * std::acos(directionInXZplane | Vector(0,0,-1));
		
	Vector rotationAxis = up ^ directionInXZplane;
	rotationAxis /= rotationAxis.nrm();

	_tiltrotation = (180.0 / M_PI) * std::acos(directionInXZplane | directionNorm);
		
	buildMatrices();
      }
      
      inline void CameraUpdate(float forward = 0, float sideways = 0, float vertical = 0)
      {
	//Forward/Backward movement
	_position[2] -= forward * std::cos(_tiltrotation * (M_PI/ 180)) 
	  * std::sin(_panrotation  * (M_PI/ 180) + M_PI * 0.5);  
	_position[0] -= forward * std::cos(_tiltrotation * (M_PI/ 180)) 
	  * std::cos(_panrotation  * (M_PI/ 180) + M_PI * 0.5);
	_position[1] += -forward * std::sin(_tiltrotation * (M_PI/ 180));
	
	//Strafe movement
	_position[2] += sideways * std::sin(_panrotation * (M_PI/ 180));
	_position[0] += sideways * std::cos(_panrotation * (M_PI/ 180));
	
	//Vertical movement
	_position[1] += vertical;

	buildMatrices();
      }
      
      inline void buildMatrices()
      {
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	//Build the projection matrix
	glLoadIdentity();
	gluPerspective(_fovY, _aspectRatio, _zNearDist, _zFarDist);
	glGetFloatv(GL_MODELVIEW_MATRIX, _projectionMatrix);

	//setup the view matrix
	glLoadIdentity();
	glRotatef(_tiltrotation, 1.0, 0.0, 0.0);
	glRotatef(_panrotation, 0.0, 1.0, 0.0);
	glTranslatef(-_position[0], -_position[1], -_position[2]);	
	glGetFloatv(GL_MODELVIEW_MATRIX, _viewMatrix);
	
	glPopMatrix();

	Matrix viewTransform = Rodrigues(Vector(0,-_panrotation * M_PI/180,0)) 
	  * Rodrigues(Vector(-_tiltrotation * M_PI/180.0,0,0));
	
	_cameraDirection =  viewTransform * Vector(0,0,-1);
	_cameraUp = viewTransform * Vector(0,1,0);

      }

      inline void loadMatrices()
      {
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf(_projectionMatrix);
	
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixf(_viewMatrix);
      }

      float _panrotation;
      float _tiltrotation;
      Vector _position;
      
      GLdouble _fovY;
      GLdouble _aspectRatio;
      GLdouble _zNearDist;
      GLdouble _zFarDist;
      
      Vector _cameraDirection, _cameraUp;
      
      MATRIX4X4 _projectionMatrix;
      MATRIX4X4 _viewMatrix;
    };
  }
}
