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

namespace magnet {
  namespace GL {    
    struct viewPort
    {
      inline viewPort():
	_rotatex(180),
	_rotatey(0),
	_cameraX(0),
	_cameraY(0),
	_cameraZ(0),
	_fovY(45),
	_aspectRatio(1),
	_zNearDist(0.001),
	_zFarDist(100)
      {}
      
      inline void CameraSetup(float forward = 0, float sideways = 0, float vertical = 0)
      {
	//Forward/Backward movement
	_cameraZ -= forward * std::cos(_rotatey * (M_PI/ 180)) 
	  * std::sin(_rotatex  * (M_PI/ 180) + M_PI * 0.5);  
	_cameraX -= forward * std::cos(_rotatey * (M_PI/ 180)) 
	  * std::cos(_rotatex  * (M_PI/ 180) + M_PI * 0.5);
	_cameraY += -forward * std::sin(_rotatey * (M_PI/ 180));
	
	//Strafe movement
	_cameraZ += sideways * std::sin(_rotatex * (M_PI/ 180));
	_cameraX += sideways * std::cos(_rotatex * (M_PI/ 180));
	
	//Vertical movement
	_cameraY += vertical;
	
	glLoadIdentity();
	//gluLookAt(-viewscale, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
	glRotatef(_rotatey, 1.0, 0.0, 0.0);
	glRotatef(_rotatex, 0.0, 1.0, 0.0);
	glTranslatef(-_cameraX,-_cameraY,-_cameraZ);
	
	//store the matricies for shadow calculations
	glGetFloatv(GL_MODELVIEW_MATRIX, _viewMatrix);
	glGetFloatv(GL_PROJECTION_MATRIX, _projectionMatrix);
	
	Matrix viewTransform = Rodrigues(Vector(0,-_rotatex * M_PI/180,0)) 
	  * Rodrigues(Vector(-_rotatey * M_PI/180.0,0,0));
	
	_cameraDirection =  viewTransform * Vector(0,0,-1);
	_cameraUp = viewTransform * Vector(0,1,0);
      }
      
      float _rotatex;
      float _rotatey;
      float _cameraX;
      float _cameraY;
      float _cameraZ;
      
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
