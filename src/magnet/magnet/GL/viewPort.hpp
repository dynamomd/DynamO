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

#include <magnet/math/matrix.hpp>
#include <magnet/clamp.hpp>

#include <magnet/exception.hpp>
#include <cmath>

namespace magnet {
  namespace GL {
    class viewPort
    {
    public:
      enum Camera_Mode
	{
	  ROTATE_VIEWPLANE,
	  ROTATE_CAMERA,
	  ROTATE_WORLD
	};

      //We need a default constructor as viewPorts may be created without GL being initialized
      inline viewPort(size_t height = 600, 
		      size_t width = 800,
		      Vector position = Vector(1,1,1), 
		      Vector lookAtPoint = Vector(0,0,0),
		      GLfloat fovY = 60.0f,
		      GLfloat zNearDist = 0.05f, GLfloat zFarDist = 20.0f,
		      Vector up = Vector(0,1,0),
		      GLfloat aspectRatio = 1
		      ):
	_height(height),
	_width(width),
	_panrotation(180),
	_tiltrotation(0),
	_position(position),
	_zNearDist(zNearDist),
	_zFarDist(zFarDist),
	_simLength(50),
	_screenWidth(41.1f / _simLength),
	_camMode(ROTATE_CAMERA)
      {
	if (_zNearDist > _zFarDist) 
	  M_throw() << "zNearDist > _zFarDist!";

	up /= up.nrm();
	
	//Now rotate about the up vector, we do tilt seperately
	Vector directionNorm = (lookAtPoint - position);
	directionNorm /= directionNorm.nrm();
	double upprojection = (directionNorm | up);
	Vector directionInXZplane = directionNorm - upprojection * up;
	directionInXZplane /= (directionInXZplane.nrm() != 0) ? directionInXZplane.nrm() : 0;
	_panrotation = -(180.0f / M_PI) * std::acos(directionInXZplane | Vector(0,0,-1));
		
	Vector rotationAxis = up ^ directionInXZplane;
	rotationAxis /= rotationAxis.nrm();

	_tiltrotation = (180.0f / M_PI) * std::acos(directionInXZplane | directionNorm);

	//We use the field of vision and the width of the screen in
	//simulation units to calculate how far back the head should
	//be at the start
	setFOVY(fovY);
      }

      inline void setFOVY(double fovY) 
      {
	//When the FOV is adjusted, we move the head position away
	//from the view plane, but we adjust the viewplane position to
	//compensate this motion
	Vector headLocationChange = Vector(0, 0, 0.5f * _screenWidth 
					   / std::tan((fovY / 180.0f) * M_PI / 2) 
					   - _headLocation[2]);

	Matrix viewTransformation 
	  = Rodrigues(Vector(0, -_panrotation * M_PI/180, 0))
	  * Rodrigues(Vector(-_tiltrotation * M_PI / 180.0, 0, 0));

	_position -= viewTransformation * headLocationChange;	
	_headLocation += headLocationChange;
      }
      
      inline double getFOVY() const
      { return 2 * std::atan2(0.5f * _screenWidth,  _headLocation[2]) * (180.0f / M_PI); }

      inline void mouseMovement(float diffX, float diffY)
      {
	switch (_camMode)
	  {
	  case ROTATE_VIEWPLANE:	
	    { //The camera is rotated and appears to rotate around the
	      //view plane
	      _panrotation += diffX;
	      _tiltrotation = magnet::clamp(diffY + _tiltrotation, -90.0f, 90.0f);
	      break;
	    }
	  case ROTATE_CAMERA:
	    { //The camera is rotated and an additional movement is
	      //added to make it appear to rotate around the head
	      //position
	      //Calculate the current camera position
	      Vector cameraLocationOld(getEyeLocation());	      
	      _panrotation += diffX;
	      _tiltrotation = magnet::clamp(diffY + _tiltrotation, -90.0f, 90.0f);	      
	      Vector cameraLocationNew(getEyeLocation());

	      _position -= cameraLocationNew - cameraLocationOld;	      
	      break;
	    }
	  case ROTATE_WORLD:
	  default:
	    M_throw() << "Bad camera mode";
	  }
      }

      inline void CameraUpdate(float forward = 0, float sideways = 0, float vertical = 0)
      {
	//Build a matrix to rotate from camera to world
	Matrix Transformation = Rodrigues(Vector(0,-_panrotation * M_PI/180,0))
	  * Rodrigues(Vector(-_tiltrotation * M_PI/180.0,0,0));
	
	//This vector is the movement vector from the camera's
	//viewpoint (not including the vertical component)
	Vector movement(sideways,0,-forward); //Strafe direction (left-right)
	
	_position += Transformation * movement + Vector(0,vertical,0);

	buildMatrices();
      }
      
      inline void buildMatrices()
      {  
	//We'll build the matricies on the modelview stack
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	//Build the projection matrix
	glLoadIdentity();

	//We will move the camera to the location of the head in sim
	//space. So we must create a viewing frustrum which, in real
	//space, cuts through the image on the screen. The trick is to
	//take the real world relative coordinates of the screen and
	//head transform them to simulation units.
	//
	//This allows us to calculate the left, right, bottom and top of
	//the frustrum as if the near plane of the frustrum was at the
	//screens location.
	//
	//Finally, all length scales are multiplied by
	//(_zNearDist/_headLocation[2]).
	//
	//This is to allow the frustrum's near plane to be placed
	//somewhere other than the screen (this factor places it at
	//_zNearDist)!
	//
	glFrustum((-0.5f * _screenWidth                - _headLocation[0]) * _zNearDist / _headLocation[2],// left
		  (+0.5f * _screenWidth                - _headLocation[0]) * _zNearDist / _headLocation[2],// right
		  (-0.5f * _screenWidth / getAspectRatio() - _headLocation[1]) * _zNearDist / _headLocation[2],// bottom 
		  (+0.5f * _screenWidth / getAspectRatio() - _headLocation[1]) * _zNearDist / _headLocation[2],// top
		  _zNearDist,//Near distance
		  _zFarDist//Far distance
		  );

	glGetFloatv(GL_MODELVIEW_MATRIX, _projectionMatrix);

	//setup the view matrix
	glLoadIdentity();
	glRotatef(_tiltrotation, 1.0, 0.0, 0.0);
	glRotatef(_panrotation, 0.0, 1.0, 0.0);

	//Now add in the movement of the head and the movement of the
	//camera
	Matrix viewTransformation 
	  = Rodrigues(Vector(0, -_panrotation * M_PI/180, 0))
	  * Rodrigues(Vector(-_tiltrotation * M_PI / 180.0, 0, 0));

	Vector cameraLocation((viewTransformation * _headLocation) + _position);

	glTranslatef(-cameraLocation[0], -cameraLocation[1], -cameraLocation[2]);

	glGetFloatv(GL_MODELVIEW_MATRIX, _viewMatrix);
	glPopMatrix();
	
	_cameraDirection = viewTransformation * Vector(0,0,-1);
	_cameraUp = viewTransformation * Vector(0,1,0);

      }

      inline void saveMatrices()
      {
	glGetFloatv(GL_PROJECTION_MATRIX, _projectionMatrix);
	glGetFloatv(GL_MODELVIEW_MATRIX, _viewMatrix);
      }

      inline void loadMatrices()
      {
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf(_projectionMatrix);
	
	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixf(_viewMatrix);
      }

      inline const GLdouble& getZNear() const { return _zNearDist; }
      inline const GLdouble& getZFar() const { return _zFarDist; }
      inline const float& getPan() const { return _panrotation; }
      inline const float& getTilt() const { return _tiltrotation; }
      inline const Vector& getViewPlanePosition() const { return _position; } 
      inline const GLfloat* getViewMatrix() const { return _viewMatrix; }

      inline const Vector 
      getEyeLocation() const 
      { 
	Matrix viewTransformation 
	  = Rodrigues(Vector(0, -_panrotation * M_PI/180, 0))
	  * Rodrigues(Vector(-_tiltrotation * M_PI / 180.0, 0, 0));

	return (viewTransformation * _headLocation) + _position;
      }

      inline void setHeightWidth(size_t height, size_t width)
      { _height = height; _width = width; }

      inline GLdouble getAspectRatio() const 
      { return ((GLdouble)_width) / _height; }

      inline const Vector& getCameraUp() const { return _cameraUp; } 
      inline const Vector& getCameraDirection() const { return _cameraDirection; }

      inline const size_t& getHeight() const { return _height; }
      inline const size_t& getWidth() const { return _width; }

    protected:
      size_t _height, _width;
      float _panrotation;
      float _tiltrotation;
      Vector _position;
      
      GLdouble _zNearDist;
      GLdouble _zFarDist;
      Vector _headLocation;
      Vector _cameraDirection, _cameraUp;
      
      GLfloat _projectionMatrix[4*4];
      GLfloat _viewMatrix[4*4];

      //!One simulation length in cm (real units)
      double _simLength;
      //!The width of the screen in simulation units
      double _screenWidth;

      Camera_Mode _camMode;
    };
  }
}
