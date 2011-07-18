/*    dynamo:- Event driven molecular dynamics simulator 
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

#include <magnet/GL/matrix.hpp>
#include <magnet/clamp.hpp>
#include <magnet/exception.hpp>

namespace magnet {
  namespace GL {
    /*! \brief An object to track the viewport (aka camera) state.
     *
     * This class can perform all the calculations required for
     * setting up the projection and modelview matricies of the
     * camera. There is also support for head tracking calculations
     * using the \ref _headLocation \ref Vector.
     */
    class ViewPort
    {
    public:
      //! \brief The mode of the mouse movement
      enum Camera_Mode
	{
	  ROTATE_VIEWPLANE,
	  ROTATE_CAMERA,
	  ROTATE_WORLD
	};

      /*! \brief The constructor
       * 
       * \param height The height of the viewport, in pixels.
       * \param width The width of the viewport, in pixels.
       * \param position The position of the screen (effectively the camera), in simulation coordinates.
       * \param lookAtPoint The location the camera is initially focussed on.
       * \param fovY The field of vision of the camera.
       * \param zNearDist The distance to the near clipping plane.
       * \param zFarDist The distance to the far clipping plane.
       * \param up A vector describing the up direction of the camera.
       */
      //We need a default constructor as viewPorts may be created without GL being initialized
      inline ViewPort(size_t height = 600, 
		      size_t width = 800,
		      Vector position = Vector(1,1,1), 
		      Vector lookAtPoint = Vector(0,0,0),
		      GLfloat fovY = 60.0f,
		      GLfloat zNearDist = 0.01f, GLfloat zFarDist = 20.0f,
		      Vector up = Vector(0,1,0)
		      ):
	_height(height),
	_width(width),
	_panrotation(180),
	_tiltrotation(0),
	_position(position),
	_zNearDist(zNearDist),
	_zFarDist(zFarDist),
	_headLocation(0, 0, 1),
	_simLength(25),
	_pixelPitch(0.025), //Measured from my screen
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

      /*! \brief Change the field of vision of the viewport/camera.
       *
       * \param fovY The field of vision in degrees.
       * \param compensate Counter the movement of the head position
       * by moving the viewing plane position.
       */
      inline void setFOVY(double fovY, bool compensate = true) 
      {
	//When the FOV is adjusted, we move the head position away
	//from the view plane, but we adjust the viewplane position to
	//compensate this motion
	Vector headLocationChange = Vector(0, 0, 0.5f * (_pixelPitch * _width / _simLength) 
					   / std::tan((fovY / 180.0f) * M_PI / 2) 
					   - _headLocation[2]);

	if (compensate)
	  {
	    ::Matrix viewTransformation 
	      = Rodrigues(Vector(0, -_panrotation * M_PI/180, 0))
	      * Rodrigues(Vector(-_tiltrotation * M_PI / 180.0, 0, 0));
	    
	    _position -= viewTransformation * headLocationChange;	
	  }

	_headLocation += headLocationChange;
      }
      
      /*! \brief Sets the OpenGL head location.
       *
       * \param head The position of the viewers head, relative to the
       * center of the near viewing plane (in cm).
       */
      inline void setHeadLocation(Vector head)
      { _headLocation = head / _simLength; }

      /*! \brief Gets the OpenGL head location (in cm).
       *
       * The position of the viewers head is relative to the center of
       * the near viewing plane (in cm).
       */
      inline const Vector getHeadLocation() const
      { return _headLocation * _simLength; }

      /*! \brief Returns the current field of vision of the viewport/camera */
      inline double getFOVY() const
      { return 2 * std::atan2(0.5f * (_pixelPitch * _width / _simLength),  _headLocation[2]) * (180.0f / M_PI); }

      /*! \brief Converts the motion of the mouse into a motion of the
       * viewport/camera.
       *
       * \param diffX The amount the mouse has moved in the x direction, in pixels.
       * \param diffY The amount the mouse has moved in the y direction, in pixels.
       */
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

      /*! \brief Converts a forward/sideways/vertical motion (e.g.,
       * obtained from keypresses) into a motion of the
       * viewport/camera.
       *
       * \param diffX The amount the mouse has moved in the x direction, in pixels.
       * \param diffY The amount the mouse has moved in the y direction, in pixels.
       */
      inline void CameraUpdate(float forward = 0, float sideways = 0, float vertical = 0)
      {
	//Build a matrix to rotate from camera to world
	::Matrix Transformation = Rodrigues(Vector(0,-_panrotation * M_PI/180,0))
	  * Rodrigues(Vector(-_tiltrotation * M_PI/180.0,0,0));
	
	//This vector is the movement vector from the camera's
	//viewpoint (not including the vertical component)
	Vector movement(sideways,0,-forward); //Strafe direction (left-right)
	
	_position += Transformation * movement + Vector(0,vertical,0);

	buildMatrices();
      }
      
      /*! \brief Constructs the OpenGL modelview and projection
       * matricies from the stored state of the viewport/camera.
       *
       * \param offset This is an offset in camera coordinates to
       * apply to the head location. It's primary use is to calculate
       * the perspective shift for the left and right eye in
       * Analygraph rendering.
       */
      inline void buildMatrices(Vector offset = Vector(0,0,0))
      {  
	//Local head location
	Vector headLoc = _headLocation + offset / _simLength;

	//Now add in the movement of the head and the movement of the
	//camera
	::Matrix viewTransformation 
	    = Rodrigues(Vector(0, -_panrotation * M_PI/180, 0))
	    * Rodrigues(Vector(-_tiltrotation * M_PI / 180.0, 0, 0));
	
	Vector cameraLocation((viewTransformation * headLoc) + _position);
	_cameraDirection = viewTransformation * Vector(0,0,-1);
	_cameraUp = viewTransformation * Vector(0,1,0);

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
	_projectionMatrix 
	  = GLMatrix::frustrum((-0.5f * getScreenPlaneWidth()  - headLoc[0]) * _zNearDist / headLoc[2],// left
			       (+0.5f * getScreenPlaneWidth()  - headLoc[0]) * _zNearDist / headLoc[2],// right
			       (-0.5f * getScreenPlaneHeight() - headLoc[1]) * _zNearDist / headLoc[2],// bottom 
			       (+0.5f * getScreenPlaneHeight() - headLoc[1]) * _zNearDist / headLoc[2],// top
			       _zNearDist,//Near distance
			       _zFarDist//Far distance)
			       );

	//setup the view matrix
	_viewMatrix
	  = GLMatrix::rotate(_tiltrotation, Vector(1,0,0))
	  * GLMatrix::rotate(_panrotation, Vector(0,1,0))
	  * GLMatrix::translate(-cameraLocation);
      }
      
      /*! \brief multiplies an inverse transformation of the
       * viewPort's modelview matrix with the current OpenGL matrix.
       *
       * This function simplifies drawing objects fixed in the camera
       * space (i.e. drawing a viewing frustrum for a given
       * viewport).
       *
       * \note This does not include any head tracking movement! (This is deliberate)
       */
      inline void applyInverseModelview() const
      {
	glTranslatef(_position[0], _position[1], _position[2]);
	glRotatef(-_panrotation, 0.0, 1.0, 0.0);
	glRotatef(-_tiltrotation, 1.0, 0.0, 0.0);
      }
      
      //! \brief Returns the screen's width (in simulation units).
      double getScreenPlaneWidth() const
      { return _pixelPitch * _width / _simLength; }

      //! \brief Returns the screen's height (in simulation units).
      double getScreenPlaneHeight() const
      { return _pixelPitch * _height / _simLength; }

      //! \brief Get the distance to the near clipping plane
      inline const GLfloat& getZNear() const { return _zNearDist; }
      //! \brief Get the distance to the far clipping plane
      inline const GLfloat& getZFar() const { return _zFarDist; }

      //! \brief Get the pan angle of the camera in degrees
      inline const float& getPan() const { return _panrotation; }

      //! \brief Get the tilt angle of the camera in degrees
      inline const float& getTilt() const { return _tiltrotation; }

      //! \brief Get the position of the viewing plane (effectively the camera position)
      inline const Vector& getViewPlanePosition() const { return _position; } 

      /*! \brief Get the modelview matrix.
       */
      inline const GLMatrix getViewMatrix() const { return _viewMatrix; }

      /*! \brief Get the projection matrix.
       */
      inline const GLMatrix getProjectionMatrix() const { return _projectionMatrix; }

      /*! \brief Fetch the location of the users eyes, in simulation
       * coordinates.
       * 
       * Useful for head tracking applications. This returns the
       * position of the eyes in simulation space by adding the head
       * location (relative to the viewing plane/screen) onto the
       * current position.
       */
      inline const Vector 
      getEyeLocation() const 
      { 
	::Matrix viewTransformation 
	  = Rodrigues(Vector(0, -_panrotation * M_PI/180, 0))
	  * Rodrigues(Vector(-_tiltrotation * M_PI / 180.0, 0, 0));

	return (viewTransformation * _headLocation) + _position;
      }

      //! \brief Set the height and width of the screen in pixels.
      inline void setHeightWidth(size_t height, size_t width)
      { _height = height; _width = width; }

      //! \brief Get the aspect ratio of the screen
      inline GLfloat getAspectRatio() const 
      { return ((GLfloat)_width) / _height; }

      //! \brief Get the up direction of the camera/viewport
      inline const Vector& getCameraUp() const { return _cameraUp; } 

      //! \brief Get the direction the camera is pointing in
      inline const Vector& getCameraDirection() const { return _cameraDirection; }

      //! \brief Get the height of the screen, in pixels.
      inline const size_t& getHeight() const { return _height; }

      //! \brief Get the width of the screen, in pixels.
      inline const size_t& getWidth() const { return _width; }
      
      /*! \brief Gets the simulation unit length (in cm). */
      inline const double& getSimUnitLength() const { return _simLength; }
      /*! \brief Sets the simulation unit length (in cm). */
      inline void setSimUnitLength(double val)  { _simLength = val; }

      /*! \brief Gets the pixel "diameter" in cm. */
      inline const double& getPixelPitch() const { return _pixelPitch; }
      /*! \brief Sets the pixel "diameter" in cm. */
      inline void setPixelPitch(double val)  { _pixelPitch = val; }

    protected:
      size_t _height, _width;
      float _panrotation;
      float _tiltrotation;
      Vector _position;
      
      GLfloat _zNearDist;
      GLfloat _zFarDist;
      Vector _headLocation;
      Vector _cameraDirection, _cameraUp;
      
      GLMatrix _projectionMatrix;
      GLMatrix _viewMatrix;

      //! \brief One simulation length in cm (real units)
      double _simLength;

      //! \brief The diameter of a pixel, in cm
      double _pixelPitch;

      Camera_Mode _camMode;
    };
  }
}
