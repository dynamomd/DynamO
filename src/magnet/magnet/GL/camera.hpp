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

#include <magnet/GL/context.hpp>
#include <magnet/GL/matrix.hpp>
#include <magnet/clamp.hpp>
#include <magnet/exception.hpp>
#include <iostream>

namespace magnet {
  namespace GL {
    /*! \brief An object to track the camera state.
     *
     * This class can perform all the calculations required for
     * setting up the projection and modelview matricies of the
     * camera. There is also support for eye tracking calculations
     * using the \ref _eyeLocation \ref math::Vector.
     */
    class Camera
    {
    public:
      //! \brief The mode of the mouse movement
      enum Camera_Mode
	{
	  ROTATE_CAMERA,
	  ROTATE_WORLD,
	  ROTATE_POINT
	};

      /*! \brief The constructor.
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
      inline Camera(size_t height = 600, 
		    size_t width = 800,
		    math::Vector position = math::Vector(0,0,-5), 
		    math::Vector lookAtPoint = math::Vector(0,0,1),
		    GLfloat fovY = 60.0f,
		    GLfloat zNearDist = 0.001f, GLfloat zFarDist = 100.0f,
		    math::Vector up = math::Vector(0,1,0)
		    ):
	_height(height),
	_width(width),
	_panrotation(180),
	_tiltrotation(0),
	_nearPlanePosition(0,0,0),
	_up(up),
	_rotatePoint(0,0,0),
	_zNearDist(zNearDist),
	_zFarDist(zFarDist),
	_simLength(25),
	_pixelPitch(0.05), //Measured from my screen
	_camMode(ROTATE_WORLD)
      {
	up /= up.nrm();

	if (_zNearDist > _zFarDist) 
	  M_throw() << "zNearDist > _zFarDist!";


	//We assume the user is around about 70cm from the screen
	setEyeLocation(math::Vector(0, 0, 70));
	setEyeLocationObjSpace(position);
	//lookAt(lookAtPoint);
      }

      inline void lookAt(math::Vector lookAtPoint)
      {
	//Generate the direction from the near plane to the object
	math::Vector directionNorm = (lookAtPoint - _nearPlanePosition);
	directionNorm /= directionNorm.nrm();

	//Now generate the direction from the eye to the near plane
	math::Matrix viewTransformation 
	    = Rodrigues(- _up * (_panrotation * M_PI/180))
	    * Rodrigues(math::Vector(-_tiltrotation * M_PI / 180.0, 0, 0));
	math::Vector eyeToPlane = viewTransformation * _eyeLocation;

	//Now generate the eye to the object vector
	math::Vector eyeToObject = (lookAtPoint - (eyeToPlane + _nearPlanePosition));
	

	//If the object lies between the eye and the plane, flip the direction
	if (((eyeToPlane | directionNorm) >= 0) && ((eyeToObject | directionNorm) <= 0))
	  directionNorm = -directionNorm;

	double upprojection = (directionNorm | _up);

	if (upprojection == 1)
	  {
	    _tiltrotation = -90;
	    _panrotation = 0;
	    return;
	  }
	else if (upprojection == -1)
	  {
	    _tiltrotation = 90;
	    _panrotation = 0;	    
	    return;
	  }
	math::Vector directionInXZplane = directionNorm - upprojection * _up;
	directionInXZplane /= (directionInXZplane.nrm() != 0) ? directionInXZplane.nrm() : 0;

	math::Vector rotationAxis = _up ^ directionInXZplane;
	rotationAxis /= rotationAxis.nrm();

	_tiltrotation = (180.0f / M_PI) * std::acos(std::min(directionInXZplane | directionNorm, 1.0));

	if (((directionNorm ^ directionInXZplane) | rotationAxis) > 0)
	  _tiltrotation = -_tiltrotation;

	_panrotation = -(180.0f / M_PI) * std::acos(std::min(directionInXZplane | math::Vector(0,0,-1), 1.0));
	if (((math::Vector(0,0,-1) ^ directionInXZplane) | _up) < 0)
	  _panrotation = -_panrotation;

#ifdef MAGNET_DEBUG
	if (std::isnan(_tiltrotation)) M_throw() << "Tilt Rotation is Nan\n" 
						 << "Tilt acos arg is " <<  (directionInXZplane | directionNorm);
	if (std::isnan(_panrotation)) M_throw() << "Pan Rotation is Nan";
#endif
      }

      inline void setRotatePoint(math::Vector vec)
      {
	if (_camMode == ROTATE_POINT)
	  _nearPlanePosition += vec - _rotatePoint;

	_rotatePoint = vec;
      }

      /*! \brief Change the field of vision of the camera.
       
        \param fovY The field of vision in degrees.
        \param compensate Counter the movement of the eye position
        by moving the viewing plane position.
       */
      inline void setFOVY(double fovY, bool compensate = true) 
      {
	//When the FOV is adjusted, we move the eye position away
	//from the view plane, but we adjust the viewplane position to
	//compensate this motion
	math::Vector eyeLocationChange = math::Vector(0, 0, 0.5f * (_pixelPitch * _width / _simLength) 
					   / std::tan((fovY / 180.0f) * M_PI / 2) 
					   - _eyeLocation[2]);

	if (compensate)
	  {
	    math::Matrix viewTransformation 
	      = Rodrigues(-_up * (_panrotation * M_PI/180))
	      * Rodrigues(math::Vector(-_tiltrotation * M_PI / 180.0, 0, 0));
	    
	    _nearPlanePosition -= viewTransformation * eyeLocationChange;	
	  }

	_eyeLocation += eyeLocationChange;
      }
      
      /*! \brief Sets the eye location.
       
        \param eye The position of the viewers eye, relative to the
        center of the near viewing plane (in cm).
       */
      inline void setEyeLocation(math::Vector eye)
      { _eyeLocation = eye / _simLength; }

      /*! \brief Gets the eye location (in cm).
       
        The position of the viewers eye is relative to the center of
        the near viewing plane (in cm).
       */
      inline const math::Vector getEyeLocation() const
      { return _eyeLocation * _simLength; }

      /*! \brief Returns the current field of vision of the camera */
      inline double getFOVY() const
      { return 2 * std::atan2(0.5f * (_pixelPitch * _width / _simLength),  _eyeLocation[2]) * (180.0f / M_PI); }

      /*! \brief Converts some inputted motion (e.g., by the mouse or keyboard) into a
        motion of the camera.


	All parameters may be negative or positive, as the sign
	defines the direction of the rotation/movement. Their name
	hints at what action they may do, depending on the camera mode
	(\ref _camMode).
       */
      inline void movement(float rotationX, float rotationY, float forwards, float sideways, float upwards)
      {
	//Build a matrix to rotate from camera to world
	math::Matrix Transformation = Rodrigues(-_up * (_panrotation * M_PI / 180.0))
	  * Rodrigues(math::Vector(- _tiltrotation * M_PI / 180.0, 0, 0));
	
	if ((_camMode == ROTATE_POINT) || (_camMode == ROTATE_WORLD))
	  {
	    if (forwards)
	      {
		//Test if the forward motion will take the eyePosition passed the viewing point, if so, don't move.
		math::Vector focus = math::Vector(0,0,0);
		if (_camMode == ROTATE_POINT)
		  focus = _rotatePoint;
		
		if (math::Vector(getEyeLocationObjSpace() - focus).nrm() > forwards)
		  _nearPlanePosition += Transformation * math::Vector(0, 0, -forwards);
	      }
	    
	    rotationX -= 10 * sideways;
	    rotationY += 10 * upwards;
	  }

	math::Vector focus = math::Vector(0,0,0);
	switch (_camMode)
	  {
	  case ROTATE_CAMERA:
	    { 
	      if (sideways || forwards)
		_nearPlanePosition += Transformation * math::Vector(sideways, 0, -forwards);
	      
	      _nearPlanePosition += math::Vector(0, upwards, 0);

	      //The camera is rotated and an additional movement is
	      //added to make it appear to rotate around the eye
	      //position
	      //Calculate the current camera position
	      math::Vector cameraLocationOld(getEyeLocationObjSpace());
	      _panrotation += rotationX;
	      _tiltrotation = magnet::clamp(rotationY + _tiltrotation, -90.0f, 90.0f);
	      math::Vector cameraLocationNew(getEyeLocationObjSpace());

	      _nearPlanePosition -= cameraLocationNew - cameraLocationOld;
	      break;
	    }
	  case ROTATE_POINT:
	    focus = _rotatePoint;
	  case ROTATE_WORLD:
	    {
	      lookAt(focus);

	      math::Vector offset = _nearPlanePosition - focus;

	      //We need to store the normal and restore it later.
	      double offset_length = offset.nrm();

	      if (rotationX)
		offset = Rodrigues(- _up * (M_PI * rotationX / 180.0f)) * offset;

	      if (rotationY)
		{

		  //We prevent the following command from returning
		  //bad axis by always using lookAt first, then
		  //getCameraUp().
		  math::Vector rotationAxis =  offset ^ getCameraUp();
		  double norm = rotationAxis.nrm();

#ifdef MAGNET_DEBUG
		  if (norm == 0)
		    M_throw() << "Bad normal on a camera rotation axis";
#endif

		  rotationAxis /= norm;
		  offset = Rodrigues(M_PI * (rotationY / 180.0f) * rotationAxis) * offset;
		}

	      offset *= offset_length / double(offset.nrm());

	      _nearPlanePosition = offset + focus;

	      lookAt(focus);
	      break;
	    }
	  default:
	    M_throw() << "Bad camera mode";
	  }
      }

      /*! \brief Tell the camera to align its view along an axis.
	
	This is useful when you want to reset the view
       */
      inline void setViewAxis(math::Vector axis)
      {
	switch (_camMode)
	  {
	  case ROTATE_CAMERA:
	    { 
	      lookAt(_nearPlanePosition + axis);
	      break;
	    }
	  case ROTATE_WORLD:
	    {
	      double origin_distance = _nearPlanePosition.nrm();
	      
	      _nearPlanePosition = - axis * origin_distance;
	      lookAt(math::Vector(0,0,0));
	      break;
	    }
	  default:
	    M_throw() << "Bad camera mode";
	  }

      }

      /*! \brief Get the modelview matrix. */
      inline GLMatrix getViewMatrix() const 
      {
	//Add in the movement of the eye and the movement of the
	//camera
	math::Matrix viewTransformation 
	  = Rodrigues(- _up * (_panrotation * M_PI/180))
	  * Rodrigues(math::Vector(-_tiltrotation * M_PI / 180.0, 0, 0));
	
	math::Vector cameraLocation((viewTransformation * _eyeLocation) + _nearPlanePosition);

	//Setup the view matrix
	return getViewRotationMatrix()
	  * GLMatrix::translate(-cameraLocation);
      }

      /*! \brief Generate a matrix that locates objects at the near
          ViewPlane (for rendering 3D objects attached to the
          screen). 
      */
      inline GLMatrix getViewPlaneMatrix() const
      {
	return getViewMatrix()
	  * GLMatrix::translate(_nearPlanePosition)
	  * GLMatrix::rotate(-_panrotation, _up)
	  * GLMatrix::rotate(-_tiltrotation, math::Vector(1, 0, 0));
      }

      /*! \brief Get the rotation part of the getViewMatrix().
       */
      inline GLMatrix getViewRotationMatrix() const 
      { 
	return GLMatrix::rotate(_tiltrotation, math::Vector(1,0,0))
	  * GLMatrix::rotate(_panrotation, _up);
      }

      /*! \brief Get the projection matrix.
       
        \param offset This is an offset in camera coordinates to apply
        to the eye location. It's primary use is to calculate the
        perspective shift for the left and right eye in Analygraph
        rendering.
       
	\param zoffset The amount to bias the depth values in the
	camera. See \ref GLMatrix::frustrum() for more information as
	the parameter is directly passed to that function.
       */
      inline GLMatrix getProjectionMatrix(GLfloat zoffset = 0) const 
      { 
	//We will move the camera to the location of the eye in sim
	//space. So we must create a viewing frustrum which, in real
	//space, cuts through the image on the screen. The trick is to
	//take the real world relative coordinates of the screen and
	//eye transform them to simulation units.
	//
	//This allows us to calculate the left, right, bottom and top of
	//the frustrum as if the near plane of the frustrum was at the
	//screens location.
	//
	//Finally, all length scales are multiplied by
	//(_zNearDist/_eyeLocation[2]).
	//
	//This is to allow the frustrum's near plane to be placed
	//somewhere other than the screen (this factor places it at
	//_zNearDist)!
	//
	return GLMatrix::frustrum((-0.5f * getScreenPlaneWidth()  - _eyeLocation[0]) * _zNearDist / _eyeLocation[2],// left
				  (+0.5f * getScreenPlaneWidth()  - _eyeLocation[0]) * _zNearDist / _eyeLocation[2],// right
				  (-0.5f * getScreenPlaneHeight() - _eyeLocation[1]) * _zNearDist / _eyeLocation[2],// bottom 
				  (+0.5f * getScreenPlaneHeight() - _eyeLocation[1]) * _zNearDist / _eyeLocation[2],// top
				  _zNearDist,//Near distance
				  _zFarDist,//Far distance
				  zoffset
				  );
      }
      
      /*! \brief Get the normal matrix.
       
        \param offset This is an offset in camera coordinates to apply
        to the eye location. It's primary use is to calculate the
        perspective shift for the left and right eye in Analygraph
        rendering.
       */
      inline math::Matrix getNormalMatrix() const 
      { return Inverse(math::Matrix(getViewMatrix())); }

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

      /*! \brief Fetch the location of the users eyes, in object space
        coordinates.
        
        Useful for eye tracking applications. This returns the
        position of the eyes in object space by adding the eye
        location (relative to the viewing plane/screen) onto the
        current position.
       */
      inline math::Vector 
      getEyeLocationObjSpace() const 
      { 
	math::Matrix viewTransformation 
	  = Rodrigues(- _up * (_panrotation * M_PI/180))
	  * Rodrigues(math::Vector(-_tiltrotation * M_PI / 180.0, 0, 0));

	return (viewTransformation * _eyeLocation) + _nearPlanePosition;
      }

      inline void
      setEyeLocationObjSpace(math::Vector Vec)
      {
	math::Matrix viewTransformation
	  = Rodrigues(- _up * (_panrotation * M_PI/180))
	  * Rodrigues(math::Vector(-_tiltrotation * M_PI / 180.0, 0, 0));
	
	_nearPlanePosition = Vec - (viewTransformation * _eyeLocation);
      }

      //! \brief Set the height and width of the screen in pixels.
      inline void setHeightWidth(size_t height, size_t width)
      { _height = height; _width = width; }

      //! \brief Get the aspect ratio of the screen
      inline GLfloat getAspectRatio() const 
      { return ((GLfloat)_width) / _height; }

      //! \brief Get the up direction of the camera.
      inline math::Vector getCameraUp() const 
      { 
	math::Matrix viewTransformation 
	  = Rodrigues(- _up * (_panrotation * M_PI/180))
	  * Rodrigues(math::Vector(-_tiltrotation * M_PI / 180.0, 0, 0));
	return viewTransformation * math::Vector(0,1,0);
      } 

      //! \brief Get the direction the camera is pointing in
      inline math::Vector getCameraDirection() const
      { 
	math::Matrix viewTransformation 
	  = Rodrigues(- _up * (_panrotation * M_PI/180))
	  * Rodrigues(math::Vector(-_tiltrotation * M_PI / 180.0, 0, 0));
	return viewTransformation * math::Vector(0,0,-1);
      } 

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

      Camera_Mode getMode() const { return _camMode; }

      void setMode(Camera_Mode val) { _camMode = val; }

      /*! \brief Used to convert mouse positions (including depth
          information) into a 3D position.
       */
      math::Vector unprojectToPosition(int windowx, int windowy, GLfloat depth)
      {
	//We need to calculate the ray from the camera
	std::tr1::array<GLfloat, 4> n = {{(2.0 * windowx) / getWidth() - 1.0,
					  1.0 - (2.0 * windowy) / getHeight(),
					  depth, 1.0}};
	//Unproject from NDC to camera coords
	std::tr1::array<GLfloat, 4> v = getProjectionMatrix().inverse() * n;
	
	//Perform the w divide
	for (size_t i(0); i < 4; ++i) v[i] /= v[3];
	
	//Unproject from camera to object space
	std::tr1::array<GLfloat, 4> w = getViewMatrix().inverse() * v;
	
	return magnet::math::Vector(w[0], w[1], w[2]);
      }

      /*! \brief Used to convert mouse positions (including depth
          information) into a 3D position.
       */
      math::Vector unprojectToDirection(int windowx, int windowy)
      {
	//We need to calculate the ray from the camera
	std::tr1::array<GLfloat, 4> n = {{(2.0 * windowx) / getWidth() - 1.0,
					  1.0 - (2.0 * windowy) / getHeight(), 
					  0.0, 1.0}};
	//Unproject from NDC to camera coords
	std::tr1::array<GLfloat, 4> v = getProjectionMatrix().inverse() * n;
	
	//Perform the w divide
	for (size_t i(0); i < 4; ++i) v[i] /= v[3];
	
	//Zero the w coordinate to stop the translations from the
	//viewmatrix affecting the vector
	v[3] = 0;

	//Unproject from camera to object space
	std::tr1::array<GLfloat, 4> w = getViewMatrix().inverse() * v;
	
	magnet::math::Vector vec(w[0], w[1], w[2]);
	vec /= vec.nrm();
	return vec;
      }

    protected:
      size_t _height, _width;
      float _panrotation;
      float _tiltrotation;
      math::Vector _nearPlanePosition;
      math::Vector _up;
      math::Vector _rotatePoint;
      
      GLfloat _zNearDist;
      GLfloat _zFarDist;
      math::Vector _eyeLocation;
      
      //! \brief One simulation length in cm (real units)
      double _simLength;

      //! \brief The diameter of a pixel, in cm
      double _pixelPitch;

      Camera_Mode _camMode;
    };
  }
}
