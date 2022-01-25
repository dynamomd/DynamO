/*    dynamo:- Event driven molecular dynamics simulator 
 *    http://www.dynamomd.org
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
#include <magnet/GL/FBO.hpp>
#include <magnet/math/quaternion.hpp>
#include <magnet/clamp.hpp>
#include <magnet/exception.hpp>
#include <iostream>

namespace magnet
{
	namespace GL
	{
		class Camera
		{
			Context::ContextPtr _context;

		public:
			inline Camera(Context::ContextPtr context, GLfloat zNearDist = 0.3f,
						  GLfloat zFarDist = 300.0f) : _context(context),
													   _renderTarget(context),
													   _Gbuffer(context),
													   _hdrBuffer(context),
													   _luminanceBuffer1(context),
													   _luminanceBuffer2(context),
													   _blurTarget1(context),
													   _blurTarget2(context),
													   _filterTarget1(context),
													   _filterTarget2(context),
													   _height(1),
													   _width(1),
													   _zNearDist(zNearDist),
													   _zFarDist(zFarDist)
			{
				if (zNearDist > zFarDist)
					M_throw() << "zNearDist > _zFarDist!";
			}

			virtual GLMatrix getViewMatrix() const = 0;
			virtual GLMatrix getProjectionMatrix() const = 0;
			virtual void setUp(math::Vector newup, math::Vector axis = math::Vector{0, 0, 0}) = 0;

			FBO _renderTarget;
			FBO _Gbuffer;

			FBO _hdrBuffer;
			FBO _luminanceBuffer1;
			FBO _luminanceBuffer2;
			FBO _blurTarget1;
			FBO _blurTarget2;
			FBO _filterTarget1;
			FBO _filterTarget2;

			virtual FBO &getResolveBuffer()
			{
				return _renderTarget;
			}

			virtual void deinit()
			{
				_renderTarget.deinit();
				_Gbuffer.deinit();
				_hdrBuffer.deinit();
				_luminanceBuffer1.deinit();
				_luminanceBuffer2.deinit();
				_filterTarget1.deinit();
				_filterTarget2.deinit();
				_blurTarget1.deinit();
				_blurTarget2.deinit();
			}

			virtual void resize(size_t width, size_t height, size_t samples)
			{
				if ((_width == width) && (_height == height))
					return;

				deinit();
				_width = width;
				_height = height;

				{
					std::shared_ptr<magnet::GL::Texture2D> colorTexture(new magnet::GL::Texture2DMultisampled(_context, samples));
					colorTexture->init(width, height, GL_RGBA16F_ARB);

					std::shared_ptr<magnet::GL::Texture2D> normalTexture(new magnet::GL::Texture2DMultisampled(_context, samples));
					normalTexture->init(width, height, GL_RGBA16F_ARB);

					std::shared_ptr<magnet::GL::Texture2D> posTexture(new magnet::GL::Texture2DMultisampled(_context, samples));
					posTexture->init(width, height, GL_RGBA16F_ARB);

					std::shared_ptr<magnet::GL::Texture2D> depthTexture(new magnet::GL::Texture2DMultisampled(_context, samples));
					depthTexture->init(width, height, GL_DEPTH_COMPONENT);

					_Gbuffer.deinit();
					_Gbuffer.init();
					_Gbuffer.attachTexture(colorTexture, 0);
					_Gbuffer.attachTexture(normalTexture, 1);
					_Gbuffer.attachTexture(posTexture, 2);
					_Gbuffer.attachTexture(depthTexture);
				}

				{
					//Build the main/left-eye render buffer
					std::shared_ptr<magnet::GL::Texture2D> colorTexture(new magnet::GL::Texture2D(_context));
					colorTexture->init(width, height, GL_RGBA8);
					colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
					colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);

					std::shared_ptr<magnet::GL::Texture2D>
						depthTexture(new magnet::GL::Texture2D(_context));
					depthTexture->init(width, height, GL_DEPTH_COMPONENT);
					depthTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
					depthTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
					depthTexture->parameter(GL_TEXTURE_COMPARE_MODE, GL_NONE);

					_renderTarget.init();
					_renderTarget.attachTexture(colorTexture, 0);
					_renderTarget.attachTexture(depthTexture);
				}

				{
					std::shared_ptr<magnet::GL::Texture2D> colorTexture(new magnet::GL::Texture2D(_context));
					colorTexture->init(width, height, GL_RGBA);
					colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
					colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
					colorTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
					colorTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

					_filterTarget1.init();
					_filterTarget1.attachTexture(colorTexture, 0);
				}

				{
					std::shared_ptr<magnet::GL::Texture2D> colorTexture(new magnet::GL::Texture2D(_context));
					colorTexture->init(width, height, GL_RGBA);
					colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
					colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
					colorTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
					colorTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

					_filterTarget2.init();
					_filterTarget2.attachTexture(colorTexture, 0);
				}

				{
					std::shared_ptr<magnet::GL::Texture2D> colorTexture(new magnet::GL::Texture2D(_context));
					colorTexture->init(width / 4, height / 4, GL_RGB16F);
					colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
					colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
					colorTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
					colorTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

					_blurTarget1.init();
					_blurTarget1.attachTexture(colorTexture, 0);
				}

				{
					std::shared_ptr<magnet::GL::Texture2D> colorTexture(new magnet::GL::Texture2D(_context));
					colorTexture->init(width / 4, height / 4, GL_RGB16F);
					colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
					colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
					colorTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
					colorTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

					_blurTarget2.init();
					_blurTarget2.attachTexture(colorTexture, 0);
				}

				{
					std::shared_ptr<magnet::GL::Texture2D>
						colorTexture(new magnet::GL::Texture2D(_context));
					colorTexture->init(width, height, GL_RGBA16F);
					colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
					colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);

					std::shared_ptr<magnet::GL::Texture2D>
						depthTexture(new magnet::GL::Texture2D(_context));
					depthTexture->init(width, height,
									   GL_DEPTH_COMPONENT);
					depthTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
					depthTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
					depthTexture->parameter(GL_TEXTURE_COMPARE_MODE, GL_NONE);

					_hdrBuffer.init();
					_hdrBuffer.attachTexture(colorTexture, 0);
					_hdrBuffer.attachTexture(depthTexture);
				}

				{
					std::shared_ptr<magnet::GL::Texture2D>
						colorTexture(new magnet::GL::Texture2D(_context));

					colorTexture->init(width, height, GL_RGBA16F);
					colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
					colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);

					_luminanceBuffer1.init();
					_luminanceBuffer1.attachTexture(colorTexture, 0);
				}

				{
					std::shared_ptr<magnet::GL::Texture2D>
						colorTexture(new magnet::GL::Texture2D(_context));

					colorTexture->init(width / 2, height / 2, GL_RGBA16F);
					colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
					colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);

					_luminanceBuffer2.init();
					_luminanceBuffer2.attachTexture(colorTexture, 0);
				}
			}

			/*! \brief Get the normal matrix.
       
        \param offset This is an offset in camera coordinates to apply
        to the eye location. It's primary use is to calculate the
        perspective shift for the left and right eye in Analygraph
        rendering.
       */
			inline math::Matrix getNormalMatrix() const
			{
				return inverse(demoteToMatrix(getViewMatrix()));
			}

			/*! \brief Get the rotation part of the getViewMatrix(). */
			GLMatrix getViewRotationMatrix() const
			{
				//We can just cut it to a normal matrix, then re-promote to
				//get the rotation piece
				return magnet::GL::promoteToGLMatrix(magnet::GL::demoteToMatrix(getViewMatrix()));
			}

			/*! \brief Fetch the location of the users eyes, in object space
        coordinates.
        
        Useful for eye tracking applications. This returns the
        position of the eyes in object space by adding the eye
        location (relative to the viewing plane/screen) onto the
        current position.
       */
			math::Vector getPosition() const
			{
				auto inv = inverse(getViewMatrix());
				return math::Vector{inv(0, 3), inv(1, 3), inv(2, 3)};
			}

			//! \brief Get the aspect ratio of the screen
			inline GLfloat getAspectRatio() const
			{
				return ((GLfloat)_width) / _height;
			}

			//! \brief Get the height of the screen, in pixels.
			inline const size_t &getHeight() const { return _height; }

			//! \brief Get the width of the screen, in pixels.
			inline const size_t &getWidth() const { return _width; }

			/*! \brief Used to convert world positions to screen coordinates (pixels).
	
	This returns y coordinates in the format that cairo and other
	image programs expect (inverted compared to OpenGL).
	
	\return An array containing the x and y pixel locations,
	followed by the depth and w value.
       */
			magnet::math::NVector<GLfloat, 4> project(math::Vector invec) const
			{
				magnet::math::NVector<GLfloat, 4> vec = {GLfloat(invec[0]), GLfloat(invec[1]), GLfloat(invec[2]), 1.0f};
				vec = getProjectionMatrix() * (getViewMatrix() * vec);

				for (size_t i(0); i < 3; ++i)
					vec[i] /= std::abs(vec[3]);

				vec[0] = (0.5 + 0.5 * vec[0]) * getWidth();
				vec[1] = (0.5 - 0.5 * vec[1]) * getHeight();
				return vec;
			}

			/*! \brief Used to convert mouse positions (including depth
          information) into a 3D position.
       */
			math::Vector unprojectToPosition(int windowx, int windowy, GLfloat depth) const
			{
				//We need to calculate the ray from the camera
				magnet::math::NVector<GLfloat, 4> n = {(2.0f * windowx) / getWidth() - 1.0f,
													   1.0f - (2.0f * windowy) / getHeight(),
													   depth, 1.0f};
				//Unproject from NDC to camera coords
				magnet::math::NVector<GLfloat, 4> v = inverse(getProjectionMatrix()) * n;

				//Perform the w divide
				for (size_t i(0); i < 4; ++i)
					v[i] /= v[3];

				//Unproject from camera to object space
				magnet::math::NVector<GLfloat, 4> w = inverse(getViewMatrix()) * v;
				return magnet::math::Vector{w[0], w[1], w[2]};
			}

			/*! \brief Used to convert mouse positions (including depth
          information) into a 3D position.
       */
			math::Vector unprojectToDirection(int windowx, int windowy) const
			{
				//We need to calculate the ray from the camera
				magnet::math::NVector<GLfloat, 4> n = {(2.0f * windowx) / getWidth() - 1.0f,
													   1.0f - (2.0f * windowy) / getHeight(),
													   0.0f, 1.0f};
				//Unproject from NDC to camera coords
				magnet::math::NVector<GLfloat, 4> v = inverse(getProjectionMatrix()) * n;

				//Perform the w divide
				for (size_t i(0); i < 4; ++i)
					v[i] /= v[3];

				//Zero the w coordinate to stop the translations from the
				//viewmatrix affecting the vector
				v[3] = 0;

				//Unproject from camera to object space
				magnet::math::NVector<GLfloat, 4> w = inverse(getViewMatrix()) * v;

				math::Vector vec{w[0], w[1], w[2]};
				vec /= vec.nrm();
				return vec;
			}

		protected:
			size_t _height, _width;
			//! \brief Distance to the near clipping plane, in cm.
			GLfloat _zNearDist;
			//! \brief Distance to the far clipping plane, in cm.
			GLfloat _zFarDist;
		};

		/*! \brief An object to track the camera state.
     
      An OpenGL camera is a mapping between the object space (rendered
      object's natural coordinate system) and the screen space.

      We take this natural connection a little further and extend it
      from the screen space to the real space, as we would like to do
      interactive things like head tracking. 

      It actually turns out to be very convenient to define certain
      properties in terms of real space. For example, the near (\ref
      _zNearDist) and far (\ref _zFarDist) clipping planes can be
      defined once in real space and they don't have to be readjusted
      for different scenes. Its very natural to say that I don't want
      objects to appear closer to my eye than 8cm, and I would like to
      see all objects up to a distance of 10m.

      We need a length scale conversion (or zoom) factor for the
      conversion between the two spaces. This is provided by \ref
      _simLength. We also need to know the size of a single pixel on
      the screen to be able to accurately render objects, given by the
      pixel pitch (\ref _pixelPitch).

      This class can perform all the calculations required for setting
      up the projection and modelview matricies of the camera. There
      is also support for eye tracking calculations using the \ref
      _eyeLocation \ref math::Vector.
     */
		class CameraHeadTracking : public Camera
		{
		public:
			//! \brief The mode of the mouse movement
			enum Camera_Mode
			{
				ROTATE_CAMERA,
				ROTATE_POINT
			};

			/*! \brief The constructor.
	
        \param height The height of the viewport, in pixels.
        \param width The width of the viewport, in pixels.
        \param position The position of the screen (effectively the camera), in simulation coordinates.
        \param lookAtPoint The location the camera is initially focussed on.
        \param zNearDist The distance to the near clipping plane, in cm.
        \param zFarDist The distance to the far clipping plane, in cm.
        \param up A vector describing the up direction of the camera.
       */
			//We need a default constructor as viewPorts may be created without GL being initialized
			inline CameraHeadTracking(Context::ContextPtr context, math::Vector position = math::Vector{0, 0, 5},
									  math::Vector lookAtPoint = math::Vector{0, 0, 0},
									  GLfloat zNearDist = 8.0f,
									  GLfloat zFarDist = 10000.0f,
									  math::Vector up = math::Vector{0, 1, 0},
									  GLfloat simLength = 30.0f,
									  math::Vector eye_location = math::Vector{0, 0, 70}) : Camera(context, zNearDist, zFarDist),
																							_up(up.normal()),
																							_nearPlanePosition({0, 0, 0}),
																							_rotatePoint({0, 0, 0}),
																							_rotation(math::Quaternion::identity()),
																							_simLength(simLength),
																							_pixelPitch(0.04653), //cm per pixel (horizontal)
																							_camMode(ROTATE_POINT)
			{
				//We assume the user is around about 70cm from the screen
				setEyeLocation(eye_location);
				setPosition(position);
				lookAt(lookAtPoint);
			}

			inline void setRenderScale(double newscale) { _simLength = newscale; }

			inline GLfloat getRenderScale() const { return _simLength; }

			inline void lookAt(math::Vector lookAtPoint)
			{
				//Generate the direction from the near plane to the object
				const math::Vector oldEyePosition = getPosition();

				//Create the vectors at, right, and up, which are the target
				//rotated versions of (0,0,-1), (1,0,0), and (0,1,0)
				//respectively.
				math::Vector at = lookAtPoint - oldEyePosition;
				at.normalise();

				math::Vector right = at ^ _up;
				right.normalise();

				math::Vector up = right ^ at;
				up.normalise();

				//We now need to find the rotation into this target
				//set. Starting with rotating the at vector into position.
				_rotation = math::Quaternion::fromToVector(at, math::Vector{0, 0, -1});

				//Figure out where the right vector now is, then rotate that
				//into the correct position. This will not incorrectly "roll"
				//the view as it is guaranteed that right is perpendicular to
				//at, so up will be used as the rotation axis and will not
				//move.
				_rotation = math::Quaternion::fromToVector(right, _rotation * math::Vector{1, 0, 0}) * _rotation;
				_rotation = _rotation.inverse();

				//Finally, readjust the head position
				setPosition(oldEyePosition);
			}

			/*! \brief Get the modelview matrix. */
			virtual inline GLMatrix getViewMatrix() const
			{
				//Add in the movement of the eye and the movement of the
				//camera
				math::Vector cameraLocation = (_rotation.inverse().toMatrix() * _eyeLocation) / _simLength + _nearPlanePosition;

				//Setup the view matrix
				return promoteToGLMatrix(_rotation.toMatrix()) * translate(-cameraLocation);
			}

			/*! \brief Converts some inputted motion (e.g., by the mouse or keyboard) into a
        motion of the camera.

	All parameters may be negative or positive, as the sign
	defines the direction of the rotation/movement. Their name
	hints at what action they may do, depending on the camera mode
	(\ref _camMode).
       */
			inline void movement(float rotationX, float rotationY, float forwards, float sideways, float upwards)
			{
				forwards /= _simLength;
				sideways /= _simLength;
				upwards /= _simLength;

				math::Vector at = _rotation.inverse() * math::Vector{0, 0, -1};
				at.normalise();
				math::Vector up = _rotation.inverse() * math::Vector{0, 1, 0};
				up.normalise();
				math::Vector right = at ^ up;
				right.normalise();

				switch (_camMode)
				{
				case ROTATE_CAMERA:
				{
					//Move the camera
					math::Vector newpos = getPosition() + up * upwards + right * sideways + at * forwards;
					//Rotate the view
					math::Vector direction = math::Quaternion::fromAngleAxis(rotationY / 180.0, right) * math::Quaternion::fromAngleAxis(rotationX / 180.0, up) * at;

					setPosition(newpos);
					lookAt(newpos + direction);
					break;
				}
				case ROTATE_POINT:
				{
					lookAt(_rotatePoint);
					if (math::Vector(getPosition() - _rotatePoint).nrm() > forwards)
						_nearPlanePosition += forwards * at;
					rotationX -= 10 * sideways;
					rotationY += 10 * upwards;

					math::Vector offset = getPosition() - _rotatePoint;

					//We need to store the normal and restore it later.
					double offset_length = offset.nrm();

					offset = math::Quaternion::fromAngleAxis(-M_PI * rotationX / 180.0f, _up) * offset;

					math::Vector rotationAxis = up ^ offset;
					rotationAxis.normalise();
					offset = math::Quaternion::fromAngleAxis(-M_PI * rotationY / 180.0f, rotationAxis) * offset;
					offset.normalise();

					setPosition(offset_length * offset + _rotatePoint);
					lookAt(_rotatePoint);
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
					lookAt(getPosition() + axis);
					break;
				}
				case ROTATE_POINT:
				{
					double focus_distance = (getPosition() - _rotatePoint).nrm();
					setPosition(_rotatePoint - focus_distance * axis);
					lookAt(_rotatePoint);
					break;
				}
				default:
					M_throw() << "Bad camera mode";
				}
			}

			inline void setPosition(math::Vector newposition)
			{
				_nearPlanePosition = newposition - (getNormalMatrix() * _eyeLocation / _simLength);
			}

			inline void setRotatePoint(math::Vector vec)
			{
				if (_rotatePoint == vec)
					return;

				math::Vector shift = vec - _rotatePoint;
				_rotatePoint = vec;

				switch (_camMode)
				{
				case ROTATE_POINT:
					_nearPlanePosition += shift;
					lookAt(_rotatePoint);
					break;
				case ROTATE_CAMERA:
					break;
				default:
					M_throw() << "Unknown camera mode";
				}
			}

			/*! \brief Sets the eye location.
       
        \param eye The position of the viewers eye, relative to the
        center of the near viewing plane (in cm).
       */
			inline void setEyeLocation(math::Vector eye)
			{
				_eyeLocation = eye;
			}

			/*! \brief Gets the eye location (in cm).
       
        The position of the viewers eye is relative to the center of
        the near viewing plane (in cm).
       */
			virtual inline const math::Vector getEyeLocation() const
			{
				return _eyeLocation;
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
			virtual inline GLMatrix getProjectionMatrix() const
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
				return frustrum((-0.5f * getScreenPlaneWidth() - _eyeLocation[0]) * _zNearDist / _eyeLocation[2],  // left
								(+0.5f * getScreenPlaneWidth() - _eyeLocation[0]) * _zNearDist / _eyeLocation[2],  // right
								(-0.5f * getScreenPlaneHeight() - _eyeLocation[1]) * _zNearDist / _eyeLocation[2], // bottom
								(+0.5f * getScreenPlaneHeight() - _eyeLocation[1]) * _zNearDist / _eyeLocation[2], // top
								_zNearDist / _simLength,														   //Near distance
								_zFarDist / _simLength,															   //Far distance
								0.0f);
			}

			//! \brief Returns the screen's width (in simulation units).
			double getScreenPlaneWidth() const
			{
				return _pixelPitch * _width / _simLength;
			}

			//! \brief Returns the screen's height (in simulation units).
			double getScreenPlaneHeight() const
			{
				return _pixelPitch * _height / _simLength;
			}

			//! \brief Get the distance to the near clipping plane
			inline const GLfloat &getZNear() const { return _zNearDist; }

			//! \brief Get the distance to the far clipping plane
			inline const GLfloat &getZFar() const { return _zFarDist; }

			//! \brief Get the up direction of the camera.
			inline math::Vector getCameraUp() const
			{
				return getNormalMatrix() * math::Vector{0, 1, 0};
			}

			//! \brief Get the direction the camera is pointing in
			inline math::Vector getCameraDirection() const
			{
				return getNormalMatrix() * math::Vector{0, 0, -1};
			}

			/*! \brief Gets the pixel "diameter" in cm. */
			inline const double &getPixelPitch() const { return _pixelPitch; }

			/*! \brief Sets the pixel "diameter" in cm. */
			inline void setPixelPitch(double val) { _pixelPitch = val; }

			Camera_Mode getMode() const { return _camMode; }

			void setMode(Camera_Mode val) { _camMode = val; }

			/*! \brief set the orientation (roll) of the camera by setting
          its up direction.
	  
	  \param newup The new up direction.  
	  
	  \param axis If set, this will cause a rotation of the camera
	  position about the axis to compensates for the rotation of
	  the camera. It will look like the system is rotating but the
	  camera is remaining fixed. The axis must be 
       */
			virtual void setUp(math::Vector newup, math::Vector axis = math::Vector{0, 0, 0})
			{
				newup.normalise();
				if (axis.nrm2() != 0)
				{
					axis.normalise();
					math::Vector to = newup - axis * (axis | newup);
					math::Vector from = _up - axis * (axis | _up);
					setPosition(math::Quaternion::fromToVector(to.normal(), from.normal()) * getPosition());
				}
				_up = newup;
				movement(0, 0, 0, 0, 0);
			}

		protected:
			math::Vector _up;
			math::Vector _nearPlanePosition;
			math::Vector _rotatePoint;

			/*! \brief The location of the viewers eye, relative to the
	screen, in cm.
      */
			math::Vector _eyeLocation;
			math::Quaternion _rotation;

			//! \brief One simulation length in cm.
			double _simLength;

			//! \brief The diameter of a pixel, in cm
			double _pixelPitch;

			Camera_Mode _camMode;
		};
	}
}
