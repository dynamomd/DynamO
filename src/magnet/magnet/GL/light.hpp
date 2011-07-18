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

#include <limits>
#include <magnet/GL/viewPort.hpp>
#include <magnet/GL/context.hpp>

namespace magnet {
  namespace GL {
    /*! \brief A specialization of the \ref ViewPort, to ease creating
     * shadow mapping light sources
     *
     * A shadow mapping light source is an OpenGL light source which
     * also requires the depth map of the scene, rendered from its
     * location, to perform shadow casting.
     *
     * This type of light is directional as it has a direction it is
     * looking at, just like the camera.
     */
    class LightInfo: public ViewPort
    {
    public:
      /*! \brief Default constructor
       *
       * We need a default constructor as ViewPort classes may be
       * created without GL being initialized.
       */
      inline LightInfo():
	ViewPort(1,1) {}

      /*! \brief Constructor
       *
       * This constructor also associates an OpenGL lightsource with this class.
       * \param lightHandle The enum of the OpenGL light source associated with this class.
       * \param position The position of the screen (effectively the camera), in simulation coordinates.
       * \param lookAtPoint The location the camera is initially focussed on.
       * \param fovY The field of vision of the camera.
       * \param zNearDist The distance to the near clipping plane.
       * \param zFarDist The distance to the far clipping plane.
       * \param up A vector describing the up direction of the camera.
       */
      inline LightInfo(Vector position, 
		       Vector lookAtPoint,
		       GLfloat fovY = 45.0f,
		       GLfloat zNearDist = 0.05f, GLfloat zFarDist = 10.0f,
		       Vector up = Vector(0,1,0)):
	ViewPort(1,1,position, lookAtPoint, fovY, zNearDist, zFarDist, up)
      {}

      /*! \brief Renders the light source as a cone in the OpenGL scene
       */
      inline void drawLight()
      {
	Context& context = Context::getContext();
	context.color(1, 1, 1);
	
	GLfloat rotationAngle 
	  = (180.0 / M_PI) * std::acos(Vector(0,0,-1) | getCameraDirection());
	
	Vector RotationAxis = Vector(0,0,-1) ^ getCameraDirection();
	float norm = RotationAxis.nrm();
	RotationAxis /= norm;
	if (norm < std::numeric_limits<double>::epsilon())
	  RotationAxis = Vector(1,0,0);

	Vector cameraLocation(getEyeLocation());
	
	GLMatrix oldviewMatrix = context.getViewMatrix();

	context.setViewMatrix(oldviewMatrix 
			      * GLMatrix::translate(cameraLocation[0], cameraLocation[1], cameraLocation[2])
			      * GLMatrix::rotate(rotationAngle, RotationAxis)
			      * GLMatrix::translate(0,0,0.0025));

	GLfloat r = 0.05f;
	double fovY = getFOVY();
	glutSolidCone(r * std::sin(fovY * M_PI / 360.0f), 
		      r * std::cos(fovY * M_PI / 360.0f), 
		      15, 15);

	context.setViewMatrix(oldviewMatrix);
      }

      /*! \brief Allow copying the lights location from a \ref ViewPort.
       * 
       * This operation does not copy all details of the \ref ViewPort. For
       * example, the field of view and aspect ratio of the light is
       * maintained.
       */
      LightInfo& operator=(const ViewPort& vp)
      {
	//The FOV and the aspect ratio of the light must be maintained
	size_t width = _width;
	size_t height = _height;
	double fovY = getFOVY();
	ViewPort::operator=(vp);
	setFOVY(fovY);
	_width = width;
	_height = height;
	
	return *this; 
      }

      /*! \brief Load the specified OpenGL texture matrix with the
       * projection required for shadow mapping.
       * 
       * \note The current OpenGL model view matrix must be the matrix
       * used for rendering.
       *
       * \param textureUnit The texture unit whose matrix is to be
       * setup for shadowmapping.
       */
      inline GLMatrix getShadowTextureMatrix()
      {
	return GLMatrix::translate(Vector(0.5, 0.5, 0.5))
	  * GLMatrix::scale(Vector(0.5, 0.5, 0.5))
	  * getProjectionMatrix()
	  * getViewMatrix()
	  * Context::getContext().getViewMatrix().inverse();
      }
    };
  }
}
