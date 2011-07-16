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
		       GLenum lightHandle = GL_LIGHT0,
		       GLfloat fovY = 45.0f,
		       GLfloat zNearDist = 0.05f, GLfloat zFarDist = 10.0f,
		       Vector up = Vector(0,1,0)):
	ViewPort(1,1,position, lookAtPoint, fovY, zNearDist, zFarDist, up),
	_lightHandle(lightHandle)
      {
	const GLfloat white[] = {1.0f, 1.0f, 1.0f, 1.0f};
	//Setup a bright light
	glLightfv(lightHandle, GL_DIFFUSE, white);
	glLightfv(lightHandle, GL_SPECULAR, white);
      }

      /*! \brief Updates the associated OpenGL light with the current light location.
       */
      inline void glUpdateLight()
      {
	Vector cameraLocation(getEyeLocation());
	
	GLfloat light_position[4] = {cameraLocation[0], cameraLocation[1], cameraLocation[2], 1.0f};
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
      }

      /*! \brief Renders the light source as a cone in the OpenGL scene
       */
      inline void drawLight()
      {
	Context::getContext().color(1.0f, 1.0f, 1.0f);
	
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
	
	buildMatrices();
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
      inline void loadShadowTextureMatrix(int textureUnit)
      {
	glActiveTextureARB(GL_TEXTURE0 + textureUnit);
	glMatrixMode(GL_TEXTURE);
	//Build the texture matrix
	glLoadIdentity();
	glTranslatef(0.5f, 0.5f, 0.5f);
	glScalef(0.5f, 0.5f, 0.5f);
	glMultMatrixf(_projectionMatrix);
	glMultMatrixf(_viewMatrix);

	GLfloat vp_viewMatrix[4*4];
	glGetFloatv(GL_MODELVIEW_MATRIX , vp_viewMatrix);

	invert(vp_viewMatrix);
	glMultMatrixf(vp_viewMatrix);
	glMatrixMode(GL_MODELVIEW);	  
      }

      /*! \brief Function to in-place invert an 4x4 OpenGL matrix.
       */
      static void invert(GLfloat* mat4x4)
      {
	GLfloat result[4*4];

	GLfloat tmp[12];											//temporary pair storage

	//calculate pairs for first 8 elements (cofactors)
	tmp[0] = mat4x4[10] * mat4x4[15];
	tmp[1] = mat4x4[11] * mat4x4[14];
	tmp[2] = mat4x4[9] * mat4x4[15];
	tmp[3] = mat4x4[11] * mat4x4[13];
	tmp[4] = mat4x4[9] * mat4x4[14];
	tmp[5] = mat4x4[10] * mat4x4[13];
	tmp[6] = mat4x4[8] * mat4x4[15];
	tmp[7] = mat4x4[11] * mat4x4[12];
	tmp[8] = mat4x4[8] * mat4x4[14];
	tmp[9] = mat4x4[10] * mat4x4[12];
	tmp[10] = mat4x4[8] * mat4x4[13];
	tmp[11] = mat4x4[9] * mat4x4[12];

	//calculate first 8 elements (cofactors)
	result[0] = tmp[0]*mat4x4[5] + tmp[3]*mat4x4[6] + tmp[4]*mat4x4[7]-tmp[1]*mat4x4[5] - tmp[2]*mat4x4[6] - tmp[5]*mat4x4[7];
	result[1] = tmp[1]*mat4x4[4] + tmp[6]*mat4x4[6] + tmp[9]*mat4x4[7]-tmp[0]*mat4x4[4] - tmp[7]*mat4x4[6] - tmp[8]*mat4x4[7];
	result[2] = tmp[2]*mat4x4[4] + tmp[7]*mat4x4[5] + tmp[10]*mat4x4[7]-tmp[3]*mat4x4[4] - tmp[6]*mat4x4[5] - tmp[11]*mat4x4[7];
	result[3] = tmp[5]*mat4x4[4] + tmp[8]*mat4x4[5] + tmp[11]*mat4x4[6]-tmp[4]*mat4x4[4] - tmp[9]*mat4x4[5] - tmp[10]*mat4x4[6];
	result[4] = tmp[1]*mat4x4[1] + tmp[2]*mat4x4[2] + tmp[5]*mat4x4[3]-tmp[0]*mat4x4[1] - tmp[3]*mat4x4[2] - tmp[4]*mat4x4[3];
	result[5] = tmp[0]*mat4x4[0] + tmp[7]*mat4x4[2] + tmp[8]*mat4x4[3]-tmp[1]*mat4x4[0] - tmp[6]*mat4x4[2] - tmp[9]*mat4x4[3];
	result[6] = tmp[3]*mat4x4[0] + tmp[6]*mat4x4[1] + tmp[11]*mat4x4[3]-tmp[2]*mat4x4[0] - tmp[7]*mat4x4[1] - tmp[10]*mat4x4[3];
	result[7] = tmp[4]*mat4x4[0] + tmp[9]*mat4x4[1] + tmp[10]*mat4x4[2]-tmp[5]*mat4x4[0] - tmp[8]*mat4x4[1] - tmp[11]*mat4x4[2];

	//calculate pairs for second 8 elements (cofactors)
	tmp[0] = mat4x4[2]*mat4x4[7];
	tmp[1] = mat4x4[3]*mat4x4[6];
	tmp[2] = mat4x4[1]*mat4x4[7];
	tmp[3] = mat4x4[3]*mat4x4[5];
	tmp[4] = mat4x4[1]*mat4x4[6];
	tmp[5] = mat4x4[2]*mat4x4[5];
	tmp[6] = mat4x4[0]*mat4x4[7];
	tmp[7] = mat4x4[3]*mat4x4[4];
	tmp[8] = mat4x4[0]*mat4x4[6];
	tmp[9] = mat4x4[2]*mat4x4[4];
	tmp[10] = mat4x4[0]*mat4x4[5];
	tmp[11] = mat4x4[1]*mat4x4[4];

	//calculate second 8 elements (cofactors)
	result[8 ] = tmp[0]*mat4x4[13] + tmp[3]*mat4x4[14] + tmp[4]*mat4x4[15]-tmp[1]*mat4x4[13] - tmp[2]*mat4x4[14] - tmp[5]*mat4x4[15];
	result[9 ] = tmp[1]*mat4x4[12] + tmp[6]*mat4x4[14] + tmp[9]*mat4x4[15]-tmp[0]*mat4x4[12] - tmp[7]*mat4x4[14] - tmp[8]*mat4x4[15];
	result[10] = tmp[2]*mat4x4[12] + tmp[7]*mat4x4[13] + tmp[10]*mat4x4[15]-tmp[3]*mat4x4[12] - tmp[6]*mat4x4[13] - tmp[11]*mat4x4[15];
	result[11] = tmp[5]*mat4x4[12] + tmp[8]*mat4x4[13] + tmp[11]*mat4x4[14]-tmp[4]*mat4x4[12] - tmp[9]*mat4x4[13] - tmp[10]*mat4x4[14];
	result[12] = tmp[2]*mat4x4[10] + tmp[5]*mat4x4[11] + tmp[1]*mat4x4[9]-tmp[4]*mat4x4[11] - tmp[0]*mat4x4[9] - tmp[3]*mat4x4[10];
	result[13] = tmp[8]*mat4x4[11] + tmp[0]*mat4x4[8] + tmp[7]*mat4x4[10]-tmp[6]*mat4x4[10] - tmp[9]*mat4x4[11] - tmp[1]*mat4x4[8];
	result[14] = tmp[6]*mat4x4[9] + tmp[11]*mat4x4[11] + tmp[3]*mat4x4[8]-tmp[10]*mat4x4[11] - tmp[2]*mat4x4[8] - tmp[7]*mat4x4[9];
	result[15] = tmp[10]*mat4x4[10] + tmp[4]*mat4x4[8] + tmp[9]*mat4x4[9]-tmp[8]*mat4x4[9] - tmp[11]*mat4x4[10] - tmp[5]*mat4x4[8];

	// calculate determinant
	GLfloat det = mat4x4[0]*result[0]+mat4x4[1]*result[1]+mat4x4[2]*result[2]
	  +mat4x4[3]*result[3];

	for (size_t i(0); i < 4*4; ++i)
	  result[i] /= det;

	if(det==0.0f)
	  for (size_t i(0); i < 4; ++i)
	    for (size_t j(0); j < 4; ++j)
	      result[4*i+j] = (i==j) ? 1 : 0;

	//Now we need to transpose the matrix as we copy out
	for (size_t i(0); i < 4; ++i)
	  for (size_t j(0); j < 4; ++j)
	    mat4x4[4*i+j] = result[4*j+i];
      }

      GLenum _lightHandle;
    };
  }
}
