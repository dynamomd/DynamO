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

#include <magnet/GL/camera.hpp>
#include <magnet/GL/multisampledFBO.hpp>
#include <limits>

namespace magnet {
  namespace GL {
    /*! \brief A specialization of the \ref Camera class, to ease
     * creating shadow mapping light sources
     *
     * A shadow mapping light source is an OpenGL light source which
     * also requires the depth map of the scene, rendered from its
     * location, to perform shadow casting.
     *
     * This type of light is directional as it has a direction it is
     * looking at, just like the camera.
     */
    class Light: public Camera
    {
    public:
      /*! \brief Default constructor
       *
       * We need a default constructor as Camera classes may be
       * created without GL being initialized.
       */
      inline Light():
	Camera(1,1) {}

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
      inline Light(math::Vector position, 
		   math::Vector lookAtPoint,
		   GLfloat fovY = 45.0f,
		   GLfloat zNearDist = 0.05f, GLfloat zFarDist = 10.0f,
		   math::Vector up = math::Vector(0,1,0)):
	Camera(1,1,position, lookAtPoint, fovY, zNearDist, zFarDist, up)
      {}

      /*! \brief Initialise the OpenGL resources of this light source. 
       */
      void init()
      {
	//Build depth buffer
	std::tr1::shared_ptr<Texture2D> depthTexture(new Texture2D);
	//We don't force GL_DEPTH_COMPONENT24 as it is likely you get
	//the best precision anyway
	depthTexture->init(1024, 1024, GL_DEPTH_COMPONENT);
	//You must select GL_NEAREST for depth data, as GL_LINEAR
	//converts the value to 8bit for interpolation (on NVidia).
	depthTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	depthTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	depthTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	depthTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	depthTexture->parameter(GL_TEXTURE_COMPARE_MODE, GL_NONE);

	//Build color texture
	std::tr1::shared_ptr<Texture2D> colorTexture(new Texture2D);
	colorTexture->init(1024, 1024, GL_RG32F);
	colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	colorTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	colorTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	_shadowFBO.init();
	_shadowFBO.attachTexture(colorTexture, 0);
	_shadowFBO.attachTexture(depthTexture);
      }

      /*! \brief Release the OpenGL resources of this light source. 
       */
      void deinit() { _shadowFBO.deinit(); }

      /*! \brief Allow copying the lights location from a \ref Camera.
       * 
       * This operation does not copy all details of the \ref Camera. For
       * example, the field of view and aspect ratio of the light is
       * maintained.
       */
      Light& operator=(const Camera& vp)
      {
	_position = vp.getViewPlanePosition();
	_panrotation = vp.getPan();
	_tiltrotation = vp.getTilt();
	return *this; 
      }

      /*! \brief Load the specified OpenGL texture matrix with the
        projection required for shadow mapping.
        
        \note The current OpenGL model view matrix must be the matrix
        used for rendering.
       
        \param textureUnit The texture unit whose matrix is to be
        setup for shadowmapping.
       */
      inline GLMatrix getShadowTextureMatrix()
      {
	return GLMatrix::translate(math::Vector(0.5, 0.5, 0.5))
	  * GLMatrix::scale(math::Vector(0.5, 0.5, 0.5))
	  * getProjectionMatrix()
	  * getViewMatrix();
      }

      /*! \brief Returns the frame buffer containing the shadow map.
       */
      FBO& shadowFBO() { return _shadowFBO; }

      /*! \brief Returns the texture containing the shadow map.
       */
      std::tr1::shared_ptr<Texture2D>& shadowTex() 
      { return _shadowFBO.getColorTexture(); }

    protected:
      FBO _shadowFBO;
    };
  }
}
