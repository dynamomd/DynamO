/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2011  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

    This program is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    version 3 as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Volume.hpp"
#include <iostream>
#include <coil/glprimatives/arrow.hpp>
#include <coil/RenderObj/console.hpp>

namespace coil {
  RVolume::RVolume(std::string name):
    RQuads(name)
  {}

  RVolume::~RVolume()
  {}

  void 
  RVolume::releaseCLGLResources()
  {
    if (_fbo.get() != NULL) _fbo->deinit();
    _fbo.release();
  }

  void 
  RVolume::initOpenGL() 
  {
    _shader.build();
    _fbo.reset(new magnet::GL::FBO);
    _fbo->init(_viewPort->getWidth(), _viewPort->getHeight());
  }

  void 
  RVolume::resize(size_t width, size_t height)
  {
    if (_fbo.get() != NULL) 
      _fbo->resize(width, height);
  }

  void 
  RVolume::initOpenCL() 
  {
    {
      float vertices[] = {-1,-1,-1,  1,-1,-1,  1, 1,-1, -1, 1,-1,
			  -1,-1, 1, -1, 1, 1,  1, 1, 1,  1,-1, 1};
      
      std::vector<float> VertexPos(vertices, vertices 
				   + sizeof(vertices) / sizeof(float));
      setGLPositions(VertexPos);
    }
    
    {
      int elements[] = {3,2,1,0, 6,7,1,2, 5,4,7,6, 3,0,4,5, 6,2,3,5, 7,4,0,1 };
      std::vector<int> ElementData(elements, elements + sizeof(elements) / sizeof(int));
      setGLElements(ElementData);
    }
  }  

  void 
  RVolume::glRender(magnet::GL::FBO& fbo)
  {
    //Before we render, we need the current depth buffer for depth testing.
    fbo.copyto(*_fbo, GL_DEPTH_BUFFER_BIT);

    //Now bind this copied depth texture to texture unit 0
    glActiveTextureARB(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, _fbo.getDepthTexture());

    //Now we can render
    GLhandleARB oldshader = glGetHandleARB(GL_PROGRAM_OBJECT_ARB);

    GLfloat FocalLength = 1.0f / std::tan(_viewPort->getFOVY() * (M_PI / 360.0f));

    _shader.attach(FocalLength, _viewPort->getWidth(), 
		   _viewPort->getHeight(), _viewPort->getEyeLocation());

    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    RQuads::glRender();
    glDisable(GL_CULL_FACE);

    glUseProgramObjectARB(oldshader);
  }
}
