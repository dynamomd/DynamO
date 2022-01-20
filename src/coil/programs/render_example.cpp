/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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

#include <magnet/arg_share.hpp>
#include <magnet/GL/context.hpp>
#include <magnet/GL/buffer.hpp>
#include <magnet/GL/shader/render.hpp>
#include <magnet/GL/shader/resolver.hpp>
#include <magnet/GL/camera.hpp>
#include <magnet/image/PNG.hpp>
#include <iostream>

#include <SDL2/SDL.h>

/* A simple function that prints a message, the error code returned by SDL,
 * and quits the application */
void sdldie(const char *msg)
{
  printf("%s: %s\n", msg, SDL_GetError());
  SDL_Quit();
  exit(1);
}

void checkSDLError(int line = -1)
{
#ifndef NDEBUG
  const char *error = SDL_GetError();
  if (*error != '\0')
  {
    printf("SDL Error: %s\n", error);
    if (line != -1)
      printf(" + line: %i\n", line);
    SDL_ClearError();
  }
#endif
}

using namespace std;

int main(int argc, char *argv[])
{
  //Setup the GL context
  magnet::ArgShare::getInstance().setArgs(argc, argv);

  ////////////////////////////////////////////////////////
  /////     SDL/PLATFORM SPECIFIC CODE  //////////////////
  ////////////////////////////////////////////////////////


  SDL_Window *mainwindow;    /* Our window handle */
  SDL_GLContext maincontext; /* Our opengl context handle */

  if (SDL_Init(SDL_INIT_VIDEO) < 0)     /* Initialize SDL's Video subsystem */
    sdldie("Unable to initialize SDL"); /* Or die on error */

  /* Request opengl 3.2 context.
     * SDL doesn't have the ability to choose which profile at this time of writing,
     * but it should default to the core profile */
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 2);

  /* Create our window centered at 512x512 resolution */
  mainwindow = SDL_CreateWindow("COIL Render Example", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                512, 512, SDL_WINDOW_OPENGL /*| SDL_WINDOW_SHOWN */);
  if (!mainwindow) /* Die if creation failed */
    sdldie("Unable to create window");

  checkSDLError(__LINE__);

  /* Create our opengl context and attach it to our window */
  maincontext = SDL_GL_CreateContext(mainwindow);
  checkSDLError(__LINE__);

  ////////////////////////////////////////////////////////
  /////     COIL SPECIFIC CODE  //////////////////////////
  ////////////////////////////////////////////////////////

  auto _glContext = magnet::GL::Context::getContext();

  //Make a camera for rendering (contains everything for the render)
  //Set the near and far rendering distances (can't be zero and inf for numerical reasons!)
  magnet::GL::CameraHeadTracking _camera(0.3f, 300.0f);
  //Set the camera resolution and number of anti-aliasing samples
  _camera.resize(200, 200, 1);

  //Make a rendering shader and compile it
  magnet::GL::shader::RenderShader _renderShader;
  _renderShader.build();

  magnet::GL::shader::ResolverShader _resolverShader;
  _resolverShader.build();

  //Make the triangles for rendering, first we list the vertice indexes
  magnet::GL::Buffer<GLuint> _elementBuff;
  {
    std::vector<GLuint> Elements = {0, 1, 2};
    //Here we have to say how many elements per thing being drawn
    _elementBuff.init(Elements, 3);
  }

  magnet::GL::Buffer<GLfloat> _posBuff;
  {
    std::vector<float> VertexPos = {
        0.0f, 0.0f, 0.0f, //x1, y1, z1
        1.0f, 0.0f, 0.0f, //x2, ....
        0.0f, 1.0f, 0.0f};
    _posBuff.init(VertexPos, 3, magnet::GL::buffer_usage::STREAM_DRAW);
  }

  magnet::GL::Buffer<GLfloat> _normBuff;
  {
    std::vector<float> VertexNormals = {
        0.0f, 0.0f, 1.0f,
        0.0f, 0.0f, 1.0f,
        0.0f, 0.0f, 1.0f};
    _normBuff.init(VertexNormals, 3);
  }

  magnet::GL::Buffer<GLubyte> _colBuff;
  {
    std::vector<GLubyte> VertexColor = {
        255, 0, 0, 255, //Red, Green, Blue, Alpha (transparency)
        0, 255, 0, 255,
        0, 0, 255, 255};
    _colBuff.init(VertexColor, 4, magnet::GL::buffer_usage::STREAM_DRAW);
  }

  //Now attach the camera ready for rendering and clear it
  _camera.setPosition(magnet::math::Vector{0, 0, 50});
  _camera.lookAt(magnet::math::Vector{0, 0, 0});
  _camera._Gbuffer.attach();
  glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
  _glContext->setDepthTest(true);
  _glContext->setBlend(false);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  _renderShader.attach();
  _renderShader["ProjectionMatrix"] = _camera.getProjectionMatrix();
  _renderShader["ViewMatrix"] = _camera.getViewMatrix();
  _colBuff.attachToColor();
  _normBuff.attachToNormal();
  _posBuff.attachToVertex();
  _elementBuff.drawElements(magnet::GL::element_type::TRIANGLES);
  _renderShader.detach();
  _camera._Gbuffer.detach();
  _glContext->cleanupAttributeArrays();

  //Now we have buffers with all the information in them. They are unfortunately, anti-aliased buffers by default, so we need to "resolve" them into normal textures.
	  std::shared_ptr<magnet::GL::Texture2D> resolveTexture_color(new magnet::GL::Texture2D);
	  resolveTexture_color->init(_camera.getWidth(), _camera.getHeight(), GL_RGBA16F_ARB);
	  resolveTexture_color->parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	  resolveTexture_color->parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	  std::shared_ptr<magnet::GL::Texture2D> resolveTexture_pos(new magnet::GL::Texture2D);
	  resolveTexture_pos->init(_camera.getWidth(), _camera.getHeight(), GL_RGBA16F_ARB);
	  resolveTexture_pos->parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	  resolveTexture_pos->parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    magnet::GL::FBO _resolveTarget;
    _resolveTarget.init();

    //The color buffer
	  _resolveTarget.attachTexture(resolveTexture_color, 0);
    _resolveTarget.attach();
    _resolverShader.attach();
    _camera._Gbuffer.getColorTexture(0)->bind(0);
    _resolverShader["inTex"] = 0;
    _resolverShader["sample"] = GLint(0);
    _resolverShader.invoke();
    _resolverShader.detach();
    _resolveTarget.detach();

    //The position buffer
    _resolveTarget.attachTexture(resolveTexture_pos, 0);
    _resolveTarget.attach();
    _resolverShader.attach();
    _camera._Gbuffer.getColorTexture(2)->bind(0);
    _resolverShader["inTex"] = 0;
    _resolverShader["sample"] = GLint(0);
    _resolverShader.invoke();
    _resolverShader.detach();
    _resolveTarget.detach();

  {
    std::vector<GLfloat> pixels;
	  pixels.resize(_camera.getWidth() * _camera.getHeight() * 4);
	  resolveTexture_color->writeto(pixels);

    std::vector<uint8_t > png_pixels;
	  png_pixels.resize(_camera.getWidth() * _camera.getHeight() * 4);
    
    for (size_t x(0); x < _camera.getWidth(); ++x)
     for (size_t y(0); y < _camera.getHeight(); ++y)
      for (size_t comp(0); comp < 4; ++comp) {
        size_t idx =  comp + 4 * ( x + _camera.getWidth() * y);
        png_pixels[idx] = uint8_t(pixels[idx] * 255);
      }

    for (size_t x(0); x < _camera.getWidth(); ++x)
     for (size_t y(0); y < _camera.getHeight(); ++y)
      if (pixels[4* _camera.getWidth() * y + 4 * x + 3] != 0)
        std::cout << x << "," << y << " <" << pixels[4* _camera.getWidth() * y + 4 * x + 0] <<  "," << pixels[4* _camera.getWidth() * y + 4 * x + 1] <<  "," << pixels[4* _camera.getWidth() * y + 4 * x + 2] << ">\n";

    magnet::image::writePNGFile("color.png", png_pixels, _camera.getWidth(), _camera.getHeight(), 4, 1, true, true);
  }

  {
    std::vector<GLfloat> pixels;
	  pixels.resize(_camera.getWidth() * _camera.getHeight() * 4);
	  resolveTexture_pos->writeto(pixels);

    for (size_t x(0); x < _camera.getWidth(); ++x)
     for (size_t y(0); y < _camera.getHeight(); ++y)
      if (pixels[4* _camera.getWidth() * y + 4 * x + 3] != 0)
        std::cout << x << "," << y << " <" << pixels[4* _camera.getWidth() * y + 4 * x + 0] <<  "," << pixels[4* _camera.getWidth() * y + 4 * x + 1] <<  "," << pixels[4* _camera.getWidth() * y + 4 * x + 2] << ">\n";
  }

  /* Delete our opengl context, destroy our window, and shutdown SDL */
  SDL_GL_DeleteContext(maincontext);
  SDL_DestroyWindow(mainwindow);
  SDL_Quit();
  return 0;
}
