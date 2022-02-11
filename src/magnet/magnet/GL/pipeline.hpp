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

//Here we have the correct order of GL includes
#include <magnet/GL/context.hpp>
#include <magnet/GL/camera.hpp>
#include <magnet/GL/shader/lightShader.hpp>
#include <magnet/GL/shader/ambientLight.hpp>
#include <magnet/GL/shader/luminance.hpp>
#include <magnet/GL/shader/blur.hpp>
#include <magnet/GL/shader/toneMap.hpp>
#include <magnet/GL/shader/depthResolver.hpp>
#include <magnet/GL/shader/copy.hpp>
#include <magnet/GL/shader/downsampler.hpp>
#include <magnet/GL/actor.hpp>
#include <magnet/GL/light.hpp>

namespace magnet
{
    namespace GL
    {
        class Pipeline
        {
        protected:
            Context::ContextPtr _context;

            shader::PointLightShader _pointLightShader;
            shader::ShadowLightShader _shadowLightShader;
            shader::AmbientLightShader _ambientLightShader;
            shader::LuminanceShader _luminanceShader;
            shader::LuminanceMipMapShader _luminanceMipMapShader;
            shader::ToneMapShader _toneMapShader;
            shader::DownsamplerShader _downsampleShader;
            shader::SeperableGaussian _blurShader;
            shader::DepthResolverShader _depthResolverShader;
            shader::CopyShader _copyShader;

            magnet::GL::FBO _shadowbuffer;

        public:
            Pipeline(Context::ContextPtr context) : _context(context),
                                                 _shadowbuffer(context)
            {
            }

            void init(int shadowbuf_size = 1024)
            {
                _copyShader.build();
                _downsampleShader.build();
                _blurShader.build();
                _pointLightShader.build();
                _shadowLightShader.build();
                _ambientLightShader.build();
                _luminanceShader.build();
                _luminanceMipMapShader.build();
                _toneMapShader.build();
                _depthResolverShader.build();

                {
                    //Build depth buffer
                    std::shared_ptr<magnet::GL::Texture2D> depthTexture(new magnet::GL::Texture2D(_context));
                    //We don't force GL_DEPTH_COMPONENT24 as it is likely you get
                    //the best precision anyway
                    depthTexture->init(shadowbuf_size, shadowbuf_size, GL_DEPTH_COMPONENT); //SIZE MUST BE THE SAME FOR THE LIGHTS
                    //You must select GL_NEAREST for depth data, as GL_LINEAR
                    //converts the value to 8bit for interpolation (on NVidia).
                    depthTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
                    depthTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
                    depthTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
                    depthTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
                    depthTexture->parameter(GL_TEXTURE_COMPARE_MODE, GL_NONE);

                    //Build color texture
                    std::shared_ptr<magnet::GL::Texture2D> colorTexture(new magnet::GL::Texture2D(_context));
                    colorTexture->init(shadowbuf_size, shadowbuf_size, GL_RG32F); //SIZE MUST BE THE SAME FOR THE LIGHTS
                    colorTexture->parameter(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
                    colorTexture->parameter(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
                    colorTexture->parameter(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
                    colorTexture->parameter(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);

                    _shadowbuffer.init();
                    _shadowbuffer.attachTexture(colorTexture, 0);
                    _shadowbuffer.attachTexture(depthTexture);
                }
            }

            void deinit()
            {
                _shadowbuffer.deinit();
                _toneMapShader.deinit();
                _depthResolverShader.deinit();
                _pointLightShader.deinit();
                _shadowLightShader.deinit();
                _ambientLightShader.deinit();
                _downsampleShader.deinit();
                _blurShader.deinit();
                _copyShader.deinit();
                _luminanceShader.deinit();
                _luminanceMipMapShader.deinit();
            }

            void render(Camera &camera, GLfloat ambient_light_level, std::vector<std::shared_ptr<Actor>> actors)
            {
                FBO &renderTarget = camera.getResolveBuffer();

                ///////////////////////Deferred Shading G-Buffer Creation /////////////////
                //We perform a deffered shading pass followed by a forward shading
                //pass for objects which cannot be deferred, like volumes etc.
                //We use the stencil buffer to track which pixels should be shaded
                //in the deferred pass.

                camera._Gbuffer.attach();
                glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
                _context->setDepthTest(true);
                _context->setBlend(false);
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

                //Enter the render ticks for all objects
                for (auto &obj : actors)
                    if (obj->visible())
                        obj->deferred_pass(camera, detail::RenderMode::DEFAULT);

                camera._Gbuffer.detach();

                ///////////////////////Lighting pass////////////////////////
                //Here we calculate the lighting of every pixel in the scene
                camera._Gbuffer.getColorTexture(0)->bind(0);
                camera._Gbuffer.getColorTexture(1)->bind(1);
                camera._Gbuffer.getColorTexture(2)->bind(2);

                //First, set up the buffers for rendering
                camera._hdrBuffer.attach();
                glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
                //We need the depth test on, to enable writes to the depth buffer
                _context->setDepthTest(true);

                //Now we need to populate the depth buffer, but nothing needs to
                //go in the color buffer
                glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
                _depthResolverShader.attach();
                _depthResolverShader["posTex"] = 2;
                _depthResolverShader["samples"] = GLint(camera.getSamples());
                _depthResolverShader["ProjectionMatrix"] = camera.getProjectionMatrix();
                _depthResolverShader.invoke();
                _depthResolverShader.detach();
                glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

                //Additive blending of all of the lights contributions, except for
                //the alpha values
                _context->setBlend(true);
                glBlendFuncSeparate(GL_ONE, GL_ONE, GL_ONE, GL_ZERO);

                //Now disable writing or testing of the depth buffer
                _context->setDepthTest(false);
                glDepthMask(GL_FALSE);

                _ambientLightShader.attach();
                _ambientLightShader["colorTex"] = 0;
                _ambientLightShader["samples"] = GLint(camera.getSamples());
                _ambientLightShader["ambientLight"] = ambient_light_level;
                _ambientLightShader.invoke();
                _ambientLightShader.detach();

                //We need to collect the lights for the forward shading pass
                std::vector<std::shared_ptr<Light>> lights;

                //Perform the shadow casting lights
                for (auto &light_obj : actors)
                {
                    std::shared_ptr<Light> light = std::dynamic_pointer_cast<Light>(light_obj);
                    if (!light)
                    continue;
                    lights.push_back(light);

                    if (!(light->shadowCasting()))
                        continue;

                    //Change from the hdr FBO
                    camera._hdrBuffer.detach();
                    //Render each light's shadow map
                    _shadowbuffer.attach();
                    _context->setDepthTest(true);
                    glDepthMask(GL_TRUE);
                    _context->setBlend(false);

                    glClearColor(light->getZFar(), light->getZFar() * light->getZFar(), 0, 0);
                    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

                    for (auto &obj : actors)
                        if (obj->visible() && obj->shadowCasting())
                            obj->deferred_pass(*light, detail::RenderMode::SHADOW);

                    _shadowbuffer.detach();

                    //_shadowbuffer.getColorTexture(0)->genMipmaps();
                    _shadowbuffer.getColorTexture(0)->bind(7);

                    //Change back to the hdr FBO
                    camera._hdrBuffer.attach();
                    _context->setDepthTest(false);
                    glDepthMask(GL_FALSE);
                    _context->setBlend(true);
                    _shadowLightShader.attach();
                    _shadowLightShader["colorTex"] = 0;
                    _shadowLightShader["normalTex"] = 1;
                    _shadowLightShader["positionTex"] = 2;
                    _shadowLightShader["shadowTex"] = 7;
                    _shadowLightShader["shadowMatrix"] = light->getShadowTextureMatrix() * inverse(camera.getViewMatrix());
                    _shadowLightShader["samples"] = GLint(camera.getSamples());
                    _shadowLightShader["lightColor"] = light->getLightColor();
                    _shadowLightShader["lightSpecularExponent"] = light->getSpecularExponent();
                    _shadowLightShader["lightSpecularFactor"] = light->getSpecularFactor();
                    _shadowLightShader["lightPosition"] = light->getEyespacePosition(camera);
                    _shadowLightShader["maxVariance"] = light->getMaxVariance();
                    _shadowLightShader["bleedReduction"] = light->getBleedReduction();
                    _shadowLightShader.invoke();
                    _shadowLightShader.detach();
                }

                //Perform the point/non-shadowing lights
                _pointLightShader.attach();
                _pointLightShader["colorTex"] = 0;
                _pointLightShader["normalTex"] = 1;
                _pointLightShader["positionTex"] = 2;
                _pointLightShader["samples"] = GLint(camera.getSamples());
                for (auto &light : lights)
                {
                    if (light->shadowCasting())
                        continue;
                    
                    _pointLightShader["lightColor"] = light->getLightColor();
                    _pointLightShader["lightSpecularExponent"] = light->getSpecularExponent();
                    _pointLightShader["lightSpecularFactor"] = light->getSpecularFactor();
                    _pointLightShader["lightPosition"] = light->getEyespacePosition(camera);
                    _pointLightShader.invoke();
                }

                _pointLightShader.detach();
            

                _context->setBlend(true);
                _context->setDepthTest(true);
                glDepthMask(GL_TRUE);
                glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

                //Enter the forward render ticks for all objects
                for (auto &obj : actors)
                    if (obj->visible())
                        obj->forward_pass(camera, lights, ambient_light_level, detail::RenderMode::DEFAULT);

                camera._hdrBuffer.detach();
                ///////////////////////Luminance Sampling//////////////////////
                //The light buffer is bound to texture unit 0 for the tone mapping too
                _context->setDepthTest(false);
                _context->setBlend(false);

                camera._hdrBuffer.getColorTexture()->bind(0);

                camera._luminanceBuffer1.attach();
                _luminanceShader.attach();
                _luminanceShader["colorTex"] = 0;
                _luminanceShader.invoke();
                _luminanceShader.detach();
                camera._luminanceBuffer1.detach();

                //    std::vector<GLfloat> data;
                //    _luminanceBuffer1.getColorTexture()->writeto(data);
                //
                //    double avg = 0, max = data[1], min = data[2], count = 0;
                //    for (size_t i(0); i < data.size() / 4; ++i)
                //      {
                //	if (data[4*i+3])
                //	  {
                //	    if (!count)
                //	      { max = data[4*i+1]; min = data[4*i+2]; }
                //
                //	    avg += data[4*i+0] * data[4*i+3];
                //	    max = std::max(double(data[4*i+1]), max);
                //	    min = std::min(double(data[4*i+2]), min);
                //	    count += data[4*i+3];
                //	  }
                //      }
                //
                //    avg /= count;
                //    std::cout << "Luminance Avg="<< avg << ", max=" << max << ", min=" << min << ", count=" << count * 10000.0 << "\n";

                FBO *luminanceSource = &camera._luminanceBuffer1;
                FBO *luminanceDestination = &camera._luminanceBuffer2;

                //Now we need to generate the mipmaps containing the scene
                //average, minimum and maximum luminances.
                {
                    GLsizei currentWidth = camera._luminanceBuffer1.getColorTexture()->getWidth();
                    GLsizei currentHeight = camera._luminanceBuffer1.getColorTexture()->getHeight();
                    GLint numLevels = camera._luminanceBuffer1.getColorTexture()->calcMipmapLevels();

                    //Attach the mipmapping shader
                    _luminanceMipMapShader.attach();
                    for (int i = 1; i < numLevels; ++i)
                    {
                        luminanceDestination->attach();
                        luminanceSource->getColorTexture()->bind(0);
                        _luminanceMipMapShader["inputTex"] = 0;

                        std::array<GLint, 2> oldSize = {{currentWidth, currentHeight}};
                        _luminanceMipMapShader["oldSize"] = oldSize;

                        //Halve the size of the textures, ensuring they never drop below 1
                        currentWidth /= 2;
                        currentWidth += !currentWidth;
                        currentHeight /= 2;
                        currentHeight += !currentHeight;
                        _context->setViewport(0, 0, currentWidth, currentHeight);

                        //Now generate the mipmap level using a shader
                        _luminanceMipMapShader.invoke();

                        luminanceDestination->detach();
                        std::swap(luminanceSource, luminanceDestination);
                    }
                    _luminanceMipMapShader.detach();
                }

                ///////////////////////Blurred Pipeline///////////////////////////////
                if (camera._bloomLighting)
                {
                    Texture2D &tex = *camera._hdrBuffer.getColorTexture();
                    tex.bind(0);

                    camera._blurTarget1.attach();
                    _downsampleShader.attach();
                    _downsampleShader["inputTex"] = 0;
                    _downsampleShader["downscale"] = GLint(4);
                    std::array<GLint, 2> oldSize = {{tex.getWidth(), tex.getHeight()}};
                    _downsampleShader["oldSize"] = oldSize;
                    _downsampleShader.invoke();
                    _downsampleShader.detach();
                    camera._blurTarget1.detach();

                    _blurShader.attach();
                    _blurShader["colorTex"] = 0;
                    std::array<GLfloat, 2> invDim = {{1.0f / (tex.getWidth() / 4),
                                                      1.0f / (tex.getHeight() / 4)}};
                    _blurShader["invDim"] = invDim;

                    for (size_t passes(0); passes < 1; ++passes)
                    {
                        camera._blurTarget1.getColorTexture()->bind(0);
                        camera._blurTarget2.attach();
                        _blurShader["direction"] = 0;
                        _blurShader.invoke();
                        _blurShader.detach();
                        camera._blurTarget2.detach();

                        camera._blurTarget2.getColorTexture()->bind(0);
                        camera._blurTarget1.attach();
                        _blurShader.attach();
                        _blurShader["direction"] = 1;
                        _blurShader.invoke();
                        camera._blurTarget1.detach();
                    }
                    _blurShader.detach();
                }

                ///////////////////////Tone Mapping///////////////////////////
                renderTarget.attach();
                camera._hdrBuffer.getColorTexture()->bind(0);
                luminanceSource->getColorTexture()->bind(1);
                if (camera._bloomLighting)
                    camera._blurTarget1.getColorTexture()->bind(2);
                _toneMapShader.attach();
                _toneMapShader["color_tex"] = 0;
                _toneMapShader["logLuma"] = 1;
                _toneMapShader["bloom_tex"] = 2;
                _toneMapShader["bloom_enable"] = camera._bloomLighting;
                _toneMapShader["bloomCompression"] = GLfloat(camera._bloomCompression);
                _toneMapShader["bloomCutoff"] = GLfloat(camera._bloomCutoff);
                _toneMapShader["Lpwhite"] = GLfloat(camera._bloomSaturation);
                _toneMapShader["scene_key"] = GLfloat(camera._sceneKey);
                _toneMapShader["background_color"] = camera._backColor;
                _toneMapShader.invoke();
                _toneMapShader.detach();
                renderTarget.detach();

                _context->setDepthTest(true);
            }
        };
    }
}