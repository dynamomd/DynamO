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
#pragma once

#include <gtkmm.h>

#include <coil/RenderObj/RenderObj.hpp>
#include <magnet/gtk/transferFunction.hpp>
#include <magnet/GL/texture.hpp>
#include <magnet/GL/buffer.hpp>
#include <magnet/GL/shader/volume.hpp>
#include <magnet/GL/shader/detail/ssshader.hpp>
#include <memory>
#include <array>

namespace coil {
  class RVolume : public RenderObj
  {
    class DepthCopyShader: public magnet::GL::shader::detail::SSShader
    {
    public:
#define STRINGIFY(A) #A
      virtual std::string initFragmentShaderSource()
      {
	return STRINGIFY(
uniform sampler2D depthTex;
void main()
{
  gl_FragDepth=texelFetch(depthTex, ivec2(gl_FragCoord.xy), 0).r;
});
      }
#undef STRINGIFY
  };

  public:
    RVolume(std::string name): RenderObj(name), _stepSizeVal(0.01), _dimensions({1,1,1}), _frameCounter(0) {}
    RVolume(stator::xml::Node);
  
    virtual bool deletable() { return true; }
    
    virtual void init(const std::shared_ptr<magnet::thread::TaskQueue>& systemQueue);
    virtual void deinit();
    virtual void forwardRender(magnet::GL::FBO& fbo,
			       const magnet::GL::Camera& cam,
			       std::vector<std::shared_ptr<RLight> >& lights,
			       GLfloat ambient,
			       RenderMode mode);

    virtual void showControls(Gtk::ScrolledWindow* win);

    void loadRawFile(std::string filename, std::array<size_t, 3> dimensions, size_t bytes, Vector dims = Vector{0,0,0});

#ifdef COIL_TIFFSUPPORT
    void loadTIFFFiles(std::vector<std::string> files);
#endif

    /*! \brief Return the icon used for the object in the render view.
     */
    virtual Glib::RefPtr<Gdk::Pixbuf> getIcon();

    virtual magnet::math::Vector getMaxCoord() const { return +0.5 * _dimensions; }
    virtual magnet::math::Vector getMinCoord() const { return -0.5 * _dimensions; }

    void loadData(const std::vector<GLubyte>& inbuffer, std::array<size_t, 3> dim, Vector scale);

    virtual void xml(stator::xml::Node);
  protected:
    void loadSphereTestPattern();

    void initGTK();
    void guiUpdate();

    void transferFunctionUpdated();

    magnet::GL::shader::VolumeShader _shader;
    magnet::GL::Buffer<GLfloat> _cubeVertices;
    magnet::GL::FBO _currentDepthFBO;
    DepthCopyShader _depthCopyShader;
    
    magnet::GL::Texture3D _data;
    magnet::GL::Texture1D _transferFuncTexture;
    magnet::GL::Texture1D _preintTransferFuncTexture;

    //GTK gui stuff
    std::unique_ptr<Gtk::VBox> _optList;
    std::unique_ptr<Gtk::Entry> _stepSize;
    std::unique_ptr<Gtk::CheckButton> _ditherRay;
    std::unique_ptr<Gtk::CheckButton> _variableDither;
    std::unique_ptr<Gtk::CheckButton> _filterData;
    std::unique_ptr<magnet::gtk::TransferFunction> _transferFunction;
    GLfloat _stepSizeVal;
    Vector _dimensions;
    size_t _frameCounter;
  };
}
