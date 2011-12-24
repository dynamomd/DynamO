/*  dynamo:- Event driven molecular dynamics simulator 
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
#pragma once

#include <gtkmm.h>

#include "RenderObj.hpp"
#include <magnet/gtk/transferFunction.hpp>
#include <magnet/GL/texture.hpp>
#include <magnet/GL/shader/volume.hpp>
#include <magnet/GL/objects/cube.hpp>
#include <memory>
#include <tr1/array>

namespace coil {
  class RVolume : public RenderObj
  {
  public:
    RVolume(std::string name): RenderObj(name), _stepSizeVal(0.01) {}
  
    virtual void init(const std::tr1::shared_ptr<magnet::thread::TaskQueue>& systemQueue);
    virtual void deinit();
    virtual void forwardRender(magnet::GL::FBO& fbo,
			       const magnet::GL::Camera& cam,
			       const magnet::GL::Camera& light,
			       RenderMode mode);

    virtual void clTick(const magnet::GL::Camera&) {}

    virtual void showControls(Gtk::ScrolledWindow* win);

    void loadRawFile(std::string filename, size_t width, size_t height, 
		     size_t depth, size_t bytes);

  protected:
    void loadRawFileWorker(std::string filename, std::tr1::array<size_t, 3> dim, 
			   size_t bytes);

    void loadSphereTestPattern();

    void loadData(const std::vector<GLubyte>& inbuffer, 
		  size_t width, size_t height, size_t depth);

    void initGTK();
    void guiUpdate();

    void transferFunctionUpdated();

    magnet::GL::shader::VolumeShader _shader;
    magnet::GL::objects::Cube _cube;
    magnet::GL::FBO _currentDepthFBO;

    magnet::GL::Texture3D _data;
    magnet::GL::Texture1D _transferFuncTexture;
    magnet::GL::Texture1D _preintTransferFuncTexture;

    //GTK gui stuff
    std::auto_ptr<Gtk::VBox> _optList;
    std::auto_ptr<Gtk::Entry> _stepSize;
    std::auto_ptr<Gtk::HScale> _diffusiveLighting;
    std::auto_ptr<Gtk::HScale> _specularLighting;
    std::auto_ptr<Gtk::HScale> _ditherRay;
    std::auto_ptr<magnet::gtk::TransferFunction> _transferFunction;
    GLfloat _stepSizeVal;
  };
}
