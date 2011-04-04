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
#pragma once
#include <magnet/CL/CLGL.hpp>
#include <magnet/thread/refPtr.hpp>
#include <magnet/thread/taskQueue.hpp>
#include <magnet/GL/viewPort.hpp>
#include <magnet/GL/FBO.hpp>

namespace Gtk {
  class ScrolledWindow;
}

namespace magnet {
  namespace GL {
    class viewPort;
  }
}

class RenderObj
{
public:
  RenderObj(std::string name):
    _name(name),
    _RenderMode(TRIANGLES),
    _renderNormals(false),
    _visible(true)
  {}
  
  ~RenderObj() {}
  
  void accessoryData(const magnet::thread::RefPtr<RenderObj>& console, 
		     const magnet::thread::RefPtr<magnet::thread::TaskQueue>& systemQueue,
		     const magnet::thread::RefPtr<magnet::CL::CLGLState>& CLState,
		     const magnet::thread::RefPtr<magnet::GL::viewPort>& viewPort)
  {
    _console = console;
    _systemQueue = systemQueue;
    _CLState = CLState;
    _viewPort = viewPort;
  }

  virtual void initGTK() {}
  virtual void initOpenGL() {}
  virtual void initOpenCL() {}
  
  virtual void clTick() {}
  virtual void glRender(magnet::GL::FBO& fbo) { glRender(); }
  virtual void glRender() {}
  virtual void interfaceRender() {}

  virtual void initPicking(cl_uint& offset) {}
  virtual void pickingRender() {}
  virtual void finishPicking(cl_uint& offset, const cl_uint val) {}

  virtual void resize(size_t width, size_t height) {}

  virtual void showControls(Gtk::ScrolledWindow* win) {}

  enum RenderModeType 
    {
      POINTS,
      LINES,
      TRIANGLES
    };

  void setRenderMode(RenderModeType rm) { _RenderMode = rm; } 
  
  inline void setDisplayNormals(bool val) { _renderNormals = val; }
  inline void setVisible(bool val) { _visible = val; }
  inline bool isVisible() const { return _visible; }

  inline const std::string& getName() const { return _name; }

  magnet::thread::RefPtr<magnet::thread::TaskQueue> getQueue() { return _systemQueue; }

  virtual void releaseCLGLResources() {}

protected:
  std::string _name;

  RenderModeType _RenderMode;
  bool _renderNormals;
  bool _visible;
  magnet::thread::RefPtr<RenderObj> _console;
  magnet::thread::RefPtr<magnet::thread::TaskQueue> _systemQueue;
  magnet::thread::RefPtr<magnet::CL::CLGLState> _CLState;
  magnet::thread::RefPtr<magnet::GL::viewPort> _viewPort;
};
