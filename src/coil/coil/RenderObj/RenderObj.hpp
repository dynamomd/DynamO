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
#include <coil/Maths/Maths.h>
#include <coil/extcode/vector2.hpp>
#include <magnet/GL/viewPort.hpp>
#include <magnet/thread/refPtr.hpp>
#include <magnet/function/delegate.hpp>

namespace Gtk {
  class ScrolledWindow;
}

class RenderObj
{
public:
  RenderObj(std::string name):
    _name(name),
    _RenderMode(TRIANGLES),
    _renderNormals(false),
    _visible(true)
  {
    _controlDelegate = magnet::function::MakeDelegate(this, &RenderObj::showControls);
  }
  
  ~RenderObj() {}
  
  void accessoryData(const magnet::thread::RefPtr<RenderObj>& console) 
  {
    _console = console;
  }

  virtual void initOpenGL() {}
  virtual void initOpenCL(magnet::CL::CLGLState&) {}
  
  virtual void clTick(magnet::CL::CLGLState&, const magnet::GL::viewPort&) {}
  virtual void glRender() {}
  virtual void interfaceRender(const magnet::GL::viewPort&) {}

  virtual void initPicking(magnet::CL::CLGLState& CLState, cl_uint& offset) {}
  virtual void pickingRender() {}
  virtual void finishPicking(magnet::CL::CLGLState& CLState, cl_uint& offset, const cl_uint val) {}

  virtual void resize(size_t width, size_t height) {}

  //!Actual entry point for window showControls
  void callShowControls(Gtk::ScrolledWindow* win) { _controlDelegate(win); }

  virtual void showControls(Gtk::ScrolledWindow* win) {}

  magnet::function::Delegate1<Gtk::ScrolledWindow*, void>&
  getControlDelegate() { return _controlDelegate; }

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

protected:
  std::string _name;

  RenderModeType _RenderMode;
  bool _renderNormals;
  bool _visible;
  magnet::thread::RefPtr<RenderObj> _console;
  
  //When the window wants to see the controls for this Object this delegate is called
  magnet::function::Delegate1<Gtk::ScrolledWindow*, void> _controlDelegate;
};
