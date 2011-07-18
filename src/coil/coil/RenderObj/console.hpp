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
#include <coil/RenderObj/RenderObj.hpp>
#include <magnet/GL/objects/axis.hpp>
#include <magnet/GL/objects/grid.hpp>
#include <tr1/array>
#include <memory>
#include <sstream>
#include <list>

namespace coil {
  class Console: public RenderObj
  {
  public:
    struct end {};

    Console(std::tr1::array<GLfloat, 3> color):
      RenderObj("Console"),
      _consoleTextColor(color)
    {}
    
    template<class T>
    Console& operator<<(const T& value) 
    {
      os << value;
      return *this;
    }

    void resize(size_t width, size_t height);

    virtual void interfaceRender();

    virtual void initOpenGL();
    virtual void initGTK();
    virtual void showControls(Gtk::ScrolledWindow* win);

    virtual void releaseCLGLResources() { _axis.deinit();}

    virtual void glRender();
    
  private:
    void guiUpdate();

    std::ostringstream os;
    typedef std::pair<float, std::string> consoleEntry;
    std::list<consoleEntry> _consoleEntries;    
    int _glutLastTime;
    std::tr1::array<GLfloat, 3> _consoleTextColor;

    magnet::GL::objects::Axis _axis;
    magnet::GL::objects::Grid _grid;

    std::auto_ptr<Gtk::VBox> _optList; 
    std::auto_ptr<Gtk::CheckButton> _showGrid;
    std::auto_ptr<Gtk::CheckButton> _showConsole;
    std::auto_ptr<Gtk::CheckButton> _showAxis;
  };

  template<>
  inline Console& Console::operator<< <Console::end>(const Console::end&)
  {
    _consoleEntries.push_front(consoleEntry(0, os.str()));
    os.str("");
    return *this;
  }
}
