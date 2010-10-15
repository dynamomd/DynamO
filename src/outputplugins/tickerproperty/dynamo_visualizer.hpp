/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifdef DYNAMO_visualizer

#include "ticker.hpp"

#include <coil/clWindow.hpp>

class CLGLWindow;
class RTSpheres;

class OPVisualizer: public OPTicker
{
 public:
  OPVisualizer(const DYNAMO::SimData*, const XMLNode&);

  ~OPVisualizer();

  virtual OutputPlugin *Clone() const
  { return new OPVisualizer(*this); }

  virtual void initialise();

  virtual void ticker();

  virtual void output(xml::XmlStream&);

  void operator<<(const XMLNode&);
  
 protected:
  void set_simlock(bool nv) { _simrun = nv; }

  void dataBuild();

  magnet::thread::RefPtr<CoilWindow> _CLWindow;
  magnet::thread::RefPtr<RenderObj> _sphereObject;
  int _lastRenderTime;

  cl::Event lastUpdate;
  std::vector<cl_float4> particleData;

  volatile bool _simrun;
};

#endif
