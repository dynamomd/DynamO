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
#include <iostream>
#include <coil/coilMaster.hpp>
#include <coil/clWindow.hpp>
#include <coil/RenderObj/DataSet.hpp>
#include <magnet/arg_share.hpp>

int main(int argc, char *argv[])
{
  magnet::ArgShare::getInstance().setArgs(argc, argv);

  coil::CoilRegister _coil;
  magnet::thread::RefPtr<coil::CoilWindow> _CLWindow(new coil::CLGLWindow("Visualizer : ", 1.0));
  
  magnet::thread::RefPtr<coil::RenderObj> simdata(new coil::DataSet("Particle Data"));
  _CLWindow.as<coil::CLGLWindow>().addRenderObj(simdata);

  _coil.getInstance().addWindow(_CLWindow);

  //Simulation loop!
  for(double t(0); ; t += 1)
    {
      if (_CLWindow.as<coil::CLGLWindow>().simupdateTick())
        {
          magnet::thread::ScopedLock 
            lock(_CLWindow.as<coil::CLGLWindow>().getDestroyLock());

          if (!_CLWindow.as<coil::CLGLWindow>().isReady()) continue;
          
          std::ostringstream os;
          os << "t:" << t;        
          _CLWindow.as<coil::CLGLWindow>().setSimStatus1(os.str());
          _CLWindow.as<coil::CLGLWindow>().flagNewData();
        }
    }
}
