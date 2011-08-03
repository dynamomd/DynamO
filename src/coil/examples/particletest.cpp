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
  size_t N = 10;

  magnet::ArgShare::getInstance().setArgs(argc, argv);
  coil::CoilRegister coil;
  std::tr1::shared_ptr<coil::CLGLWindow> window(new coil::CLGLWindow("Visualizer : ", 1.0));
  std::tr1::shared_ptr<coil::DataSet> data(new coil::DataSet("Particle Data", N));
  window->addRenderObj(data);
  coil.getInstance().addWindow(window);

  data->addAttribute("Positions", coil::Attribute::COORDINATE, 3);
  data->addAttribute("3Radii", coil::Attribute::INTENSIVE, 3);
  data->addAttribute("1Radii", coil::Attribute::INTENSIVE, 1);

  /* Your simulation loop */
  for(double t(0); ; t += 1)
    {
      /* Run your simulation timestep here*/
      

      /* Now update the visualisation */
      if (window->simupdateTick())
        {
          magnet::thread::ScopedLock lock(window->getDestroyLock());
          if (!window->isReady()) continue;
         
	  //Perform your updating of the attributes here
	  
	  //Here we update the positions
	  std::vector<GLfloat>& posdata = (*data)["Positions"].getData();
	  std::vector<GLfloat>& radii_vec = (*data)["3Radii"].getData();
	  for (size_t i(0); i < N; ++i)
	    {
	      radii_vec[3 * i + 0] = radii_vec[3 * i + 1] = radii_vec[3 * i + 2] = 0.5;
	      posdata[3 * i + 0] = std::sin(t * 0.01 + i);
	      posdata[3 * i + 1] = std::cos(t * 0.01 + i);
	      posdata[3 * i + 2] = i;
	    }
	  (*data)["Positions"].flagNewData();
	  (*data)["3Radii"].flagNewData();

	  std::vector<GLfloat>& radii_scalar = (*data)["1Radii"].getData();
	  for (size_t i(0); i < N; ++i)
	    radii_scalar[i] = 0.5;
	  (*data)["1Radii"].flagNewData();

          std::ostringstream os;
          os << "t:" << t;        
          window->setSimStatus1(os.str());
          window->flagNewData();
        }
    }
}
