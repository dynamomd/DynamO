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

//A prototype for the function which will update the data in the renderer
void updateRenderer();

//Some shared data between updateRenderer() and main()
//Number of particles
const size_t N = 10;
//The current simulation time
double t = 0;
std::tr1::shared_ptr<coil::CLGLWindow> window;
std::tr1::shared_ptr<coil::DataSet> data;

int main(int argc, char *argv[])
{
  //Initialise your simulation here
  //...
  //...

  //Now to initialise coil.

  //Share the command line arguments with coil (and any other
  //magnet-dependent library)
  magnet::ArgShare::getInstance().setArgs(argc, argv);

  //Register that this thread (your program) wants to use coil, this
  //must only be destroyed when you have finished using coil
  coil::CoilRegister coil;

  //Make our window
  window.reset(new coil::CLGLWindow("Visualizer : ", 1.0));

  //Make our data set and add it to the window, you can do this
  //multiple times for several different data sets
  data.reset(new coil::DataSet("Particle Data", N));
  window->addRenderObj(data);

  //Add our window to the coil library, the window will run in the
  //main coil thread.
  coil.getInstance().addWindow(window);

  //Now add some data (aka Attributes) to our dataset
  data->addAttribute("Positions", coil::Attribute::COORDINATE, 3);
  data->addAttribute("1 Component values", coil::Attribute::INTENSIVE, 1);
  data->addAttribute("2 Component values", coil::Attribute::INTENSIVE, 2);
  data->addAttribute("3 Component values", coil::Attribute::INTENSIVE, 3);
  data->addAttribute("4 Component values", coil::Attribute::INTENSIVE, 4);

  //Inform the window of a function it needs to call to update the
  //simulation data
  window->signal_data_update().connect(&updateRenderer);  
  //And call it once to get the initial data loaded into coil.
  updateRenderer();

  /* Your simulation loop */
  while (true)
    {
      /* Run your simulation timestep here*/
      //Here we just increment the simulation time, very boring
      t += 1;

      /* Now update the visualisation */
      window->simupdateTick(t);
    }
}

/* This function is an example of how to update the data in the
   simulation */
void updateRenderer()
{
  //Here we update the positions
  std::vector<GLfloat>& posdata = (*data)["Positions"].getData();
  std::vector<GLfloat>& data_1 = (*data)["1 Component values"].getData();
  std::vector<GLfloat>& data_2 = (*data)["2 Component values"].getData();
  std::vector<GLfloat>& data_3 = (*data)["3 Component values"].getData();
  std::vector<GLfloat>& data_4 = (*data)["4 Component values"].getData();
  for (size_t i(0); i < N; ++i)
    {
      posdata[3 * i + 0] = std::sin(t * 0.01 + i);
      posdata[3 * i + 1] = std::cos(t * 0.01 + i);
      posdata[3 * i + 2] = i;
      
      data_1[i] = data_2[2 * i + 0] = data_3[3 * i + 0] = data_4[4 * i + 0] = std::sin(t * 0.01 + i);
      data_2[2 * i + 1] = data_3[3 * i + 1] = data_4[4 * i + 1] = std::cos(t * 0.01 + i);
      data_3[3 * i + 2] = data_4[4 * i + 2] = std::cos(t * 0.01 + 13.131 * M_PI * i);
      data_4[4 * i + 3] = std::sin(t * 0.01 + 12304.123 * M_PI * i);
    }
  
  (*data)["Positions"].flagNewData();
  (*data)["1 Component values"].flagNewData();
  (*data)["2 Component values"].flagNewData();
  (*data)["3 Component values"].flagNewData();
  (*data)["4 Component values"].flagNewData();
}
