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

#include <coil/clWindow.hpp>
#include <coil/RenderObj/Function.hpp>
#include <coil/RenderObj/Spheres.hpp>
#include <magnet/arg_share.hpp>

int
main(int argc, char *argv[])
{
  magnet::ArgShare::getInstance().setArgs(argc, argv); //Share the command line args with the library

  double tickTime = 0.5; //This is a variable you want the visualizer
			 //to set. It's typically how often your
			 //simulation tries to update the
			 //visualization.

  //Build a window, ready to display it
  magnet::thread::RefPtr<CoilWindow> Window = new CLGLWindow(800, 600, //height, width
							     0, 0, //initPosition (x,y)
							     "Visualizer Test", //Window name
							     tickTime);


//  //////////////////////////////VISUALIZING FUNCTIONS///////////////////////////////
//  Vector axis1 = Vector(1,-1,0);
//  axis1 /= axis1.nrm();
//  Vector orthvec = Vector(-1,0,-0.3);
//  orthvec /= orthvec.nrm();
//  Vector axis2 = (orthvec ^ axis1);
//  axis2 /= axis2.nrm();
//  Vector axis3 = axis2 ^ axis1;
//  axis3 /= axis3.nrm();
//  axis3 *= 0.1f;
//
//  std::string _function = 
//    "const float decayrate = 2.5f;\n"
//    "const float invWaveLength = 40.0f;\n"
//    "const float freq = -4;\n"
//    "float r = native_sqrt(pos.x * pos.x + pos.y * pos.y);\n"
//    "f = native_exp( - decayrate * r) * native_sin(invWaveLength * r + freq * t);\n";
//
//  std::string _normalCalc = 
//    "float dfodr = native_exp(- decayrate * r)\n"
//    "  * (invWaveLength * native_cos(r * invWaveLength + freq * t)\n"
//    "  + decayrate * native_sin(r * invWaveLength + freq * t));\n"
//    "normal = normalize((float3)(-dfodr * pos.x / r, -dfodr * pos.y / r,1));\n"
//    "if (r==0) normal = (float3)(0,0,1);";
//
//  std::string _colorCalc = "colors[0] = (uchar4)(255*f,255,255,255);";
//
//  Window.as<CLGLWindow>()
//    .addRenderObj(new RFunction((size_t)100,
//				Vector(-1, 0.7, -1),
//				axis1, axis2, axis3, //Axis of the function, x,y,z
//				-0.7, -0.7,//Start point of the functions evaluation (x,y)
//				1.4, 1.4,//Range of the function to evaluate (xrange,yrange
//				false, //Render a set of Axis as well?
//				false, //Is the shape static, i.e. is there no time dependence
//                              "Test Function",
//				_function,
//				_normalCalc,
//				_colorCalc
//				));
  
  //////////////////////////Visualizing Spheres///////////////////////////////////////////////////
  size_t N = 10;

  magnet::thread::RefPtr<RenderObj> Spheres = new RTSpheres(N, "Spheres");
  Window.as<CLGLWindow>().addRenderObj(Spheres);


  std::vector<cl_float4> particleData(N);
  std::vector<cl_uchar4> particleColorData(N);

  float diam = 1.0f/N;
  
  for (size_t i(0); i < N; ++i)
    {
      float pos[3] = {i/float(N), 0, 0};

      for (size_t dim(0); dim < 3; ++dim)
	particleData[i].s[dim] = pos[dim]; //position
      
      particleData[i].s[3] = diam * 0.5; //radius
    }


  //////////////////////////Finished adding objects to draw///////////////////////////////////////
  CoilMaster::getInstance().addWindow(Window);

  //////////////////////////Main loop of the "simulation"/////////////////////////////////////////
  for(;;) {//Our infinite simulation loop
    {

      if (Window.as<CLGLWindow>().simupdateTick())
	{
	  particleData[0].s[1] += 0.0001;
	  if (particleData[0].s[1] > 1.0) particleData[0].s[1] = 0;
	  
	  particleData[1].s[1] += 0.001;
	  if (particleData[1].s[1] > 1.0) particleData[1].s[1] = 0;
	  
	  particleData[2].s[1] -= 0.001;
	  if (particleData[2].s[1] < 0.0) particleData[2].s[1] = 1;
	  
	  Window.as<CLGLWindow>().getCLState().getCommandQueue().enqueueWriteBuffer
	    (static_cast<RTSpheres&>(*Spheres).getSphereDataBuffer(),
	     false, 0, N * sizeof(cl_float4), &particleData[0]);
	}
    }
  }
}
