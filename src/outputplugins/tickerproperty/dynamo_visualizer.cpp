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

#ifdef DYNAMO_visualizer

#include <boost/foreach.hpp>

#include <algorithm>

#include "../../dynamics/include.hpp"
#include "dynamo_visualizer.hpp"

#include <coil/clWindow.hpp>
#include <coil/RenderObj/TestWaves.hpp>
#include <coil/RenderObj/Spheres.hpp>

OPVisualizer::OPVisualizer(const DYNAMO::SimData* tmp, const XMLNode& XML):
  OPTicker(tmp,"Magnet"),
  _CLWindow(NULL),
  _sphereObject(NULL)
{
  operator<<(XML);
}

OPVisualizer::~OPVisualizer()
{
  if (_CLWindow != NULL)
    {
      CoilMaster::getInstance().shutdownRenderer();
      CoilMaster::getInstance().waitForRendererShutdown();
    }
}

void 
OPVisualizer::operator<<(const XMLNode& XML)
{
//  try {
//    if (XML.isAttributeSet("binwidth"))
//      binWidth = Vector 
//	(boost::lexical_cast<Iflt>(XML.getAttribute("binwidth")),
//	 boost::lexical_cast<Iflt>(XML.getAttribute("binwidth")),
//	 boost::lexical_cast<Iflt>(XML.getAttribute("binwidth")));
//    
//    if (XML.isAttributeSet("Snapshots")) snapshots = true;
//    if (XML.isAttributeSet("Fields")) fields = true;
//    if (XML.isAttributeSet("CollisionStats")) CollisionStats = true;
//    
//      }
//  catch (std::exception& excep)
//    {
//      M_throw() << "Error while parsing " << name << "options\n"
//		<< excep.what();
//    }
}

void
OPVisualizer::initialise()
{
  std::vector<cl::Platform> platforms;
  cl::Platform::get(&platforms);
  
  cl::Platform clplatform = platforms[0];
    
  int argc = 1;
  const char* argv = "dynarun";

  CoilMaster::getInstance(argc, (char**)&argv);
    
  _CLWindow = new CLGLWindow(1024, 1024,//height, width
			     0, 0,//initPosition (x,y)
			     "GLCLWindow",//title
			     clplatform);
    
//  //CLWindow.addRenderObj<RTTestWaves>((size_t)1000, 0.0f);
//
  std::vector<RTSpheres::SphereDetails> sphereDetailLevels;

  //Work computer test render
  size_t spheres_rendered = 0;

  size_t stage_spheres = std::min(10ul, Sim->N - spheres_rendered);
  if (stage_spheres)
    {
      sphereDetailLevels.push_back(RTSpheres::SphereDetails(Sphere::icosahedron, 2, stage_spheres));
      spheres_rendered += stage_spheres;
    }

  stage_spheres = std::min(1000ul, Sim->N - spheres_rendered);
  if (stage_spheres)
    {
      sphereDetailLevels.push_back(RTSpheres::SphereDetails(Sphere::icosahedron, 1, stage_spheres));
      spheres_rendered += stage_spheres;
    }

  stage_spheres = std::min(10000ul, Sim->N - spheres_rendered);
  if (stage_spheres)
    {
      sphereDetailLevels.push_back(RTSpheres::SphereDetails(Sphere::icosahedron, 0, stage_spheres));
      spheres_rendered += stage_spheres;
    }

  stage_spheres = std::min(200000ul, Sim->N - spheres_rendered);
  if (stage_spheres)
    {
      sphereDetailLevels.push_back(RTSpheres::SphereDetails(Sphere::octahedron, 0, stage_spheres));
      spheres_rendered += stage_spheres;
    }

  stage_spheres = Sim->N - spheres_rendered;
  if (stage_spheres)
    sphereDetailLevels.push_back(RTSpheres::SphereDetails(Sphere::tetrahedron, 0, Sim->N - spheres_rendered));

  _sphereObject = &(_CLWindow->addRenderObj<RTSpheres>((size_t)Sim->N, sphereDetailLevels));
  
  CoilMaster::getInstance().addWindow(_CLWindow);

  //Start the render thread
  CoilMaster::getInstance().bootRenderer();

  _lastRenderTime = _CLWindow->getLastFrameTime();

  //Place the initial radii into the visualizer

  cl_float4* sphereDataPtr = _sphereObject->writePositionData(_CLWindow->getCommandQueue());

  BOOST_FOREACH(const ClonePtr<Species>& spec, Sim->dynamics.getSpecies())
    {
      Iflt diam = spec->getIntPtr()->hardCoreDiam();

      BOOST_FOREACH(unsigned long ID, *(spec->getRange()))
	{
	  Vector pos = Sim->particleList[ID].getPosition();
	  
	  Sim->dynamics.BCs().applyBC(pos);

	  for (size_t i(0); i < NDIM; ++i)
	    sphereDataPtr[ID].s[i] = pos[i];

	  sphereDataPtr[ID].w = diam * 0.5;
	}
    }
  
  //Now update all the particle data
  //sphereDataPtr[0].w = (edit % 2) ? 0.01 : 0.05;
  
  //Return it
  _sphereObject->returnPositionData(_CLWindow->getCommandQueue(), sphereDataPtr);

}

void 
OPVisualizer::ticker()
{
  //Now for the update test
  if (_lastRenderTime == _CLWindow->getLastFrameTime()) return;

  //The screen was redrawn! Lets continue
  
  //Now try getting access to the sphere position data
  cl_float4* sphereDataPtr = _sphereObject->writePositionData(_CLWindow->getCommandQueue());

  BOOST_FOREACH(const Particle& part, Sim->particleList)
    {
      Vector pos = part.getPosition();
      
      Sim->dynamics.BCs().applyBC(pos);
      
      for (size_t i(0); i < NDIM; ++i)
	sphereDataPtr[part.getID()].s[i] = pos[i];
    }
  
  //Now update all the particle data
  //sphereDataPtr[0].w = (edit % 2) ? 0.01 : 0.05;
  
  //Return it
  _sphereObject->returnPositionData(_CLWindow->getCommandQueue(), sphereDataPtr);

  //Mark when the last update was
  _lastRenderTime = _CLWindow->getLastFrameTime();
}

void 
OPVisualizer::output(xml::XmlStream& XML)
{
}

#endif
