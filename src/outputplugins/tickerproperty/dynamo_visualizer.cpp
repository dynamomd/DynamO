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

#include <magnet/CL/CLGL.hpp>

OPVisualizer::OPVisualizer(const DYNAMO::SimData* tmp, const XMLNode& XML):
  OPTicker(tmp,"Magnet"),
  _CLWindow(NULL),
  _sphereObject(NULL)
{
  operator<<(XML);
}

OPVisualizer::~OPVisualizer()
{}

void 
OPVisualizer::operator<<(const XMLNode& XML)
{}

void
OPVisualizer::initialise()
{
  //Build a window, ready to display it
  _CLWindow = new CLGLWindow(1024, 1024,//height, width
			     0, 0,//initPosition (x,y)
			     "Visualizer"//title
			     );
    
//  //CLWindow.addRenderObj<RTTestWaves>((size_t)1000, 0.0f);
//
  std::vector<RTSpheres::SphereDetails> sphereDetailLevels;

  //Work computer test render
  size_t spheres_rendered = 0;

  size_t stage_spheres = std::min(10ul, Sim->N - spheres_rendered);
  if (stage_spheres)
    {
      sphereDetailLevels.push_back(RTSpheres::SphereDetails(magnet::GL::primatives::Sphere::icosahedron, 2, stage_spheres));
      spheres_rendered += stage_spheres;
    }

  stage_spheres = std::min(1000ul, Sim->N - spheres_rendered);
  if (stage_spheres)
    {
      sphereDetailLevels.push_back(RTSpheres::SphereDetails(magnet::GL::primatives::Sphere::icosahedron, 1, stage_spheres));
      spheres_rendered += stage_spheres;
    }

  stage_spheres = std::min(10000ul, Sim->N - spheres_rendered);
  if (stage_spheres)
    {
      sphereDetailLevels.push_back(RTSpheres::SphereDetails(magnet::GL::primatives::Sphere::icosahedron, 0, stage_spheres));
      spheres_rendered += stage_spheres;
    }

  stage_spheres = std::min(200000ul, Sim->N - spheres_rendered);
  if (stage_spheres)
    {
      sphereDetailLevels.push_back(RTSpheres::SphereDetails(magnet::GL::primatives::Sphere::octahedron, 0, stage_spheres));
      spheres_rendered += stage_spheres;
    }

  stage_spheres = Sim->N - spheres_rendered;
  if (stage_spheres)
    sphereDetailLevels.push_back(RTSpheres::SphereDetails(magnet::GL::primatives::Sphere::tetrahedron, 0, Sim->N - spheres_rendered));

  _sphereObject = &(_CLWindow->addRenderObj<RTSpheres>((size_t)Sim->N, sphereDetailLevels));
  
  CoilMaster::getInstance().addWindow(_CLWindow);

  //Second window test
  CLGLWindow *_CLWindow2 = new CLGLWindow(1024, 1024,//height, width
					  0, 0,//initPosition (x,y)
					  "Visualizer2"//title
					  );

  _CLWindow2->addRenderObj<RTSpheres>((size_t)Sim->N, sphereDetailLevels);
  
  CoilMaster::getInstance().addWindow(_CLWindow2);


  //We must lock coil before doing anything with the window
  {
    magnet::thread::ScopedLock lock(CoilMaster::getInstance()._coilLock);
    
    if (!CoilMaster::getInstance().isRunning())
      M_throw() << "Coil aborted before initial data was loaded into the rendered spheres";

    _lastRenderTime = _CLWindow->getLastFrameTime();
    
    //Place the initial radii into the visualizer
    cl_float4* sphereDataPtr = _sphereObject->writePositionData(_CLWindow->getCLState());
    
    BOOST_FOREACH(const magnet::ClonePtr<Species>& spec, Sim->dynamics.getSpecies())
      {
	double diam = spec->getIntPtr()->hardCoreDiam();
	
	BOOST_FOREACH(unsigned long ID, *(spec->getRange()))
	  {
	    Vector pos = Sim->particleList[ID].getPosition();
	    
	    Sim->dynamics.BCs().applyBC(pos);
	    
	    for (size_t i(0); i < NDIM; ++i)
	      sphereDataPtr[ID].s[i] = pos[i];
	    
	    sphereDataPtr[ID].w = diam * 0.5;
	  }
      }

    //Return it
    _sphereObject->returnPositionData(_CLWindow->getCLState(), sphereDataPtr);
  }

}

void 
OPVisualizer::ticker()
{
  return;

  //Now for the update test
  if (_lastRenderTime == _CLWindow->getLastFrameTime()) return;
  //The screen was redrawn! Lets continue

  //First, lock coil from killing itself
  magnet::thread::ScopedLock lock(CoilMaster::getInstance()._coilLock);
    
  //Check it is actually still running
  if (!CoilMaster::getInstance().isRunning()) return;
  //Still running, so lets do an update
  
  //Now try getting access to the sphere position data
  cl_float4* sphereDataPtr = _sphereObject->writePositionData(_CLWindow->getCLState());

  BOOST_FOREACH(const Particle& part, Sim->particleList)
    {
      Vector pos = part.getPosition();
      
      Sim->dynamics.BCs().applyBC(pos);
      
      for (size_t i(0); i < NDIM; ++i)
	sphereDataPtr[part.getID()].s[i] = pos[i];
    }
    
  //Return it
  _sphereObject->returnPositionData(_CLWindow->getCLState(), sphereDataPtr);

  //Mark when the last update was
  _lastRenderTime = _CLWindow->getLastFrameTime();
}

void 
OPVisualizer::output(xml::XmlStream& XML)
{
}

#endif
