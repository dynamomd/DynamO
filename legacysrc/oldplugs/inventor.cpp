/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "inventor.hpp"
#include "../error.hpp"

#include <Inventor/Xt/viewers/SoXtExaminerViewer.h>
#include <Inventor/fields/SoSFFloat.h>
#include <Inventor/fields/SoMFColor.h>
#include <Inventor/fields/SoSFVec3f.h>
#include <Inventor/fields/SoSFBool.h>
#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/nodes/SoCube.h>
#include <Inventor/nodes/SoCone.h>
#include <Inventor/actions/SoWriteAction.h>

SoSeparator *COPInventor::makeMoleculeScene() 
{
  SoSeparator *root = new SoSeparator();  
  {
    SoMaterial *molMat = new SoMaterial();
    molMat->diffuseColor.setValue(1.0, 0.0, 0.0);
    root->addChild(molMat);
  }
  
  for (std::vector<CParticle>::const_iterator iPtr = particleList.begin ();
       iPtr  != particleList.end (); iPtr++)
    {
      //Prep the vectors for plotting
      Vector  pos = iPtr->getPosition(), vel = iPtr->getVelocity();
      CIntEventClass collclass = dynamics->findCollClass(iPtr->getType(), iPtr->getType());      
      dynamics->setPBC(pos);
      vel = vel.unitVector() * collclass.diameter * 1.5;
      vel += pos;
    
      {
	SoSeparator *molecule = new SoSeparator();
	{
	  SoTransform *molTrans = new SoTransform();
	  molTrans->translation.setValue(pos[0], pos[1], pos[2]);
	  molecule->addChild(molTrans);
	  
	  SoSphere *molSph = new SoSphere();
	  molSph->radius = collclass.diameter / 2.0;
	  molecule->addChild(molSph);
	}
	root->addChild(molecule);
      }
    }

  return root;
}


COPInventor::COPInventor(const std::vector<CParticle> &pList, const CDynamics * const dyn):
  COutputPlugin(pList,dyn),
  frameCount(0)
{
  std::cout << "COPInventor: Loaded\n";

  // Init the system
  // Initialize Inventor and Xt  
  myWindow = SoXt::init("DYNAMO");

  if (myWindow == NULL)
    throw CException() << "Could not initialise Inventor and SoXT";

  //To view the current scene in a window (live as it were)
  //viewScene();  

}

COPInventor::~COPInventor()
{
  viewScene();  
  std::cout << "COPInventor: Unloaded\n";
}

void 
COPInventor::viewScene()
{
  std::cout << "COPInventor: Rendering scene\n";
  //Opens the current scene in a window, and freezes the simulation permanently
  SoSeparator *root = makeMoleculeScene();

  SoXtExaminerViewer *myViewer = new SoXtExaminerViewer(myWindow);
  myViewer->setSceneGraph(root);
  myViewer->setTitle("DYNAMO");
  myViewer->show();

  SoXt::show(myWindow);
  SoXt::mainLoop();
}

void 
COPInventor::collisionUpdate(const CIntEvent &collision, const CIntEventData &preColl)
{
  //This prints a binary inventor file of the current collision
  /*  SoWriteAction wa;
  SoSeparator *root = makeMoleculeScene();

  SoOutput *out = wa.getOutput();

  char *fileName;
  asprintf(&fileName,"collision%d.iv", frameCount++);
  
  if (out->openFile(fileName) !=NULL)
    {
      out->setBinary(TRUE);
      wa.apply(root);
      }*/
}
