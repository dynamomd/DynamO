/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#include "cluster.hpp"

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
#include <Inventor/nodes/SoCylinder.h>
#include <Inventor/actions/SoWriteAction.h>
#include <Inventor/SbLinear.h>
#include <math.h>

#include <list>

SoSeparator *
COPCluster::makeScene(Iflt length)
{
  SoSeparator *root = new SoSeparator();  
  
  //System color
  {
    SoMaterial *molMat = new SoMaterial();
    molMat->diffuseColor.setValue(1.0, 0.0, 0.0);
    root->addChild(molMat);
  }
  
   
  //Particle points
  for (std::vector<CParticle>::const_iterator iPtr = particleList.begin ();
       iPtr  != particleList.end (); iPtr++)
    {
      //Prep the vectors for plotting
      Vector  pos = iPtr->getPosition(), vel = iPtr->getVelocity();
      CIntEventClass collclass = dynamics->findCollClass(iPtr->getType(), iPtr->getType());      
      dynamics->setPBC(pos);
      
      {
	SoSeparator *molecule = new SoSeparator();
	{
	  SoTransform *molTrans = new SoTransform();
	  molTrans->translation.setValue(pos[0], pos[1], pos[2]);
	  molecule->addChild(molTrans);
	  
	  SoSphere *molSph = new SoSphere();
	  //Quarter the size of the molecules
	  molSph->radius = collclass.diameter / 16.0;
	  molecule->addChild(molSph);
	}
	root->addChild(molecule);
      }
    }
  
  Iflt distance = length * dynamics->unitClass().diameter;
  
  Vector  vec, rij, pos1, pos2;
  //Now draw the links in
  for (std::vector<CParticle>::const_iterator iPtr = particleList.begin(); iPtr != particleList.end(); iPtr++)
    {
      pos1 = iPtr->getPosition();
      dynamics->setPBC(pos1);
      
      std::vector<CParticle>::const_iterator jPtr = iPtr;
      jPtr++;
      
      for (; jPtr != particleList.end(); jPtr++)
	{
	  pos2 = jPtr->getPosition();
	  dynamics->setPBC(pos2);
	  
	  rij = pos1 - pos2;
	  dynamics->setPBC(rij);
	  
	  if (rij.length() < distance)
	    {
	      //ok we have to form a link between iPtr and jPtr particles
	      SoSeparator *link = new SoSeparator();
	      {
		SbVec3f Ipos1(pos1[0], pos1[1], pos1[2]),
		  Ipos2(pos2[0],pos2[1], pos2[2]);
		
		{
		  //Move to the centre of the first particle
		  SoTransform *linkTrans = new SoTransform();
		  dynamics->setPBC(vec);
		  linkTrans->pointAt(Ipos2, Ipos1);
		  link->addChild(linkTrans);
		}      
		
		{
		  //Move to the centre of the link
		  SoTransform *linkTrans = new SoTransform();
		  linkTrans->translation.setValue(0.0, 0.0, -rij.length()/2.0);
		  link->addChild(linkTrans);
		}
		
		{
		  //Rotate about the x axis
		  SoTransform *linkTrans = new SoTransform();
		  linkTrans->rotation.setValue(SbVec3f(1, 0, 0), 1.5707963f);
		  link->addChild(linkTrans);
		}
		
		SoCylinder *linkCyl = new SoCylinder();
		linkCyl->height = rij.length();
		linkCyl->radius = dynamics->unitClass().diameter / 16.0;
		link->addChild(linkCyl);
	      }
	      root->addChild(link);
	    }
	}
    }
  
  //  root->unrefNoDelete();
  return root;
}

void
COPCluster::viewClusters(Iflt length)
{
  // Init the system
  // Initialize Inventor and Xt  
  Widget myWindow;
  
  myWindow = SoXt::init("DYNAMO");
  
  if (myWindow == NULL)
    throw CException() << "Could not initialise Inventor and SoXT";
  
  SoSeparator *root = makeScene(length);
  
  SoXtExaminerViewer *myViewer = new SoXtExaminerViewer(myWindow);
  myViewer->setSceneGraph(root);
  myViewer->setTitle("DYNAMO");
  myViewer->show();
  
  SoXt::show(myWindow);
  SoXt::mainLoop();
}

COPCluster::COPCluster(const std::vector<CParticle> &pList, const CDynamics * const dyn):
  COutputPlugin(pList,dyn)
{
  std::cout << "COPCluster: Loaded\n";


  Widget myWindow;
  myWindow = SoXt::init("DYNAMO");  
  //writeToFile(makeScene(1.001), "L=1.001.inv");
  //writeToFile(makeScene(1.002), "L=1.002.inv");
  //writeToFile(makeScene(1.004), "L=1.004.inv");
  //writeToFile(makeScene(1.008), "L=1.008.inv");
  //writeToFile(makeScene(1.01), "L=1.010.inv");
}

COPCluster::~COPCluster()
{
  std::cout << "COPCluster: Unloaded\n";
}

void
COPCluster::writeToFile(SoSeparator *root, char *filename)
{
  Widget myWindow;
  myWindow = SoXt::init("DYNAMO");

  //This prints a binary inventor file of the current collision
  SoWriteAction wa;
  
  SoOutput *out = wa.getOutput();
  
  if (out->openFile(filename) != NULL)
    {
      out->setBinary(TRUE);
      wa.apply(root);
    }
}


void 
COPCluster::collisionUpdate(const CIntEvent &collision, const CIntEventData &preColl)
{}

void
COPCluster::output(xmlw::XmlStream &XML)
{
  //This does the links vs distance graph

  XML << xmlw::tag("Cluster")
      << xmlw::chardata();

  Vector  rij;
  long linkCount;
  Iflt linkMax = (Iflt) particleList.size();
  linkMax *= linkMax -1;
  linkMax /= 2;
  
  for (Iflt distance = dynamics->unitClass().diameter; distance < 0.3; distance += 0.001)
    {
      linkCount = 0;
      
      for (std::vector<CParticle>::const_iterator iPtr = particleList.begin(); iPtr != particleList.end(); iPtr++)
	{
	  std::vector<CParticle>::const_iterator jPtr = iPtr;
	  jPtr++;
	  for (; jPtr != particleList.end(); jPtr++)
	    {
	      rij = iPtr->getPosition() - jPtr->getPosition();
	      dynamics->setPBC(rij);
	      if (rij.length() < distance)
		linkCount++;
	    }
	}
      XML << distance/(dynamics->unitClass().diameter) << " " << linkCount << "\n";
    }

  XML << xmlw::endtag("Cluster");
  
  for (Iflt distance = 1.0; distance < 1.1; distance += 0.01)
    std::cout << distance << " " << orderParameter(distance) << "\n";
}

void
COPCluster::periodicOutput()
{}

Iflt
COPCluster::orderParameter(Iflt length)
{
  Iflt distance = length * dynamics->unitClass().diameter;
  
  std::list< Vector  > rijlist;
  
  Vector  rij;

  for (std::vector<CParticle>::const_iterator iPtr = particleList.begin(); iPtr != particleList.end(); iPtr++)
    {
      std::vector<CParticle>::const_iterator jPtr = iPtr;
      jPtr++;
      for (; jPtr != particleList.end(); jPtr++)
	{
	  rij = iPtr->getPosition() - jPtr->getPosition();
	  dynamics->setPBC(rij);
	  if (rij.length() < distance)
	    rijlist.push_back(rij.unitVector());
	}
    }
  
  //iterate over all pairs of joins
  Iflt avgcos = 0.0;
  long sampleCount= 0;
  for (std::list<Vector  >::const_iterator iPtr = rijlist.begin(); iPtr != rijlist.end(); iPtr++)
    {
      std::list<Vector  >::const_iterator jPtr = iPtr;
      jPtr++;
      for (; jPtr != rijlist.end(); jPtr++)
	{
	  avgcos += fabs((*iPtr) % (*jPtr));
	  sampleCount++;
	}
    }
  
  return avgcos/((Iflt) sampleCount);
}

