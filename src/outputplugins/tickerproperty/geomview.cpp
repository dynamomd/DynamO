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

#include "geomview.hpp"
#include <fstream>
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_exception.hpp"
#include "../../base/is_simdata.hpp"
#include "../../base/is_colormap.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"
#include "../../dynamics/interactions/squarebond.hpp"
#include "../../dynamics/ranges/2RList.hpp"

COPGeomview::COPGeomview(const DYNAMO::SimData* tmp, const XMLNode&):
  COPTicker(tmp,"Geomview"),
  frameCount(0)
{
  printImage();
}

void 
COPGeomview::ticker()
{
  printImage();
}

void
COPGeomview::printImage()
{
  char *fileName;

  //Dont let this fill up your hard drive!
  if (frameCount > 1000)
    return;

  Sim->Dynamics.Liouvillean().updateAllParticles();

  if ( asprintf(&fileName, "geomview.frame%05d.list", frameCount++) < 0)
    D_throw() << "asprintf error in geomview";
  
  std::ofstream of(fileName);
  
  if (!of.is_open())
    D_throw() << "Could not open geomview file for writing";

  of << "{LIST\n";
  
  DYNAMO::ColorMap<double> colmap(0, Sim->Dynamics.getSpecies().size());
  DYNAMO::ColorMap<double>::RGB tmpCol(0,0,0);

  unsigned int i = 0;

  BOOST_FOREACH(const CSpecies& spec, Sim->Dynamics.getSpecies())
    {      

      of << "{LIST\n";

      BOOST_FOREACH(unsigned long ID, *spec.getRange())
	{
	  const CParticle& part = Sim->vParticleList[ID];
	  CVector<> pos = part.getPosition();
	  Sim->Dynamics.BCs().setPBC(pos);

	  tmpCol = colmap.getColor(i + Sim->Dynamics.getInteraction
				   (Sim->vParticleList[ID], 
				    Sim->vParticleList[ID])->getColourFraction
				   (Sim->vParticleList[ID]));
	  
	  of << "appearance {\n"
	     << "material {\ndiffuse "<< tmpCol.R  <<" "<< tmpCol.G <<" " 
	     << tmpCol.B <<" }\n" << "}\n" 
	     << "SPHERE " << 
	    Sim->Dynamics.getInteraction(part, part)->hardCoreDiam()/2.0 
	     << " " << pos[0] << " " << pos[1] << " " << pos[2] << " \n";
	}
      of << "\n}\n";
      ++i;
    }
  
  //Now add the bonds
  typedef std::pair<const unsigned long, std::list<unsigned long> > mypair;

  BOOST_FOREACH(const smrtPlugPtr<CInteraction> intPtr, Sim->Dynamics.getInteractions())
    if (dynamic_cast<const CISquareBond *>(intPtr.get_ptr()) != NULL)
      if (dynamic_cast<const C2RList*>(intPtr->getRange().get_ptr()) != NULL)
	BOOST_FOREACH(const mypair& mp, dynamic_cast<const C2RList*>(intPtr->getRange().get_ptr())->getPairMap())
	  BOOST_FOREACH(const unsigned int ID2, mp.second)
	  {
	    
	    CVector<> pos = Sim->vParticleList[mp.first].getPosition();
	    CVector<> rij = Sim->vParticleList[ID2].getPosition() - pos;
	    Sim->Dynamics.BCs().setPBC(pos);
	    Sim->Dynamics.BCs().setPBC(rij);

	    of << "{VECT 1 2 1 \n 2 \n 1 \n " << pos[0] << " " << pos[1] 
	       << " " << pos[2] << "\n" << pos[0] + rij[0] << " " << pos[1] + rij[1] << " " 
	       << pos[2] + rij[2] << " \n 0.0 0.0 1.0 1.0 }\n";
	  }
	 
  //This is for walls
  /*Iflt x = 0.3;
  of << "{POLY\n " << x << " -0.5 -0.5  " << x << " -0.5 0.5 " 
     << x << " 0.5 0.5 " << x << " 0.5 -0.5}";

  x = -0.3;
  of << "{POLY\n " << x << " -0.5 -0.5  " << x << " -0.5 0.5 " 
     << x << " 0.5 0.5 " << x << " 0.5 -0.5}";
  */
  
  
  of << "}\n";
  of.close();
}
