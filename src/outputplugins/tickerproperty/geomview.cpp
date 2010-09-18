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

#include "geomview.hpp"
#include <fstream>
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_simdata.hpp"
#include "../../base/is_colormap.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"
#include "../../dynamics/interactions/squarebond.hpp"
#include "../../dynamics/ranges/2RList.hpp"
#include "../../dynamics/liouvillean/OrientationL.hpp"

OPGeomview::OPGeomview(const DYNAMO::SimData* tmp, const XMLNode&):
  OPTicker(tmp,"Geomview"),
  frameCount(0)
{}

void 
OPGeomview::ticker()
{
  printImage();
}

void
OPGeomview::printImage()
{
  char *fileName;

  //Dont let this fill up your hard drive!
  if (frameCount > 1000)
    return;

  if ( asprintf(&fileName, "geomview.frame%05d.list", frameCount++) < 0)
    M_throw() << "asprintf error in geomview";
  
  std::ofstream of(fileName);
  
  if (!of.is_open())
    M_throw() << "Could not open geomview file for writing";

  of << "{LIST\n";


  Vector  corner(0,0,0);
  for (size_t iDim(0); iDim < NDIM; ++iDim)
    corner[iDim] -= 0.5 * Sim->aspectRatio[iDim];

  /*
  for (size_t iDim(0); iDim < NDIM; ++iDim)
    {
      Vector  point(corner);

      point[iDim] = 0.5 * Sim->aspectRatio[iDim];

      of << "{VECT 1 2 1 \n 2 \n 1 \n " << corner[0] << " " << corner[1]
	 << " " << corner[2] << "\n" << point[0]  << " " << point[1] << " " 
	 << point[2] << " \n 1.0 1.0 1.0 1.0 }\n";
    }
  */

  DYNAMO::ColorMap<Iflt> colmap(0, Sim->dynamics.getSpecies().size());
  DYNAMO::RGB tmpCol(0,0,0);

  unsigned int i = 0;

  BOOST_FOREACH(const ClonePtr<Species>& spec, Sim->dynamics.getSpecies())
    {      

      of << "{LIST\n";

      if (dynamic_cast<const LNOrientation*>(&Sim->dynamics.getLiouvillean()) != NULL)
	BOOST_FOREACH(unsigned long ID, *spec->getRange())
	  {
	    const Particle& part = Sim->particleList[ID];
	    Vector  pos = part.getPosition();
	    Sim->dynamics.BCs().applyBC(pos);
	   
	    const LNOrientation::rotData& 
	      rdat(static_cast<const LNOrientation&>
		   (Sim->dynamics.getLiouvillean()).getRotData(part));	    
 
	    tmpCol = colmap.getColor(i + Sim->dynamics.getInteraction
				     (Sim->particleList[ID], 
				      Sim->particleList[ID])->getColourFraction
				     (Sim->particleList[ID]));
	    
	    Vector  point1 = pos - 0.5 * spec->getIntPtr()->maxIntDist() * rdat.orientation;
	    Vector  point2 = pos + 0.5 * spec->getIntPtr()->maxIntDist() * rdat.orientation;

	    of << "{VECT 1 2 1 \n 2 \n 1 \n " << point1[0] << " " << point1[1] 
	       << " " << point1[2] << "\n" << point2[0]  << " " << point2[1] << " " 
	       << point2[2] << " \n " << tmpCol.R << " " <<  tmpCol.G << " " << tmpCol.B << " 1.0 }\n";
	  }
      else
	BOOST_FOREACH(unsigned long ID, *spec->getRange())
	  {
	    const Particle& part = Sim->particleList[ID];
	    Vector  pos = part.getPosition();
	    Sim->dynamics.BCs().applyBC(pos);
	    
	    tmpCol = colmap.getColor(i + Sim->dynamics.getInteraction
				   (Sim->particleList[ID], 
				    Sim->particleList[ID])->getColourFraction
				     (Sim->particleList[ID]));
	    
	    of << "appearance {\n"
	       << "material {\ndiffuse "<< tmpCol.R  <<" "<< tmpCol.G <<" " 
	       << tmpCol.B <<" }\n" << "}\n" 
	       << "SPHERE " << 
	      Sim->dynamics.getInteraction(part, part)->hardCoreDiam()/2.0 
	       << " " << pos[0] << " " << pos[1] << " " << pos[2] << " \n";
	  }
      of << "\n}\n";
      ++i;
    }
  
  //Now add the bonds
  typedef std::pair<const unsigned long, std::list<unsigned long> > mypair;

  BOOST_FOREACH(const ClonePtr<Interaction> intPtr, Sim->dynamics.getInteractions())
    if (dynamic_cast<const ISquareBond *>(intPtr.get_ptr()) != NULL)
      if (dynamic_cast<const C2RList*>(intPtr->getRange().get_ptr()) != NULL)
	BOOST_FOREACH(const mypair& mp, dynamic_cast<const C2RList*>(intPtr->getRange().get_ptr())->getPairMap())
	  BOOST_FOREACH(const unsigned int ID2, mp.second)
	  {
	    
	    Vector  pos = Sim->particleList[mp.first].getPosition();
	    Vector  rij = Sim->particleList[ID2].getPosition() - pos;
	    Sim->dynamics.BCs().applyBC(pos);
	    Sim->dynamics.BCs().applyBC(rij);

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
