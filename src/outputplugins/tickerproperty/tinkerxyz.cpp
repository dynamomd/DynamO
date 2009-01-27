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

#include "tinkerxyz.hpp"
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
#include "radiusGyration.hpp"
#include "../../dynamics/topology/chain.hpp"

COPTinkerXYZ::COPTinkerXYZ(const DYNAMO::SimData* tmp, const XMLNode&):
  COPTicker(tmp,"TinkerXYZ"),
  frameCount(0)
{}

void 
COPTinkerXYZ::ticker()
{
  printImage();
}

void
COPTinkerXYZ::printImage()
{
  char *fileName;

  //Dont let this fill up your hard drive!
  if (frameCount > 1000)
    return;

  std::vector<COPRGyration::molGyrationDat> gyrationData;

  BOOST_FOREACH(const smrtPlugPtr<CTopology>& plugPtr, Sim->Dynamics.getTopology())
    if (dynamic_cast<const CTChain*>(plugPtr.get_ptr()) != NULL)
      BOOST_FOREACH(const smrtPlugPtr<CRange>& range, static_cast<const CTChain*>(plugPtr.get_ptr())->getMolecules())
	gyrationData.push_back(COPRGyration::getGyrationEigenSystem(range,Sim));	    

  if ( asprintf(&fileName, "tinker.frame%05d.xyz", frameCount) < 0)
    D_throw() << "asprintf error in tinkerXYZ";
  
  std::ofstream of(fileName);
  
  free(fileName);

  if ( asprintf(&fileName, "tinker.frame%05d.r3d", frameCount++) < 0)
    D_throw() << "asprintf error in tinkerXYZ";

  std::ofstream obj_of(fileName);

  free(fileName);
 
  if (!of.is_open())
    D_throw() << "Could not open file for writing";

  if (!obj_of.is_open())
    D_throw() << "Could not open object file for writing";

  of << Sim->lN << "\nDYNAMO Tinker TXYZ file\n";

  CVector<> tmpVec;
  BOOST_FOREACH (const CParticle& part, Sim->vParticleList)
    {
      tmpVec = part.getPosition();
      Sim->Dynamics.BCs().setPBC(tmpVec);
      of << "C ";
      for (int iDim = 0; iDim < NDIM; iDim++)
	of << tmpVec[iDim] * 3.4 
	  / Sim->Dynamics.units().unitLength() << " ";
      of << "\n";
    }


  obj_of << "r3d input script\n"
    "167 139          tiles in x,y                         \n"
    "4 6          computing pixels per tile		   \n"
    "4              alti-aliasing scheme 4; 3x3 -> 2x2     \n"
    "0.00 0.00 0.00 background color		           \n"
    "T              shadows on			           \n"
    "20             Phong power			           \n"
    "1.00           secondary light contribution	   \n"
    "0.10           ambient light contribution	           \n"
    "0.50           specular reflection component	   \n"
    "  0.83         Eye position			   \n"
    "1 0 0          main light source position	           \n"
    "1 0 0 0        global xform matrix		           \n"
    "0 1 0 0					           \n"
    "0 0 1 0					           \n"
    "0 0 0 2.406					   \n"
    "3						           \n"
    "*\n*\n*\n";
  CVector<> tmpVec2;
  BOOST_FOREACH(const COPRGyration::molGyrationDat& mDat, gyrationData)
    {
      tmpVec = mDat.MassCentre;
      Sim->Dynamics.BCs().setPBC(tmpVec);

      tmpVec2 = ((tmpVec/Sim->Dynamics.units().unitLength()) + 0.2 * mDat.EigenVec[NDIM-1]) * 3.4;
      tmpVec =  ((tmpVec/Sim->Dynamics.units().unitLength()) - 0.2 * mDat.EigenVec[NDIM-1]) * 3.4;

      obj_of << "5\n";
      for (int iDim = 0; iDim < NDIM; iDim++)
	obj_of << tmpVec[iDim] << " ";
      obj_of << " 0.05 ";

      for (int iDim = 0; iDim < NDIM; iDim++)
	obj_of << tmpVec2[iDim] << " ";
      obj_of << " 0.05 1.0 0.0 0.0\n";
    }

  BOOST_FOREACH(const smrtPlugPtr<CTopology>& plugPtr, Sim->Dynamics.getTopology())
    if (dynamic_cast<const CTChain*>(plugPtr.get_ptr()) != NULL)
      BOOST_FOREACH(const smrtPlugPtr<CRange>& range, static_cast<const CTChain*>(plugPtr.get_ptr())->getMolecules())
	for (CRange::const_iterator iPtr = range->begin() + 1; iPtr != range->end(); ++iPtr)
	  {
	    CVector<> pos1 = Sim->vParticleList[*iPtr].getPosition();
	    CVector<> pos2 = Sim->vParticleList[*(iPtr-1)].getPosition();
	    CVector<> rij(pos1);
	    rij -= pos2;
	    
	    Sim->Dynamics.BCs().setPBC(pos1);

	    Sim->Dynamics.BCs().setPBC(pos2);

	    Sim->Dynamics.BCs().setPBC(rij);	    	    

	    //Check theres no periodic wrap around, 1.01 is a fudge factor
	    if ((pos1 - pos2).square() < 1.01 * rij.square())
	      {
		pos1 *= 3.4 / Sim->Dynamics.units().unitLength();
		
		pos2 *= 3.4 / Sim->Dynamics.units().unitLength();

		obj_of << "5\n";
		for (int iDim = 0; iDim < NDIM; iDim++)
		  obj_of << pos1[iDim] << " ";
		obj_of << " 0.05 ";
		
		
		for (int iDim = 0; iDim < NDIM; iDim++)
		  obj_of << pos2[iDim] << " ";
		obj_of << " 0.05 1.0 1.0 1.0\n";
	      }
	  }
	    
  obj_of.close();
}
