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

#include "tinkerxyz.hpp"
#include <fstream>
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_simdata.hpp"
#include "../../base/is_colormap.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"
#include "../../dynamics/interactions/squarebond.hpp"
#include "../../dynamics/ranges/2RList.hpp"
#include "radiusGyration.hpp"
#include "../../dynamics/topology/chain.hpp"

#include "vmd_imd/vmdsock.h"
#include "vmd_imd/imd.h"

#include "../../dynamics/liouvillean/CompressionL.hpp"

static const size_t HEADERSIZE = 8;


OPTinkerXYZ::OPTinkerXYZ(const DYNAMO::SimData* tmp, const XMLNode& XML):
  OPTicker(tmp,"TinkerXYZ"),
  frameCount(0),
  fileOutput(true),
  liveOutput(false),
  blockForVMD(true),
  P1track(false),
  clientsock(NULL),
  sock(NULL),
  port(3333)
{
  operator<<(XML);
}

OPTinkerXYZ::~OPTinkerXYZ()
{
  if (clientsock)
    {
      imd_disconnect(clientsock);
      vmdsock_shutdown(clientsock);
      vmdsock_destroy(clientsock);
      clientsock = NULL;
    }
}

void 
OPTinkerXYZ::ticker()
{
  if (fileOutput) printFileImage();
  if (liveOutput) printLiveImage();
}

void 
OPTinkerXYZ::operator<<(const XMLNode& XML)
{
  try 
    {
      if (XML.isAttributeSet("LiveVMD")) liveOutput = true;
      if (XML.isAttributeSet("File")) fileOutput = true;
      if (XML.isAttributeSet("NoFile")) fileOutput = false;
      if (XML.isAttributeSet("NoBlock")) blockForVMD = false;
      if (XML.isAttributeSet("P1Track")) P1track = true;
      if (XML.isAttributeSet("Port")) 
	port = boost::lexical_cast<int>(XML.getAttribute("Port"));
    }
  catch (std::exception& excep)
    {
      M_throw() << "Error while parsing " << name << "options\n"
		<< excep.what();
    }
}

void 
OPTinkerXYZ::initialise()
{ 
  printFileImage();
  
  if (liveOutput) 
    {
      coords.resize(NDIM * Sim->N + (HEADERSIZE / sizeof(float)));
      fill_header((IMDheader *)&coords[0], IMD_FCOORDS, Sim->N);


      I_cout() << "Setting up incoming socket of VMD";
      vmdsock_init();
      sock = vmdsock_create();
      vmdsock_bind(sock, port);
      vmdsock_listen(sock);
      I_cout() << "Listening for VMD connection on port 3333";

      printLiveImage();
    }
}

void
OPTinkerXYZ::printLiveImage()
{
  if (!clientsock)
    {
      if (blockForVMD) 
	{
	  I_cout() << "Blocking simulation till VMD connects";
	  std::cout.flush();
	}

      do {
	if (vmdsock_selread(sock, blockForVMD ? -1 : 0) > 0)
	  {
	    clientsock = vmdsock_accept(sock);
	    if (imd_handshake(clientsock))
	      clientsock = NULL;
	    else
	      {
		I_cout() << "VMD port active, blocking for a handshake";
		int bytes_avail = vmdsock_selread(clientsock, -1);
		if (bytes_avail != 1)
		  {
		    clientsock = NULL;
		    I_cout() << "VMD handshake failed"
			     << "\nFound " << bytes_avail;
		  }
		else
		  {
		    int length;
		    IMDType shakeType = imd_recv_header(clientsock, &length);
		    if (shakeType != IMD_GO) 
		      {
			I_cout() << "VMD handshake failed"
				 << "\nRecieved a shake of " << shakeType
				 << "\nNot an IMD_GO"
				 << "\nIgnoring, these handshakes seem broken on 32bit";
		      }
		    else
		      I_cout() << "Connected to VMD session";
		  }
	      }
	    std::cout.flush();
	  }
      } while ((!clientsock) && blockForVMD);
    }

  if (clientsock)
    {
      double coeff = 3.4 / Sim->dynamics.units().unitLength();
      
      if (Sim->dynamics.liouvilleanTypeTest<LCompression>())
	coeff /= 1.0 + static_cast<const LCompression&>(Sim->dynamics.getLiouvillean()).getGrowthRate() * Sim->dSysTime;

      Vector offset(0,0,0);
      if (P1track)
	offset = Sim->particleList.front().getPosition();

      for (size_t ID(0); ID < Sim->N; ++ID)
	{
	  Vector pos = Sim->particleList[ID].getPosition() - offset;
	  Sim->dynamics.BCs().applyBC(pos);
	  //The plus two is the header offset
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    coords[ID * NDIM + iDim + (HEADERSIZE / sizeof(float))] = coeff * pos[iDim];
	}
      
      int32 size = HEADERSIZE + 12 * Sim->N;

      if (imd_writen(clientsock, (const char*) &coords[0], size) != size) 
	{
	  clientsock = NULL;
	  I_cout() << "VMD session disconnected";
	}
    }
}

void
OPTinkerXYZ::printFileImage()
{
  char *fileName;
  
  //Dont let this fill up your hard drive!
  if (frameCount > 1000)
    return;
  
  std::vector<OPRGyration::molGyrationDat> gyrationData;

  BOOST_FOREACH(const magnet::ClonePtr<Topology>& plugPtr, Sim->dynamics.getTopology())
    if (dynamic_cast<const CTChain*>(plugPtr.get_ptr()) != NULL)
      BOOST_FOREACH(const magnet::ClonePtr<CRange>& range, static_cast<const CTChain*>(plugPtr.get_ptr())->getMolecules())
	gyrationData.push_back(OPRGyration::getGyrationEigenSystem(range,Sim));	    

  if ( asprintf(&fileName, "tinker.frame%05d.xyz", frameCount) < 0)
    M_throw() << "asprintf error in tinkerXYZ";
  
  std::ofstream of(fileName);
  
  free(fileName);

  if ( asprintf(&fileName, "tinker.frame%05d.r3d", frameCount++) < 0)
    M_throw() << "asprintf error in tinkerXYZ";

  std::ofstream obj_of(fileName);

  free(fileName);
 
  if (!of.is_open())
    M_throw() << "Could not open file for writing";

  if (!obj_of.is_open())
    M_throw() << "Could not open object file for writing";

  of << Sim->N << "\nDYNAMO Tinker TXYZ file\n";

  Vector  tmpVec;
  BOOST_FOREACH (const Particle& part, Sim->particleList)
    {
      tmpVec = part.getPosition();
      Sim->dynamics.BCs().applyBC(tmpVec);
      of << "C ";
      for (size_t iDim = 0; iDim < NDIM; iDim++)
	of << tmpVec[iDim] * 3.4 
	  / Sim->dynamics.units().unitLength() << " ";
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
  Vector  tmpVec2;
  BOOST_FOREACH(const OPRGyration::molGyrationDat& mDat, gyrationData)
    {
      tmpVec = mDat.MassCentre;
      Sim->dynamics.BCs().applyBC(tmpVec);

      tmpVec2 = ((tmpVec/Sim->dynamics.units().unitLength()) + 0.2 * mDat.EigenVec[NDIM-1]) * 3.4;
      tmpVec =  ((tmpVec/Sim->dynamics.units().unitLength()) - 0.2 * mDat.EigenVec[NDIM-1]) * 3.4;

      obj_of << "5\n";
      for (size_t iDim = 0; iDim < NDIM; iDim++)
	obj_of << tmpVec[iDim] << " ";
      obj_of << " 0.05 ";

      for (size_t iDim = 0; iDim < NDIM; iDim++)
	obj_of << tmpVec2[iDim] << " ";
      obj_of << " 0.05 1.0 0.0 0.0\n";
    }

  BOOST_FOREACH(const magnet::ClonePtr<Topology>& plugPtr, Sim->dynamics.getTopology())
    if (dynamic_cast<const CTChain*>(plugPtr.get_ptr()) != NULL)
      BOOST_FOREACH(const magnet::ClonePtr<CRange>& range, static_cast<const CTChain*>(plugPtr.get_ptr())->getMolecules())
	for (CRange::const_iterator iPtr = range->begin() + 1; iPtr != range->end(); ++iPtr)
	  {
	    Vector  pos1 = Sim->particleList[*iPtr].getPosition();
	    Vector  pos2 = Sim->particleList[*(iPtr-1)].getPosition();
	    Vector  rij(pos1);
	    rij -= pos2;
	    
	    Sim->dynamics.BCs().applyBC(pos1);

	    Sim->dynamics.BCs().applyBC(pos2);

	    Sim->dynamics.BCs().applyBC(rij);	    	    

	    //Check theres no periodic wrap around, 1.01 is a fudge factor
	    if ((pos1 - pos2).nrm2() < 1.01 * rij.nrm2())
	      {
		pos1 *= 3.4 / Sim->dynamics.units().unitLength();
		
		pos2 *= 3.4 / Sim->dynamics.units().unitLength();

		obj_of << "5\n";
		for (size_t iDim = 0; iDim < NDIM; iDim++)
		  obj_of << pos1[iDim] << " ";
		obj_of << " 0.05 ";
		
		
		for (size_t iDim = 0; iDim < NDIM; iDim++)
		  obj_of << pos2[iDim] << " ";
		obj_of << " 0.05 1.0 1.0 1.0\n";
	      }
	  }
	    
  obj_of.close();
}
