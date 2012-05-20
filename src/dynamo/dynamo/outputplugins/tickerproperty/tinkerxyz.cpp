/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.dynamomd.org
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

#include <dynamo/outputplugins/tickerproperty/tinkerxyz.hpp>
#include <dynamo/outputplugins/tickerproperty/radiusGyration.hpp>
#include <dynamo/dynamics/include.hpp>
#include <dynamo/base/is_simdata.hpp>
#include <dynamo/dynamics/liouvillean/liouvillean.hpp>
#include <dynamo/dynamics/interactions/squarebond.hpp>
#include <dynamo/dynamics/ranges/2RList.hpp>
#include <dynamo/dynamics/topology/chain.hpp>
#include <dynamo/dynamics/liouvillean/CompressionL.hpp>
#include <dynamo/outputplugins/tickerproperty/vmd_imd/vmdsock.h>
#include <dynamo/outputplugins/tickerproperty/vmd_imd/imd.h>
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>
#include <fstream>

namespace dynamo {
  static const size_t HEADERSIZE = 8;

  OPTinkerXYZ::OPTinkerXYZ(const dynamo::SimData* tmp, const magnet::xml::Node& XML):
    OPTicker(tmp,"TinkerXYZ"),
    frameCount(0),
    fileOutput(true),
    liveOutput(false),
    blockForVMD(true),
    max_frame_count(1000),
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
  OPTinkerXYZ::operator<<(const magnet::xml::Node& XML)
  {
    try 
      {
	if (XML.hasAttribute("LiveVMD")) liveOutput = true;
	if (XML.hasAttribute("File")) fileOutput = true;
	if (XML.hasAttribute("NoFile")) fileOutput = false;
	if (XML.hasAttribute("NoBlock")) blockForVMD = false;
	if (XML.hasAttribute("P1Track")) P1track = true;
      
	if (XML.hasAttribute("Port"))
	  port = XML.getAttribute("Port").as<int>();

	if (XML.hasAttribute("MaxFrames"))
	  max_frame_count = XML.getAttribute("MaxFrames").as<size_t>();
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


	dout << "Setting up incoming socket of VMD" << std::endl;
	vmdsock_init();
	sock = vmdsock_create();
	vmdsock_bind(sock, port);
	vmdsock_listen(sock);
	dout << "Listening for VMD connection on port 3333" << std::endl;

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
	    dout << "Blocking simulation till VMD connects" << std::endl;
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
		  dout << "VMD port active, blocking for a handshake" << std::endl;
		  int bytes_avail = vmdsock_selread(clientsock, -1);
		  if (bytes_avail != 1)
		    {
		      clientsock = NULL;
		      dout << "VMD handshake failed"
			   << "\nFound " << bytes_avail << std::endl;
		    }
		  else
		    {
		      int length;
		      IMDType shakeType = imd_recv_header(clientsock, &length);
		      if (shakeType != IMD_GO) 
			{
			  dout << "VMD handshake failed"
			       << "\nRecieved a shake of " << shakeType
			       << "\nNot an IMD_GO"
			       << "\nIgnoring, these handshakes seem broken on 32bit" << std::endl;
			}
		      else
			dout << "Connected to VMD session" << std::endl;
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
	    dout << "VMD session disconnected" << std::endl;
	  }
      }
  }

  void
  OPTinkerXYZ::printFileImage()
  {
    char *fileName;
  
    //Dont let this fill up your hard drive!
    if (frameCount > max_frame_count)
      return;
  
    if ( asprintf(&fileName, "tinker.frame%05d.xyz", frameCount) < 0)
      M_throw() << "asprintf error in tinkerXYZ";
  
    std::ofstream of(fileName);
  
    free(fileName);
 
    if (!of.is_open())
      M_throw() << "Could not open file for writing";

    of << Sim->N << "\ndynamo Tinker TXYZ file, t = " 
       << Sim->dSysTime / Sim->dynamics.units().unitLength() 
       << ", NOTE: All units here have been scaled by 3.4 "
      "(the van-der-Walls radius of Carbon!)\n";

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
  }
}
