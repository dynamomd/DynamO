/*  dynamo:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
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

#pragma once
#include <dynamo/inputplugins/cells/cell.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/chain.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/device/back_inserter.hpp>

#include <boost/lexical_cast.hpp>

namespace dynamo {
  struct CUFile: public UCell
  {
    CUFile(Vector  ndimensions, std::string fn, UCell* nextCell):
      UCell(nextCell),
      dimensions(ndimensions),
      fileName(fn)
    {}

    Vector  dimensions;
    std::string fileName;
    std::vector<Vector  > particleCache;
 
    virtual void initialise() 
    {
      { 
	uc->initialise();

	magnet::xml::Document doc;

	namespace io = boost::iostreams;
    
	if (!boost::filesystem::exists(fileName))
	  M_throw() << "Could not find the XML file named " << fileName
		    << "\nPlease check the file exists.";
	{ //This scopes out the file objects
      
	  //We use the boost iostreams library to load the file into a
	  //string which may be compressed.
      
	  //We make our filtering iostream
	  io::filtering_istream inputFile;
      
	  //Now check if we should add a decompressor filter
	  if (std::string(fileName.end()-8, fileName.end()) == ".xml.bz2")
	    inputFile.push(io::bzip2_decompressor());
	  else if (!(std::string(fileName.end()-4, fileName.end()) == ".xml"))
	    M_throw() << "Unrecognized extension for xml file";

	  //Finally, add the file as a source
	  inputFile.push(io::file_source(fileName));
	  
	  //Force the copy to occur
	  io::copy(inputFile, io::back_inserter(doc.getStoredXMLData()));
	}

	doc.parseData();

	magnet::xml::Node xSubNode = doc.getNode("dynamoconfig").getNode("ParticleData");
      
	if (xSubNode.hasAttribute("AttachedBinary")
	    && (std::toupper(xSubNode.getAttribute("AttachedBinary")[0]) == 'Y'))
	  M_throw() << "This packer only works on XML config files without binary data,"
		    << " please convert to plain xml using \"dynamod --text\"";
      
	for (magnet::xml::Node node = xSubNode.fastGetNode("Pt");
	     node.valid(); ++node)
	  {
	    Vector  posVector;
	    posVector << node.getNode("P");
	    particleCache.push_back(posVector);
	  }
      }

      //The file has been loaded now, just clean up the positions
      Vector  centreOPoints(0,0,0);
      BOOST_FOREACH(const Vector & vec, particleCache)
	centreOPoints += vec;

      centreOPoints /= particleCache.size();

      BOOST_FOREACH(Vector & vec, particleCache)
	vec -= centreOPoints;

      BOOST_FOREACH(Vector & vec, particleCache)
	for (size_t iDim(0); iDim < NDIM; ++iDim)
	  vec[iDim] *= dimensions[iDim];
    }

    virtual std::vector<Vector  > placeObjects(const Vector & centre)
    {
      std::vector<Vector  > retval;

      BOOST_FOREACH(const Vector & position, particleCache)
	BOOST_FOREACH(const Vector & vec, uc->placeObjects(position + centre))
	retval.push_back(vec);

      return retval;
    }
  };
}
