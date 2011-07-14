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
#include <boost/lexical_cast.hpp>

struct CUFile: public CUCell
{
  CUFile(Vector  ndimensions, std::string fn, CUCell* nextCell):
    CUCell(nextCell),
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
      //Open the file for XML parsing
      magnet::xml::Document doc(fileName);  

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
