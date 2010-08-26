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

#include "vtk.hpp"
#include <fstream>
#include <boost/foreach.hpp>
#include <iomanip>

#include "../../extcode/xmlwriter.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_exception.hpp"
#include "../../base/is_simdata.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"
#include "../../dynamics/liouvillean/OrientationL.hpp"
#include "../../base/is_stream_op.hpp"
#include "../../dynamics/systems/rescale.hpp"

OPVTK::OPVTK(const DYNAMO::SimData* tmp, const XMLNode& XML):
  OPTicker(tmp,"VTK"),
  binWidth(1,1,1),
  imageCounter(0),
  snapshots(false),
  fields(false),
  CollisionStats(false),
  eventCounter(0),
  collstatsfilecounter(0)
{
  operator<<(XML);
}

void 
OPVTK::operator<<(const XMLNode& XML)
{
  try {
    if (XML.isAttributeSet("binwidth"))
      binWidth = Vector 
	(boost::lexical_cast<Iflt>(XML.getAttribute("binwidth")),
	 boost::lexical_cast<Iflt>(XML.getAttribute("binwidth")),
	 boost::lexical_cast<Iflt>(XML.getAttribute("binwidth")));
    
    if (XML.isAttributeSet("Snapshots")) snapshots = true;
    if (XML.isAttributeSet("Fields")) fields = true;
    if (XML.isAttributeSet("CollisionStats")) CollisionStats = true;
    
      }
  catch (std::exception& excep)
    {
      D_throw() << "Error while parsing " << name << "options\n"
		<< excep.what();
    }
}

void 
OPVTK::eventUpdate(const IntEvent& IEvent, const PairEventData& PDat)
{
  if (CollisionStats)
    {
      ++collCounter[getCellID(PDat.particle1_.getParticle().getPosition())];
      ++collCounter[getCellID(PDat.particle2_.getParticle().getPosition())];

      if (!(++eventCounter % 50000))
	{
	  char *fileName;
	  if ( asprintf(&fileName, "%05ld", ++collstatsfilecounter) < 0)
	    D_throw() << "asprintf error in tinkerXYZ";
	  
	  std::ofstream of((std::string("CollStats") + fileName + std::string(".vtu")).c_str());
	  
	  free(fileName);

	  xmlw::XmlStream XML(of);
	  
	  
	  XML << xmlw::tag("VTKFile")
	      << xmlw::attr("type") << "ImageData"
	      << xmlw::attr("version") << "0.1"
	      << xmlw::attr("byte_order") << "LittleEndian"
	      << xmlw::attr("compressor") << "vtkZLibDataCompressor"
	      << xmlw::tag("ImageData")
	      << xmlw::attr("WholeExtent");
	  
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    XML << " " << "0 " << nBins[iDim] - 1;
	  
	  XML << xmlw::attr("Origin");
	  
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    XML << (Sim->aspectRatio[iDim] * (-0.5))
	      / Sim->dynamics.units().unitLength()
		<< " ";
	  
	  XML << xmlw::attr("Spacing");
	  
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    XML << binWidth[iDim] / Sim->dynamics.units().unitLength() << " ";
	  
	  XML << xmlw::tag("Piece")
	      << xmlw::attr("Extent");
	  
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    XML << " " << "0 " << nBins[iDim] - 1;
	  
	  XML << xmlw::tag("PointData");
	  
	  
	  //////////////////////////HERE BEGINS THE OUTPUT OF THE FIELDS
	  DYNAMO::Line_Breaker lb(6);
	  
	  ////////////SAMPLE COUNTS
	  XML << xmlw::tag("DataArray")
	      << xmlw::attr("type") << "Int32"
	      << xmlw::attr("Name") << "Collisions Per Snapshot"
	      << xmlw::attr("format") << "ascii"
	      << xmlw::chardata();
	  
	  BOOST_FOREACH(const unsigned long& val, collCounter)
	    XML << val << lb;
	  
	  XML << "\n" << xmlw::endtag("DataArray");
	  
	  BOOST_FOREACH(unsigned long& val, collCounter)
	    val = 0;
	  
	  std::vector<size_t> density(nBins[0] * nBins[1] * nBins[2], 0);

	  BOOST_FOREACH(const Particle& part, Sim->particleList)
	    ++density[getCellID(part.getPosition())];
	  
	  XML << xmlw::tag("DataArray")
	      << xmlw::attr("type") << "Float32"
	      << xmlw::attr("Name") << "Density"
	      << xmlw::attr("format") << "ascii"
	      << xmlw::chardata();
	  
	  lb.reset();
	  BOOST_FOREACH(const size_t& val, density)
	    XML << (val / binVol) << lb;
	  
	  XML << "\n" << xmlw::endtag("DataArray");

	  ////////////Postamble
	  XML << xmlw::endtag("PointData")
	      << xmlw::tag("CellData")
	      << xmlw::endtag("CellData")
	      << xmlw::endtag("Piece")
	      << xmlw::endtag("ImageData")
	      << xmlw::endtag("VTKFile");
	}
    }
}

void
OPVTK::initialise()
{
  size_t vecSize = 1;
  
  for (size_t iDim(0); iDim < NDIM; ++iDim)
    {
      binWidth[iDim] *= Sim->dynamics.units().unitLength();
      
      if (binWidth[iDim] > 0.5 * Sim->aspectRatio[iDim])
	D_throw() << "Your bin width is too large for the " << iDim 
		  << " dimension";
      
      nBins[iDim] = static_cast<size_t>
	(Sim->aspectRatio[iDim] / binWidth[iDim]);
      
      //This is just to ensure the bin width fits an integer number of
      //times into the simulation
      binWidth[iDim] = Sim->aspectRatio[iDim] / nBins[iDim];
      
      invBinWidth[iDim] = 1.0 / binWidth[iDim];
      
      vecSize *= nBins[iDim];
    }
  
  binVol = binWidth[0] * binWidth[1] * binWidth[2];
  
  
  if (CollisionStats)
    {
      collCounter.clear();
      collCounter.resize(vecSize, 0);
    }
  
  if (fields)
    {      
      mVsquared.resize(vecSize, 0.0);
      SampleCounter.resize(vecSize, 0);
      Momentum.resize(vecSize, Vector (0,0,0));
      
      std::string tmp("< ");
      
      for (size_t iDim(0); iDim < NDIM; ++iDim)
	tmp += boost::lexical_cast<std::string>(nBins[iDim]) + " ";
      
      I_cout() << "Number of bins " << tmp << ">";
      
      tmp = std::string("< ");
      
      for (size_t iDim(0); iDim < NDIM; ++iDim)
	tmp +=boost::lexical_cast<std::string>
	  (binWidth[iDim]/Sim->dynamics.units().unitLength()) + " ";
      
      I_cout() << "Bin width " << tmp << ">";  
    } 

  ticker();
}

size_t 
OPVTK::getCellID(Vector  pos)
{
  size_t retval(0);
  size_t factor(1);
  
  Sim->dynamics.BCs().applyBC(pos);

  for (size_t iDim(0); iDim < NDIM; ++iDim)
    {
      retval += factor 
	* static_cast<size_t>((pos[iDim] + 0.5 * Sim->aspectRatio[iDim]) 
			      * invBinWidth[iDim]);
      factor *= nBins[iDim];
    }

  return retval;
}

void 
OPVTK::ticker()
{
  ++imageCounter;

  if (fields)
    BOOST_FOREACH(const Particle & Part, Sim->particleList)
      {
	Vector  position = Part.getPosition(),
	  velocity = Part.getVelocity();
	  
	Sim->dynamics.BCs().applyBC(position, velocity);
	  
	size_t id(getCellID(position));
	  
	//Samples
	++SampleCounter[id];
	  
	//Velocity Vectors
	Momentum[id] += velocity 
	  * Sim->dynamics.getSpecies(Part).getMass();
	  
	//Energy Field
	mVsquared[id] += velocity.nrm2()
	  * Sim->dynamics.getSpecies(Part).getMass();
      }

  if (snapshots)
    {
      char *fileName;
      if ( asprintf(&fileName, "%05ld", imageCounter) < 0)
	D_throw() << "asprintf error in tinkerXYZ";
      
      std::ofstream of((std::string("paraview") + fileName + std::string(".vtu")).c_str());
      
      free(fileName);

      xmlw::XmlStream XML(of);
      
      XML //<< std::scientific
	//This has a minus one due to the digit in front of the decimal
	//An extra one is added if we're rounding
	<< std::setprecision(std::numeric_limits<Iflt>::digits10 - 1)
	<< xmlw::prolog() << xmlw::tag("VTKFile")
	<< xmlw::attr("type") << "UnstructuredGrid"
	<< xmlw::attr("version") << "0.1"
	<< xmlw::attr("byte_order") << "LittleEndian"
	<< xmlw::tag("UnstructuredGrid")
	<< xmlw::tag("Piece") 
	<< xmlw::attr("NumberOfPoints") << Sim->N
	<< xmlw::attr("NumberOfCells") << 0
	<< xmlw::tag("Points")
	<< xmlw::tag("DataArray")
	<< xmlw::attr("type") << "Float32"
      	<< xmlw::attr("format") << "ascii"
      	<< xmlw::attr("NumberOfComponents") << "3"
	<< xmlw::chardata();
      
      BOOST_FOREACH(const Particle& part, Sim->particleList)
	XML << part.getPosition()[0] / Sim->dynamics.units().unitLength() << " "
	    << part.getPosition()[1] / Sim->dynamics.units().unitLength() << " "
	    << part.getPosition()[2] / Sim->dynamics.units().unitLength() << "\n";
      
      XML << xmlw::endtag("DataArray")
	  << xmlw::endtag("Points")
	  << xmlw::tag("Cells") 

	  << xmlw::tag("DataArray")
	  << xmlw::attr("type") << "Int32" 
	  << xmlw::attr("Name") << "connectivity" 
	  << xmlw::attr("format") << "ascii" 
	  << xmlw::endtag("DataArray") 

	  << xmlw::tag("DataArray") 
	  << xmlw::attr("type") << "Int32" 
	  << xmlw::attr("Name") << "offsets" 
	  << xmlw::attr("format") << "ascii" 
	  << xmlw::endtag("DataArray") 

	  << xmlw::tag("DataArray") 
	  << xmlw::attr("type") << "UInt8" 
	  << xmlw::attr("Name") << "types" 
	  << xmlw::attr("format") << "ascii" 
	  << xmlw::endtag("DataArray") 

	  << xmlw::endtag("Cells")
	  << xmlw::tag("CellData") << xmlw::endtag("CellData")
	  << xmlw::tag("PointData"); 

      //Velocity data    
      XML << xmlw::tag("DataArray")
	  << xmlw::attr("type") << "Float32"
	  << xmlw::attr("Name") << "Velocities"
	  << xmlw::attr("NumberOfComponents") << "3"
	  << xmlw::attr("format") << "ascii"
	  << xmlw::chardata();
    
      BOOST_FOREACH(const Particle& part, Sim->particleList)
	XML << part.getVelocity()[0] / Sim->dynamics.units().unitVelocity() << " "
	    << part.getVelocity()[1] / Sim->dynamics.units().unitVelocity() << " "
	    << part.getVelocity()[2] / Sim->dynamics.units().unitVelocity() << "\n";
    
      XML << xmlw::endtag("DataArray");

      if (Sim->dynamics.liouvilleanTypeTest<LNOrientation>())
	{
	  //Orientation data
	  XML << xmlw::tag("DataArray")
	      << xmlw::attr("type") << "Float32"
	      << xmlw::attr("Name") << "Orientations"
	      << xmlw::attr("NumberOfComponents") << "3"
	      << xmlw::attr("format") << "ascii"
	      << xmlw::chardata();
    
	  BOOST_FOREACH(const Particle& part, Sim->particleList)
	    {
	      const Vector& tmp = static_cast<const LNOrientation&>
		(Sim->dynamics.getLiouvillean()).getRotData(part).orientation;
	     
	      XML << tmp[0] << " " << tmp[1] << " " << tmp[2] << "\n";
	    }

	  XML << xmlw::endtag("DataArray");
	}

      XML << xmlw::endtag("PointData")
	  << xmlw::endtag("Piece")
	  << xmlw::endtag("UnstructuredGrid")
	  << xmlw::endtag("VTKFile")
	;
    }
}

void 
OPVTK::output(xmlw::XmlStream& XML)
{
  XML << xmlw::tag("VTK")
      << xmlw::attr("ImagesTaken") << imageCounter
      << xmlw::tag("VTKFile")
      << xmlw::attr("type") << "ImageData"
      << xmlw::attr("version") << "0.1"
      << xmlw::attr("byte_order") << "LittleEndian"
      << xmlw::attr("compressor") << "vtkZLibDataCompressor"
      << xmlw::tag("ImageData")
      << xmlw::attr("WholeExtent");
  
  for (size_t iDim(0); iDim < NDIM; ++iDim)
    XML << " " << "0 " << nBins[iDim] - 1;
   
  XML << xmlw::attr("Origin");

  for (size_t iDim(0); iDim < NDIM; ++iDim)
    XML << (Sim->aspectRatio[iDim] * (-0.5))
      / Sim->dynamics.units().unitLength()
	<< " ";
  
  XML << xmlw::attr("Spacing");
  
  for (size_t iDim(0); iDim < NDIM; ++iDim)
    XML << binWidth[iDim] / Sim->dynamics.units().unitLength() << " ";
  
  XML << xmlw::tag("Piece")
      << xmlw::attr("Extent");
  
  for (size_t iDim(0); iDim < NDIM; ++iDim)
    XML << " " << "0 " << nBins[iDim] - 1;

  XML << xmlw::tag("PointData");


  //////////////////////////HERE BEGINS THE OUTPUT OF THE FIELDS
  DYNAMO::Line_Breaker lb(6);

  ////////////SAMPLE COUNTS
  XML << xmlw::tag("DataArray")
      << xmlw::attr("type") << "Int32"
      << xmlw::attr("Name") << "Samples per cell"
      << xmlw::attr("format") << "ascii"
      << xmlw::chardata();

  for (size_t id(0); id < SampleCounter.size(); ++id)
    XML << SampleCounter[id] << lb;

  XML << "\n" << xmlw::endtag("DataArray");

  ////////////Momentum field
  lb.reset();

  XML << xmlw::tag("DataArray")
      << xmlw::attr("type") << "Float32"
      << xmlw::attr("Name") << "Avg Particle Momentum"
      << xmlw::attr("NumberOfComponents") << NDIM   
      << xmlw::attr("format") << "ascii"
      << xmlw::chardata();

  for (size_t id(0); id < Momentum.size(); ++id)
    {
      for (size_t iDim(0); iDim < NDIM; ++iDim)
	//Nans are not tolerated by paraview
	if (SampleCounter[id])	  
	  XML << Momentum[id][iDim] 
	    / (SampleCounter[id] * Sim->dynamics.units().unitMomentum()) 
	      << lb;
	else
	  XML << 0.0 << lb;
    }

  XML << "\n" << xmlw::endtag("DataArray");
  

  ////////////Energy
  lb.reset();

  XML << xmlw::tag("DataArray")
      << xmlw::attr("type") << "Float32"
      << xmlw::attr("Name") << "Avg Particle Energy"
      << xmlw::attr("format") << "ascii"
      << xmlw::chardata();

  for (size_t id(0); id < SampleCounter.size(); ++id)
    //Nans are not tolerated by paraview
    if (SampleCounter[id])
      XML << mVsquared[id] * 0.5 
	/ (SampleCounter[id] * Sim->dynamics.units().unitEnergy()) << lb;
    else
      XML << 0.0 << lb;

  XML << "\n" << xmlw::endtag("DataArray");

  ////////////Postamble
  XML << xmlw::endtag("PointData")
      << xmlw::tag("CellData")
      << xmlw::endtag("CellData")
      << xmlw::endtag("Piece")
      << xmlw::endtag("ImageData")
      << xmlw::endtag("VTKFile")
      << xmlw::endtag("VTK");
}
