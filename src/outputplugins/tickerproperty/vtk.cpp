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
	(boost::lexical_cast<double>(XML.getAttribute("binwidth")),
	 boost::lexical_cast<double>(XML.getAttribute("binwidth")),
	 boost::lexical_cast<double>(XML.getAttribute("binwidth")));
    
    if (XML.isAttributeSet("Snapshots")) snapshots = true;
    if (XML.isAttributeSet("Fields")) fields = true;
    if (XML.isAttributeSet("CollisionStats")) CollisionStats = true;
    
      }
  catch (std::exception& excep)
    {
      M_throw() << "Error while parsing " << name << "options\n"
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
	    M_throw() << "asprintf error in tinkerXYZ";
	  
	  std::ofstream of((std::string("CollStats") + fileName + std::string(".vtu")).c_str());
	  
	  free(fileName);

	  xml::XmlStream XML(of);
	  
	  
	  XML << xml::tag("VTKFile")
	      << xml::attr("type") << "ImageData"
	      << xml::attr("version") << "0.1"
	      << xml::attr("byte_order") << "LittleEndian"
	      << xml::attr("compressor") << "vtkZLibDataCompressor"
	      << xml::tag("ImageData")
	      << xml::attr("WholeExtent");
	  
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    XML << " " << "0 " << nBins[iDim] - 1;
	  
	  XML << xml::attr("Origin");
	  
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    XML << (Sim->aspectRatio[iDim] * (-0.5))
	      / Sim->dynamics.units().unitLength()
		<< " ";
	  
	  XML << xml::attr("Spacing");
	  
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    XML << binWidth[iDim] / Sim->dynamics.units().unitLength() << " ";
	  
	  XML << xml::tag("Piece")
	      << xml::attr("Extent");
	  
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    XML << " " << "0 " << nBins[iDim] - 1;
	  
	  XML << xml::tag("PointData");
	  
	  
	  //////////////////////////HERE BEGINS THE OUTPUT OF THE FIELDS
	  DYNAMO::Line_Breaker lb(6);
	  
	  ////////////SAMPLE COUNTS
	  XML << xml::tag("DataArray")
	      << xml::attr("type") << "Int32"
	      << xml::attr("Name") << "Collisions Per Snapshot"
	      << xml::attr("format") << "ascii"
	      << xml::chardata();
	  
	  BOOST_FOREACH(const unsigned long& val, collCounter)
	    XML << val << lb;
	  
	  XML << "\n" << xml::endtag("DataArray");
	  
	  BOOST_FOREACH(unsigned long& val, collCounter)
	    val = 0;
	  
	  std::vector<size_t> density(nBins[0] * nBins[1] * nBins[2], 0);

	  BOOST_FOREACH(const Particle& part, Sim->particleList)
	    ++density[getCellID(part.getPosition())];
	  
	  XML << xml::tag("DataArray")
	      << xml::attr("type") << "Float32"
	      << xml::attr("Name") << "Density"
	      << xml::attr("format") << "ascii"
	      << xml::chardata();
	  
	  lb.reset();
	  BOOST_FOREACH(const size_t& val, density)
	    XML << (val / binVol) << lb;
	  
	  XML << "\n" << xml::endtag("DataArray");

	  ////////////Postamble
	  XML << xml::endtag("PointData")
	      << xml::tag("CellData")
	      << xml::endtag("CellData")
	      << xml::endtag("Piece")
	      << xml::endtag("ImageData")
	      << xml::endtag("VTKFile");
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
	M_throw() << "Your bin width is too large for the " << iDim 
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
	M_throw() << "asprintf error in tinkerXYZ";
      
      std::ofstream of((std::string("paraview") + fileName + std::string(".vtu")).c_str());
      
      free(fileName);

      xml::XmlStream XML(of);
      
      XML //<< std::scientific
	//This has a minus one due to the digit in front of the decimal
	//An extra one is added if we're rounding
	<< std::setprecision(std::numeric_limits<double>::digits10 - 1)
	<< xml::prolog() << xml::tag("VTKFile")
	<< xml::attr("type") << "UnstructuredGrid"
	<< xml::attr("version") << "0.1"
	<< xml::attr("byte_order") << "LittleEndian"
	<< xml::tag("UnstructuredGrid")
	<< xml::tag("Piece") 
	<< xml::attr("NumberOfPoints") << Sim->N
	<< xml::attr("NumberOfCells") << 0
	<< xml::tag("Points")
	<< xml::tag("DataArray")
	<< xml::attr("type") << "Float32"
      	<< xml::attr("format") << "ascii"
      	<< xml::attr("NumberOfComponents") << "3"
	<< xml::chardata();
      
      BOOST_FOREACH(const Particle& part, Sim->particleList)
	XML << part.getPosition()[0] / Sim->dynamics.units().unitLength() << " "
	    << part.getPosition()[1] / Sim->dynamics.units().unitLength() << " "
	    << part.getPosition()[2] / Sim->dynamics.units().unitLength() << "\n";
      
      XML << xml::endtag("DataArray")
	  << xml::endtag("Points")
	  << xml::tag("Cells") 

	  << xml::tag("DataArray")
	  << xml::attr("type") << "Int32" 
	  << xml::attr("Name") << "connectivity" 
	  << xml::attr("format") << "ascii" 
	  << xml::endtag("DataArray") 

	  << xml::tag("DataArray") 
	  << xml::attr("type") << "Int32" 
	  << xml::attr("Name") << "offsets" 
	  << xml::attr("format") << "ascii" 
	  << xml::endtag("DataArray") 

	  << xml::tag("DataArray") 
	  << xml::attr("type") << "UInt8" 
	  << xml::attr("Name") << "types" 
	  << xml::attr("format") << "ascii" 
	  << xml::endtag("DataArray") 

	  << xml::endtag("Cells")
	  << xml::tag("CellData") << xml::endtag("CellData")
	  << xml::tag("PointData"); 

      //Velocity data    
      XML << xml::tag("DataArray")
	  << xml::attr("type") << "Float32"
	  << xml::attr("Name") << "Velocities"
	  << xml::attr("NumberOfComponents") << "3"
	  << xml::attr("format") << "ascii"
	  << xml::chardata();
    
      BOOST_FOREACH(const Particle& part, Sim->particleList)
	XML << part.getVelocity()[0] / Sim->dynamics.units().unitVelocity() << " "
	    << part.getVelocity()[1] / Sim->dynamics.units().unitVelocity() << " "
	    << part.getVelocity()[2] / Sim->dynamics.units().unitVelocity() << "\n";
    
      XML << xml::endtag("DataArray");

      if (Sim->dynamics.liouvilleanTypeTest<LNOrientation>())
	{
	  //Orientation data
	  XML << xml::tag("DataArray")
	      << xml::attr("type") << "Float32"
	      << xml::attr("Name") << "Orientations"
	      << xml::attr("NumberOfComponents") << "3"
	      << xml::attr("format") << "ascii"
	      << xml::chardata();
    
	  BOOST_FOREACH(const Particle& part, Sim->particleList)
	    {
	      const Vector& tmp = static_cast<const LNOrientation&>
		(Sim->dynamics.getLiouvillean()).getRotData(part).orientation;
	     
	      XML << tmp[0] << " " << tmp[1] << " " << tmp[2] << "\n";
	    }

	  XML << xml::endtag("DataArray");
	}

      XML << xml::endtag("PointData")
	  << xml::endtag("Piece")
	  << xml::endtag("UnstructuredGrid")
	  << xml::endtag("VTKFile")
	;
    }
}

void 
OPVTK::output(xml::XmlStream& XML)
{
  XML << xml::tag("VTK")
      << xml::attr("ImagesTaken") << imageCounter
      << xml::tag("VTKFile")
      << xml::attr("type") << "ImageData"
      << xml::attr("version") << "0.1"
      << xml::attr("byte_order") << "LittleEndian"
      << xml::attr("compressor") << "vtkZLibDataCompressor"
      << xml::tag("ImageData")
      << xml::attr("WholeExtent");
  
  for (size_t iDim(0); iDim < NDIM; ++iDim)
    XML << " " << "0 " << nBins[iDim] - 1;
   
  XML << xml::attr("Origin");

  for (size_t iDim(0); iDim < NDIM; ++iDim)
    XML << (Sim->aspectRatio[iDim] * (-0.5))
      / Sim->dynamics.units().unitLength()
	<< " ";
  
  XML << xml::attr("Spacing");
  
  for (size_t iDim(0); iDim < NDIM; ++iDim)
    XML << binWidth[iDim] / Sim->dynamics.units().unitLength() << " ";
  
  XML << xml::tag("Piece")
      << xml::attr("Extent");
  
  for (size_t iDim(0); iDim < NDIM; ++iDim)
    XML << " " << "0 " << nBins[iDim] - 1;

  XML << xml::tag("PointData");


  //////////////////////////HERE BEGINS THE OUTPUT OF THE FIELDS
  DYNAMO::Line_Breaker lb(6);

  ////////////SAMPLE COUNTS
  XML << xml::tag("DataArray")
      << xml::attr("type") << "Int32"
      << xml::attr("Name") << "Samples per cell"
      << xml::attr("format") << "ascii"
      << xml::chardata();

  for (size_t id(0); id < SampleCounter.size(); ++id)
    XML << SampleCounter[id] << lb;

  XML << "\n" << xml::endtag("DataArray");

  ////////////Momentum field
  lb.reset();

  XML << xml::tag("DataArray")
      << xml::attr("type") << "Float32"
      << xml::attr("Name") << "Avg Particle Momentum"
      << xml::attr("NumberOfComponents") << NDIM   
      << xml::attr("format") << "ascii"
      << xml::chardata();

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

  XML << "\n" << xml::endtag("DataArray");
  

  ////////////Energy
  lb.reset();

  XML << xml::tag("DataArray")
      << xml::attr("type") << "Float32"
      << xml::attr("Name") << "Avg Particle Energy"
      << xml::attr("format") << "ascii"
      << xml::chardata();

  for (size_t id(0); id < SampleCounter.size(); ++id)
    //Nans are not tolerated by paraview
    if (SampleCounter[id])
      XML << mVsquared[id] * 0.5 
	/ (SampleCounter[id] * Sim->dynamics.units().unitEnergy()) << lb;
    else
      XML << 0.0 << lb;

  XML << "\n" << xml::endtag("DataArray");

  ////////////Postamble
  XML << xml::endtag("PointData")
      << xml::tag("CellData")
      << xml::endtag("CellData")
      << xml::endtag("Piece")
      << xml::endtag("ImageData")
      << xml::endtag("VTKFile")
      << xml::endtag("VTK");
}
