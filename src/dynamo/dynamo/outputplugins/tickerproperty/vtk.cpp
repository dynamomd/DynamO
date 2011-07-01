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

#include "vtk.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_simdata.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"
#include "../../dynamics/systems/rescale.hpp"
#include <magnet/xmlwriter.hpp>
#include <magnet/xmlreader.hpp>
#include <boost/foreach.hpp>
#include <fstream>
#include <iomanip>

OPVTK::OPVTK(const dynamo::SimData* tmp, const magnet::xml::Node& XML):
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
OPVTK::operator<<(const magnet::xml::Node& XML)
{
  try {
    float bw = 1;
    if (XML.hasAttribute("binwidth"))
      bw = XML.getAttribute("binwidth").as<double>();

    binWidth = Vector(bw, bw, bw);
    
    if (XML.hasAttribute("Snapshots")) snapshots = true;
    if (XML.hasAttribute("Fields")) fields = true;
    if (XML.hasAttribute("CollisionStats")) CollisionStats = true;
    
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
	    M_throw() << "asprintf error in VTK";
	  
	  std::ofstream of((std::string("CollStats") + fileName + std::string(".vtu")).c_str());
	  
	  free(fileName);

	  magnet::xml::XmlStream XML(of);
	  
	  
	  XML << magnet::xml::tag("VTKFile")
	      << magnet::xml::attr("type") << "ImageData"
	      << magnet::xml::attr("version") << "0.1"
	      << magnet::xml::attr("byte_order") << "LittleEndian"
	      << magnet::xml::attr("compressor") << "vtkZLibDataCompressor"
	      << magnet::xml::tag("ImageData")
	      << magnet::xml::attr("WholeExtent");
	  
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    XML << " " << "0 " << nBins[iDim] - 1;
	  
	  XML << magnet::xml::attr("Origin");
	  
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    XML << (Sim->primaryCellSize[iDim] * (-0.5))
	      / Sim->dynamics.units().unitLength()
		<< " ";
	  
	  XML << magnet::xml::attr("Spacing");
	  
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    XML << binWidth[iDim] / Sim->dynamics.units().unitLength() << " ";
	  
	  XML << magnet::xml::tag("Piece")
	      << magnet::xml::attr("Extent");
	  
	  for (size_t iDim(0); iDim < NDIM; ++iDim)
	    XML << " " << "0 " << nBins[iDim] - 1;
	  
	  XML << magnet::xml::tag("PointData");
	  
	  
	  //////////////////////////HERE BEGINS THE OUTPUT OF THE FIELDS
	  ////////////SAMPLE COUNTS
	  XML << magnet::xml::tag("DataArray")
	      << magnet::xml::attr("type") << "Int32"
	      << magnet::xml::attr("Name") << "Collisions Per Snapshot"
	      << magnet::xml::attr("format") << "ascii"
	      << magnet::xml::chardata();
	  
	  BOOST_FOREACH(const unsigned long& val, collCounter)
	    XML << val;
	  
	  XML << "\n" << magnet::xml::endtag("DataArray");
	  
	  BOOST_FOREACH(unsigned long& val, collCounter)
	    val = 0;
	  
	  std::vector<size_t> density(nBins[0] * nBins[1] * nBins[2], 0);

	  BOOST_FOREACH(const Particle& part, Sim->particleList)
	    ++density[getCellID(part.getPosition())];
	  
	  XML << magnet::xml::tag("DataArray")
	      << magnet::xml::attr("type") << "Float32"
	      << magnet::xml::attr("Name") << "Density"
	      << magnet::xml::attr("format") << "ascii"
	      << magnet::xml::chardata();
	  
	  BOOST_FOREACH(const size_t& val, density)
	    XML << (val / binVol);
	  
	  XML << "\n" << magnet::xml::endtag("DataArray");

	  ////////////Postamble
	  XML << magnet::xml::endtag("PointData")
	      << magnet::xml::tag("CellData")
	      << magnet::xml::endtag("CellData")
	      << magnet::xml::endtag("Piece")
	      << magnet::xml::endtag("ImageData")
	      << magnet::xml::endtag("VTKFile");
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
      
      if (binWidth[iDim] > 0.5 * Sim->primaryCellSize[iDim])
	M_throw() << "Your bin width is too large for the " << iDim 
		  << " dimension";
      
      nBins[iDim] = static_cast<size_t>
	(Sim->primaryCellSize[iDim] / binWidth[iDim]);
      
      //This is just to ensure the bin width fits an integer number of
      //times into the simulation
      binWidth[iDim] = Sim->primaryCellSize[iDim] / nBins[iDim];
      
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
      
      dout << "Number of bins " << tmp << ">" << std::endl;
      
      tmp = std::string("< ");
      
      for (size_t iDim(0); iDim < NDIM; ++iDim)
	tmp +=boost::lexical_cast<std::string>
	  (binWidth[iDim]/Sim->dynamics.units().unitLength()) + " ";
      
      dout << "Bin width " << tmp << ">" << std::endl;  
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
	* static_cast<size_t>((pos[iDim] + 0.5 * Sim->primaryCellSize[iDim]) 
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
    {
      BOOST_FOREACH(const Particle & Part, Sim->particleList)
	{
	  Vector  position = Part.getPosition(),
	    velocity = Part.getVelocity();
	  
	  Sim->dynamics.BCs().applyBC(position, velocity);
	  
	  size_t id(getCellID(position));
	  
	  //Samples
	  ++SampleCounter[id];
	  
	  double mass = Sim->dynamics.getSpecies(Part).getMass(Part.getID());

	  //Velocity Vectors
	  Momentum[id] += velocity * mass;
	  
	  //Energy Field
	  mVsquared[id] += velocity.nrm2() * mass;
	}
    }

  if (snapshots)
    {
      char *fileName;
      if ( asprintf(&fileName, "%05ld", imageCounter) < 0)
	M_throw() << "asprintf error in tinkerXYZ";
      
      std::ofstream of((std::string("paraview") + fileName + std::string(".vtu")).c_str());
      
      free(fileName);

      magnet::xml::XmlStream XML(of);
      
      XML //<< std::scientific
	//This has a minus one due to the digit in front of the decimal
	//An extra one is added if we're rounding
	<< std::setprecision(std::numeric_limits<double>::digits10 - 1)
	<< magnet::xml::prolog() << magnet::xml::tag("VTKFile")
	<< magnet::xml::attr("type") << "UnstructuredGrid"
	<< magnet::xml::attr("version") << "0.1"
	<< magnet::xml::attr("byte_order") << "LittleEndian"
	<< magnet::xml::tag("UnstructuredGrid")
	<< magnet::xml::tag("Piece") 
	<< magnet::xml::attr("NumberOfPoints") << Sim->N
	<< magnet::xml::attr("NumberOfCells") << 0
	<< magnet::xml::tag("Points")
	<< magnet::xml::tag("DataArray")
	<< magnet::xml::attr("type") << "Float32"
      	<< magnet::xml::attr("format") << "ascii"
      	<< magnet::xml::attr("NumberOfComponents") << "3"
	<< magnet::xml::chardata();
      
      BOOST_FOREACH(const Particle& part, Sim->particleList)
	XML << part.getPosition()[0] / Sim->dynamics.units().unitLength() << " "
	    << part.getPosition()[1] / Sim->dynamics.units().unitLength() << " "
	    << part.getPosition()[2] / Sim->dynamics.units().unitLength() << "\n";
      
      XML << magnet::xml::endtag("DataArray")
	  << magnet::xml::endtag("Points")
	  << magnet::xml::tag("Cells") 

	  << magnet::xml::tag("DataArray")
	  << magnet::xml::attr("type") << "Int32" 
	  << magnet::xml::attr("Name") << "connectivity" 
	  << magnet::xml::attr("format") << "ascii" 
	  << magnet::xml::endtag("DataArray") 

	  << magnet::xml::tag("DataArray") 
	  << magnet::xml::attr("type") << "Int32" 
	  << magnet::xml::attr("Name") << "offsets" 
	  << magnet::xml::attr("format") << "ascii" 
	  << magnet::xml::endtag("DataArray") 

	  << magnet::xml::tag("DataArray") 
	  << magnet::xml::attr("type") << "UInt8" 
	  << magnet::xml::attr("Name") << "types" 
	  << magnet::xml::attr("format") << "ascii" 
	  << magnet::xml::endtag("DataArray") 

	  << magnet::xml::endtag("Cells")
	  << magnet::xml::tag("CellData") << magnet::xml::endtag("CellData")
	  << magnet::xml::tag("PointData"); 

      //Velocity data    
      XML << magnet::xml::tag("DataArray")
	  << magnet::xml::attr("type") << "Float32"
	  << magnet::xml::attr("Name") << "Velocities"
	  << magnet::xml::attr("NumberOfComponents") << "3"
	  << magnet::xml::attr("format") << "ascii"
	  << magnet::xml::chardata();
    
      BOOST_FOREACH(const Particle& part, Sim->particleList)
	XML << part.getVelocity()[0] / Sim->dynamics.units().unitVelocity() << " "
	    << part.getVelocity()[1] / Sim->dynamics.units().unitVelocity() << " "
	    << part.getVelocity()[2] / Sim->dynamics.units().unitVelocity() << "\n";
    
      XML << magnet::xml::endtag("DataArray");

      XML << magnet::xml::endtag("PointData")
	  << magnet::xml::endtag("Piece")
	  << magnet::xml::endtag("UnstructuredGrid")
	  << magnet::xml::endtag("VTKFile")
	;
    }
}

void 
OPVTK::output(magnet::xml::XmlStream& XML)
{
  XML << magnet::xml::tag("VTK")
      << magnet::xml::attr("ImagesTaken") << imageCounter
      << magnet::xml::tag("VTKFile")
      << magnet::xml::attr("type") << "ImageData"
      << magnet::xml::attr("version") << "0.1"
      << magnet::xml::attr("byte_order") << "LittleEndian"
      << magnet::xml::attr("compressor") << "vtkZLibDataCompressor"
      << magnet::xml::tag("ImageData")
      << magnet::xml::attr("WholeExtent");
  
  for (size_t iDim(0); iDim < NDIM; ++iDim)
    XML << " " << "0 " << nBins[iDim] - 1;
   
  XML << magnet::xml::attr("Origin");

  for (size_t iDim(0); iDim < NDIM; ++iDim)
    XML << (Sim->primaryCellSize[iDim] * (-0.5))
      / Sim->dynamics.units().unitLength()
	<< " ";
  
  XML << magnet::xml::attr("Spacing");
  
  for (size_t iDim(0); iDim < NDIM; ++iDim)
    XML << binWidth[iDim] / Sim->dynamics.units().unitLength() << " ";
  
  XML << magnet::xml::tag("Piece")
      << magnet::xml::attr("Extent");
  
  for (size_t iDim(0); iDim < NDIM; ++iDim)
    XML << " " << "0 " << nBins[iDim] - 1;

  XML << magnet::xml::tag("PointData");


  ////////////SAMPLE COUNTS
  XML << magnet::xml::tag("DataArray")
      << magnet::xml::attr("type") << "Int32"
      << magnet::xml::attr("Name") << "Samples per cell"
      << magnet::xml::attr("format") << "ascii"
      << magnet::xml::chardata();

  for (size_t id(0); id < SampleCounter.size(); ++id)
    XML << SampleCounter[id];

  XML << "\n" << magnet::xml::endtag("DataArray");

  ////////////Momentum field
  XML << magnet::xml::tag("DataArray")
      << magnet::xml::attr("type") << "Float32"
      << magnet::xml::attr("Name") << "Avg Particle Momentum"
      << magnet::xml::attr("NumberOfComponents") << NDIM   
      << magnet::xml::attr("format") << "ascii"
      << magnet::xml::chardata();

  for (size_t id(0); id < Momentum.size(); ++id)
    {
      for (size_t iDim(0); iDim < NDIM; ++iDim)
	//Nans are not tolerated by paraview
	if (SampleCounter[id])	  
	  XML << Momentum[id][iDim] 
	    / (SampleCounter[id] * Sim->dynamics.units().unitMomentum());
	else
	  XML << 0.0;
    }

  XML << "\n" << magnet::xml::endtag("DataArray");
  

  ////////////Energy
  XML << magnet::xml::tag("DataArray")
      << magnet::xml::attr("type") << "Float32"
      << magnet::xml::attr("Name") << "Avg Particle Energy"
      << magnet::xml::attr("format") << "ascii"
      << magnet::xml::chardata();

  for (size_t id(0); id < SampleCounter.size(); ++id)
    //Nans are not tolerated by paraview
    if (SampleCounter[id])
      XML << mVsquared[id] * 0.5 
	/ (SampleCounter[id] * Sim->dynamics.units().unitEnergy());
    else
      XML << 0.0;

  XML << "\n" << magnet::xml::endtag("DataArray");

  ////////////Postamble
  XML << magnet::xml::endtag("PointData")
      << magnet::xml::tag("CellData")
      << magnet::xml::endtag("CellData")
      << magnet::xml::endtag("Piece")
      << magnet::xml::endtag("ImageData")
      << magnet::xml::endtag("VTKFile")
      << magnet::xml::endtag("VTK");
}
