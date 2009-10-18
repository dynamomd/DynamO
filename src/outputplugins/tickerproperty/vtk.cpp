/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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
#include "../../extcode/xmlwriter.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_exception.hpp"
#include "../../base/is_simdata.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"
#include "../../base/is_stream_op.hpp"

OPVTK::OPVTK(const DYNAMO::SimData* tmp, const XMLNode& XML):
  OPTicker(tmp,"VTK"),
  binWidth(1,1,1),
  imageCounter(0)
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
      }
  catch (std::exception& excep)
    {
      D_throw() << "Error while parsing " << name << "options\n"
		<< excep.what();
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
  
  BOOST_FOREACH(const CParticle & Part, Sim->vParticleList)
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
