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

#include <dynamo/outputplugins/tickerproperty/vtk.hpp>
#include <dynamo/include.hpp>
#include <dynamo/simulation.hpp>
#include <dynamo/dynamics/dynamics.hpp>
#include <magnet/xmlwriter.hpp>
#include <fstream>
#include <sstream>

namespace dynamo {
  OPVTK::OPVTK(const dynamo::Simulation* tmp, const magnet::xml::Node& XML):
    OPTicker(tmp,"VTK"),
    imageCount(0),
    _fields(true)
  {
    operator<<(XML);
  }

  void 
  OPVTK::operator<<(const magnet::xml::Node& XML) {
    double minBinWidth = 1;
    if (XML.hasAttribute("MinBinWidth"))
      minBinWidth = XML.getAttribute("MinBinWidth").as<size_t>();
    
    _binWidths = Vector{minBinWidth, minBinWidth, minBinWidth};

    if (XML.hasAttribute("NoFields"))
      _fields = false;
  }


  void 
  OPVTK::initialise()
  {
    if (_fields) {
      size_t vecSize(1);
      for (size_t iDim(0); iDim < NDIM; ++iDim)
	{
	  _binWidths[iDim] *= Sim->units.unitLength();
	  
	  if (_binWidths[iDim] > 0.5 * Sim->primaryCellSize[iDim])
	    M_throw() << "Your bin width is too large for dimension " << iDim;
	  
	  _binCounts[iDim] = static_cast<size_t>(Sim->primaryCellSize[iDim] / _binWidths[iDim]);
	  _binWidths[iDim] = Sim->primaryCellSize[iDim] / _binCounts[iDim];
	  
	  vecSize *= _binCounts[iDim];
	}
      	_numberField.resize(vecSize, 0);
      	_massField.resize(vecSize, 0);
	_kineticEnergyField.resize(vecSize, 0);
	_momentumField.resize(vecSize, Vector{0,0,0});

	dout << "Number of bins: ";
	for (size_t iDim(0); iDim < NDIM; ++iDim)
	  dout << _binCounts[iDim] << " ";
	
	dout << std::endl << "Bin width: ";
	for (size_t iDim(0); iDim < NDIM; ++iDim)
	  dout << _binWidths[iDim] / Sim->units.unitLength() << " ";
	
	dout << std::endl;
    }
    
    ticker();
  }

  void 
  OPVTK::ticker()
  {
    using namespace magnet::xml;
    XmlStream XML;
    
    XML << prolog()
	<< tag("VTKFile")
	<< attr("type") << "UnstructuredGrid"
	<< attr("version") << "0.1"
	<< attr("byte_order") << "LittleEndian"
	<< tag("UnstructuredGrid")
	<< tag("Piece") 
	<< attr("NumberOfPoints") << Sim->particles.size()
	<< attr("NumberOfCells") << 0
	<< tag("Points")
	<< tag("DataArray")
	<< attr("type") << "Float32"
	<< attr("format") << "ascii"
	<< attr("NumberOfComponents") << "3"
	<< chardata();
    
    for (const Particle& part: Sim->particles) {
      Vector r = part.getPosition();
      Sim->BCs->applyBC(r);
      r = r / Sim->units.unitLength();
      XML << r[0]  << " " << r[1] << " " << r[2] << "\n";
    }
      
    XML << endtag("DataArray")
	<< endtag("Points")
	<< tag("Cells") 

	<< tag("DataArray")
	<< attr("type") << "Int32" 
	<< attr("Name") << "connectivity" 
	<< attr("format") << "ascii" 
	<< endtag("DataArray") 

	<< tag("DataArray") 
	<< attr("type") << "Int32" 
	<< attr("Name") << "offsets" 
	<< attr("format") << "ascii" 
	<< endtag("DataArray")

	<< tag("DataArray") 
	<< attr("type") << "UInt8" 
	<< attr("Name") << "types" 
	<< attr("format") << "ascii" 
	<< endtag("DataArray") 

	<< endtag("Cells")
	<< tag("CellData")
	<< endtag("CellData")
	<< tag("PointData"); 

    //Velocity data    
    XML << tag("DataArray")
	<< attr("type") << "Float32"
	<< attr("Name") << "Velocities"
	<< attr("NumberOfComponents") << "3"
	<< attr("format") << "ascii"
	<< chardata();
    
    for (const Particle& part: Sim->particles)
      XML << part.getVelocity()[0] / Sim->units.unitVelocity() << " "
	  << part.getVelocity()[1] / Sim->units.unitVelocity() << " "
	  << part.getVelocity()[2] / Sim->units.unitVelocity() << "\n";
    
    XML << endtag("DataArray");
    
    XML << endtag("PointData")
	<< endtag("Piece")
	<< endtag("UnstructuredGrid")
	<< endtag("VTKFile");

    std::ostringstream filename_oss;
    filename_oss << "particles_" << std::setw(5) << std::setfill('0') << imageCount << ".vtu";
    XML.write_file(filename_oss.str());

    if (_fields) {
      std::fill(_numberField.begin(), _numberField.end(), 0);
      std::fill(_massField.begin(), _massField.end(), 0.0);
      std::fill(_momentumField.begin(), _momentumField.end(), Vector{0,0,0});
      std::fill(_kineticEnergyField.begin(), _kineticEnergyField.end(), 0.0);

      for (const Particle& p : Sim->particles) {
	Vector  position = p.getPosition(),
	  velocity = p.getVelocity();
	Sim->BCs->applyBC(position, velocity);
	
	size_t cellID(0);
	size_t factor(1);
	for (size_t iDim(0); iDim < NDIM; ++iDim)
	  {
	    cellID += factor * static_cast<size_t>((position[iDim] + 0.5 * Sim->primaryCellSize[iDim]) / _binWidths[iDim]);
	    factor *= _binCounts[iDim];
	  }

	const double mass = Sim->species[p]->getMass(p.getID());
	++_numberField[cellID];
	_massField[cellID] += mass;
	_momentumField[cellID] += mass * velocity;
	_kineticEnergyField[cellID] += mass * velocity.nrm2() / 2;
      }
      
      XmlStream XML;
      XML << magnet::xml::tag("VTKFile")
	  << magnet::xml::attr("type") << "ImageData"
	  << magnet::xml::attr("version") << "0.1"
	  << magnet::xml::attr("byte_order") << "LittleEndian"
	  << magnet::xml::attr("compressor") << "vtkZLibDataCompressor"
	  << magnet::xml::tag("ImageData")
	  << magnet::xml::attr("WholeExtent");
  
      for (size_t iDim(0); iDim < NDIM; ++iDim)
	XML << " " << "0 " << _binCounts[iDim] - 1;
   
      XML << magnet::xml::attr("Origin");

      for (size_t iDim(0); iDim < NDIM; ++iDim)
	XML << (Sim->primaryCellSize[iDim] * (-0.5))
	  / Sim->units.unitLength()
	    << " ";
  
      XML << magnet::xml::attr("Spacing");
  
      for (size_t iDim(0); iDim < NDIM; ++iDim)
	XML << _binWidths[iDim] / Sim->units.unitLength() << " ";
  
      XML << magnet::xml::tag("Piece")
	  << magnet::xml::attr("Extent");
  
      for (size_t iDim(0); iDim < NDIM; ++iDim)
	XML << " " << "0 " << _binCounts[iDim] - 1;

      double cellVol = 1;
      for (size_t iDim(0); iDim < NDIM; ++iDim)
	cellVol *= _binWidths[iDim];
      
      XML << magnet::xml::tag("PointData");

      ////////////Mass field
      XML << magnet::xml::tag("DataArray")
	  << magnet::xml::attr("type") << "Float32"
	  << magnet::xml::attr("Name") << "Number density"
	  << magnet::xml::attr("NumberOfComponents") << NDIM
	  << magnet::xml::attr("format") << "ascii"
	  << magnet::xml::chardata();

      for (size_t id(0); id < _numberField.size(); ++id)
	for (size_t iDim(0); iDim < NDIM; ++iDim)
	  XML << _numberField[id] * Sim->units.unitVolume() / cellVol << " ";

      XML << "\n" << magnet::xml::endtag("DataArray");

      ////////////Mass field
      XML << magnet::xml::tag("DataArray")
	  << magnet::xml::attr("type") << "Float32"
	  << magnet::xml::attr("Name") << "Mass density"
	  << magnet::xml::attr("NumberOfComponents") << NDIM   
	  << magnet::xml::attr("format") << "ascii"
	  << magnet::xml::chardata();

      for (size_t id(0); id < _massField.size(); ++id)
	for (size_t iDim(0); iDim < NDIM; ++iDim)
	  XML << _massField[id] * Sim->units.unitVolume() / (cellVol * Sim->units.unitMass()) << " ";

      XML << "\n" << magnet::xml::endtag("DataArray");

      ////////////Momentum field
      XML << magnet::xml::tag("DataArray")
	  << magnet::xml::attr("type") << "Float32"
	  << magnet::xml::attr("Name") << "Momentum density"
	  << magnet::xml::attr("NumberOfComponents") << NDIM   
	  << magnet::xml::attr("format") << "ascii"
	  << magnet::xml::chardata();

      for (size_t id(0); id < _momentumField.size(); ++id)
	for (size_t iDim(0); iDim < NDIM; ++iDim)
	  XML << _momentumField[id][iDim]  * Sim->units.unitVolume() / (cellVol * Sim->units.unitMomentum()) << " ";

      XML << "\n" << magnet::xml::endtag("DataArray");
  

      ////////////Energy
      XML << magnet::xml::tag("DataArray")
	  << magnet::xml::attr("type") << "Float32"
	  << magnet::xml::attr("Name") << "Temperature"
	  << magnet::xml::attr("format") << "ascii"
	  << magnet::xml::chardata();

      for (size_t id(0); id < _kineticEnergyField.size(); ++id)
	XML << 2 * _kineticEnergyField[id] / (NDIM * (_numberField[id] + (_numberField[id] == 1)) * Sim->units.unitEnergy()) << " ";

      XML << "\n" << magnet::xml::endtag("DataArray");

      ////////////Postamble
      XML << magnet::xml::endtag("PointData")
	  << magnet::xml::tag("CellData")
	  << magnet::xml::endtag("CellData")
	  << magnet::xml::endtag("Piece")
	  << magnet::xml::endtag("ImageData")
	  << magnet::xml::endtag("VTKFile");

      std::ostringstream filename_oss;
      filename_oss << "fields_" << std::setw(5) << std::setfill('0') << imageCount << ".vti";
      XML.write_file(filename_oss.str());
    }
    ++imageCount;
  }

  namespace {
    void writePVDfile(std::string prefix, std::string ext, size_t imgCount, double dt) {
      using namespace magnet::xml;
      //Write out the unified file at the end
      XmlStream XML;
      XML << prolog()
	  << tag("VTKFile")
	  << attr("type") << "Collection" 
	  << attr("version") << "0.1"
	  << attr("byte_order") << "LittleEndian"
	  << attr("compressor") << "vtkZLibDataCompressor"
	  << tag("Collection");
      for (size_t i(0); i < imgCount; ++i) {
	std::ostringstream filename_oss;
	filename_oss << prefix << "_" << std::setw(5) << std::setfill('0') << i << "." << ext;
	XML << tag("DataSet")
	    << attr("timestep") << i * dt
	    << attr("group") << ""
	    << attr("part") << "0"
	    << attr("file") << filename_oss.str()
	    << endtag("DataSet");
      }
      XML << endtag("Collection")
	  << endtag("VTKFile");

      XML.write_file(prefix+".pvd");
    }
  }
  
  void 
  OPVTK::output(magnet::xml::XmlStream&)
  {
    const double dt = getTickerTime();
    writePVDfile("particles", "vtu", imageCount, dt);
    writePVDfile("fields", "vti", imageCount, dt);
  }
}
