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
    imageCount(0)
  {
    operator<<(XML);
  }

  void 
  OPVTK::operator<<(const magnet::xml::Node& XML)
  {}


  void 
  OPVTK::initialise()
  {
    ticker();
  }

  void 
  OPVTK::ticker()
  {
    printImage();
  }

  void
  OPVTK::printImage()
  {
    using namespace magnet::xml;
    XmlStream XML;
    
    XML //<< std::scientific
      //This has a minus one due to the digit in front of the decimal
      //An extra one is added if we're rounding
      << prolog() << tag("VTKFile")
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
	<< endtag("VTKFile")
      ;

    std::ostringstream filename_oss;
    filename_oss << "paraview" << std::setw(5) << std::setfill('0') << i << ".vtu";
    XML.write_file(filename_oss.str());
    ++imageCount;
  }

  void 
  OPVTK::output(magnet::xml::XmlStream&)
  {
    using namespace magnet::xml;
    //Write out the unified file at the end
    XmlStream XML;

    double dt = getTickerTime();

    XML //<< std::scientific
      //This has a minus one due to the digit in front of the decimal
      //An extra one is added if we're rounding
      << std::setprecision(std::numeric_limits<double>::digits10 - 1)
      << prolog()
      << tag("VTKFile")
      << attr("type") << "Collection" 
      << attr("version") << "0.1"
      << attr("byte_order") << "LittleEndian"
      << attr("compressor") << "vtkZLibDataCompressor"
      << tag("Collection");
    for (size_t i(0); i < imageCount; ++i) {
      std::ostringstream filename_oss;
      filename_oss << "paraview" << std::setw(5) << std::setfill('0') << i << ".vtu";
      XML << tag("DataSet")
	  << attr("timestep") << i * dt
	  << attr("group") << ""
	  << attr("part") << "0"
	  << attr("file") << filename.str()
	  << endtag("DataSet");
    }
    XML << endtag("Collection")
	<< endtag("VTKFile");

    XML.write_file("paraview.pvd");
  }
}
