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

#include "3d_field.hpp"
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/foreach.hpp>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkImageData.h>
#include <vtkRectilinearGrid.h>
#include <vtkPointData.h>
#include <vtkXMLImageDataWriter.h>
#include "../simulation/particle.hpp"
#include "../dynamics/include.hpp"
#include "../extcode/xmlwriter.hpp"

namespace fs = boost::filesystem;
  
OP3DField::OP3DField(DYNAMO::SimData* tmp):
  OutputPlugin(tmp,"3dField"),
  Density(tmp), 
  Vsquared(tmp),
  SampleCounter(tmp),
  Velocity(tmp),
  imageCounter(0)
{}

OP3DField::~OP3DField()
{}

void 
OP3DField::collisionUpdate(const IntEvent &collision, 
			    const CIntEventData &preColl)
{
  if (!(Sim->eventCount % 100))
    {
      imageCounter++;

      BOOST_FOREACH( const Particle & Part, Sim->particleList)
	{
	  Vector  position = Part.getPosition();

	  Sim->dynamics.BCs().setPBC(position);
	  //Samples
	  SampleCounter[position]++;

	  //Velocity Vectors
	  Velocity[position] += Sim->dynamics.getLabVelocity(Part);
	  
	  //Density field
	  Density[position] += 1.0;

	  //Energy Field
	  Vsquared[position] += (Part.getVelocity()).square();
	}
    }
}

void 
OP3DField::output(xml::XmlStream &XML)
{
  //Create a image
  vtkImageData* image = Density.getVTKImage();

  //Add the density
  image->GetPointData()->AddArray(Density.getVTKField("Density", Sim->particleList.size() * imageCounter / (NBins * NBins * NBins * getNumberDensity())));

  //Add the sample counts for each cell
  image->GetPointData()->AddArray(SampleCounter.getVTKField("Samples per Cell"));

  //Average the velocity and output
  for (long z = 0; z < NBins; z++)
    for (long y = 0; y < NBins; y++)
      for (long x = 0; x < NBins; x++)
	Velocity[x][y][z] /= SampleCounter[x][y][z];
  image->GetPointData()->AddArray(Velocity.getVTKField("Velocity Field"));

  //Averaged velocity without the X component
  CFieldArray<Vector  > VelocityNoX(Velocity);
  for (long z = 0; z < NBins; z++)
    for (long y = 0; y < NBins; y++)
      for (long x = 0; x < NBins; x++)
	VelocityNoX[x][y][z][0] = 0.0;
  image->GetPointData()->AddArray(VelocityNoX.getVTKField("Velocity with no x component"));


  //Average the velocity over the x direction
  CFieldArray<Vector  > VelocityYZPlane(Sim);
  for (long z = 0; z < NBins; z++)
    for (long y = 0; y < NBins; y++)
      {
	for (long x = 0; x < NBins; x++)
	  VelocityYZPlane[((long)0)][y][z] += Velocity[x][y][z];
	
	VelocityYZPlane[((long)0)][y][z] /= NBins;
	
	//Now copy the plane to fill the array
	for (long x = 1; x < NBins; x++)
	  VelocityYZPlane[x][y][z] = VelocityYZPlane[((long)0)][y][z];
      }
  image->GetPointData()->AddArray(VelocityYZPlane.getVTKField("Velocity Field avg. over the x direction"));
  
  //Velocity avgd over the x dimension with no X component
  CFieldArray<Vector  > VelocityYZPlaneNoX(VelocityYZPlane);
  for (long z = 0; z < NBins; z++)
    for (long y = 0; y < NBins; y++)
      for (long x = 0; x < NBins; x++)
	VelocityYZPlaneNoX[x][y][z][0] = 0.0;
  image->GetPointData()->AddArray(VelocityYZPlaneNoX.getVTKField("Velocity Field avg. over the x direction, no X component"));
  
  //Energy field
  CFieldArray<Vector  > Energy(Sim);
  for (long z = 0; z < NBins; z++)
    for (long y = 0; y < NBins; y++)
      for (long x = 0; x < NBins; x++)
	Energy[x][y][z] = (Vsquared[x][y][z]/SampleCounter[x][y][z]) 
	  - Velocity[x][y][z].square();
  image->GetPointData()->AddArray(Energy.getVTKField("Vsquared (Energy)"));
  
  //Output this vtkImage to an xml file
  vtkXMLImageDataWriter *writer = vtkXMLImageDataWriter::New();
  writer->SetInput(image);
  writer->SetFileName("paraview.vti");
  writer->Write();
}
