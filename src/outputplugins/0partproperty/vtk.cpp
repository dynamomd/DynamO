/*  DYNAMO:- Event driven molecular dynamics simulator 
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2008  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

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

#ifdef DYNAMO_VTK
#include "vtk.hpp"
#include <fstream>
#include <boost/foreach.hpp>
#include "../../extcode/xmlwriter.hpp"
#include "../../dynamics/include.hpp"
#include "../../base/is_exception.hpp"
#include "../../base/is_simdata.hpp"
#include "../../dynamics/liouvillean/liouvillean.hpp"
#include "../../datatypes/field_array_vtk.hpp"

COPVTK::COPVTK(const DYNAMO::SimData* tmp, const XMLNode&):
  COPCollTicker(tmp,"VTK"),
  Density(tmp), 
  Vsquared(tmp),
  SampleCounter(tmp),
  Velocity(tmp),
  imageCounter(0)
{}

void 
COPVTK::stream(Iflt)
{}

void 
COPVTK::ticker()
{
  if (!(Sim->lNColl % 100))
    {
      imageCounter++;

      Sim->Dynamics.Liouvillean().updateAllParticles();

      BOOST_FOREACH(const CParticle & Part, Sim->vParticleList)
	{
	  CVector<> position = Part.getPosition();

	  Sim->Dynamics.BCs().setPBC(position);

	  //Samples
	  ++SampleCounter[position];

	  //Velocity Vectors
	  Velocity[position] += Part.getVelocity();
	  
	  //Density field
	  Density[position] += 1.0;

	  //Energy Field
	  Vsquared[position] += (Part.getVelocity()).square();
	}
    }
}

void 
COPVTK::output(xmlw::XmlStream&)
{
  //Create a image
  vtkImageData* image = DYNAMO::getVTKImage(Sim);

  //Add the density
  image->GetPointData()->AddArray(getVTKFloatField(Density, "Density", Sim->Dynamics.units().simVolume() * Sim->lN * imageCounter / (NBins * NBins * NBins)));

  //Add the sample counts for each cell
  image->GetPointData()->AddArray(getVTKIntField(SampleCounter,"Samples per Cell"));

  //Average the velocity and output
  for (long z = 0; z < NBins; z++)
    for (long y = 0; y < NBins; y++)
      for (long x = 0; x < NBins; x++)
	Velocity[x][y][z] /= SampleCounter[x][y][z];
  image->GetPointData()->AddArray(getVTKField(Velocity, "Velocity Field"));

  //Averaged velocity without the X component
  CFieldArray<CVector<> > VelocityNoX(Velocity);
  for (long z = 0; z < NBins; z++)
    for (long y = 0; y < NBins; y++)
      for (long x = 0; x < NBins; x++)
	VelocityNoX[x][y][z][0] = 0.0;
  image->GetPointData()->AddArray(getVTKField(VelocityNoX,"Velocity with no x component"));


  //Average the velocity over the x direction
  CFieldArray<CVector<> > VelocityYZPlane(Sim);
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
  image->GetPointData()->AddArray(getVTKField(VelocityYZPlane, "Velocity Field avg. over the x direction"));
  
  //Velocity avgd over the x dimension with no X component
  CFieldArray<CVector<> > VelocityYZPlaneNoX(VelocityYZPlane);
  for (long z = 0; z < NBins; z++)
    for (long y = 0; y < NBins; y++)
      for (long x = 0; x < NBins; x++)
	VelocityYZPlaneNoX[x][y][z][0] = 0.0;
  image->GetPointData()->AddArray(getVTKField(VelocityYZPlaneNoX, "Velocity Field avg. over the x direction, no X component"));
  
  //Energy field
  CFieldArray<Iflt > Energy(Sim);
  for (long z = 0; z < NBins; z++)
    for (long y = 0; y < NBins; y++)
      for (long x = 0; x < NBins; x++)
	Energy[x][y][z] = (Vsquared[x][y][z] / SampleCounter[x][y][z]) 
	  - Velocity[x][y][z].square();
  image->GetPointData()->AddArray(getVTKFloatField(Energy,"Vsquared (Energy)"));
  
  //Output this vtkImage to an xml file
  vtkXMLImageDataWriter *writer = vtkXMLImageDataWriter::New();
  writer->SetInput(image);
  writer->SetFileName("paraview.vti");
  writer->Write();
}
#endif
