/*  DYNAMO:- Event driven molecular dynamics simulator 
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
#ifdef DYNAMO_VTK
#include "field_array.hpp"

#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkImageData.h>
#include <vtkRectilinearGrid.h>
#include <vtkPointData.h>
#include <vtkXMLImageDataWriter.h>

namespace DYNAMO
{
  inline vtkImageData* getVTKImage(const DYNAMO::SimData* const Sim)
  {
    vtkImageData *vol = vtkImageData::New();
    vol->SetDimensions(NBins,NBins,NBins);
    vol->SetOrigin(-0.5*Sim->primaryCellSize[0],-0.5*Sim->primaryCellSize[1],-0.5*Sim->primaryCellSize[2]);
    vol->SetSpacing(Sim->primaryCellSize[0]/NBins, Sim->primaryCellSize[1]/NBins, Sim->primaryCellSize[2]/NBins);
    
    return vol;
  }
  
  inline vtkRectilinearGrid* getVTKRectilinearGrid(const DYNAMO::SimData* const Sim)
  {
    vtkRectilinearGrid *rgrid = vtkRectilinearGrid::New();
    vtkFloatArray *xCoords = vtkFloatArray::New();
    vtkFloatArray *yCoords = vtkFloatArray::New();
    vtkFloatArray *zCoords = vtkFloatArray::New();
    for (int i=0; i < NBins; i++) xCoords->InsertNextValue(i*(Sim->primaryCellSize[0])/NBins - 0.5*(Sim->primaryCellSize[0]));  
    for (int i=0; i < NBins; i++) yCoords->InsertNextValue(i*Sim->primaryCellSize[1]/NBins - 0.5*Sim->primaryCellSize[1]);  
    for (int i=0; i < NBins; i++) zCoords->InsertNextValue(i*Sim->primaryCellSize[2]/NBins - 0.5*Sim->primaryCellSize[2]);
    
    rgrid->SetDimensions(NBins,NBins,NBins);
    rgrid->SetXCoordinates(xCoords);
    rgrid->SetYCoordinates(yCoords);
    rgrid->SetZCoordinates(zCoords);
    
    return rgrid;
  }

  template<class T> 
  vtkFloatArray* getVTKFloatField(CFieldArray<T>& FA, const char *fieldName, double scale = 1.0)
  {
    vtkFloatArray *scalars = vtkFloatArray::New();
    scalars->SetName(fieldName);
    
    for (long z = 0; z < NBins; z++)
      for (long y = 0; y < NBins; y++)
	for (long x = 0; x < NBins; x++)
	  {
	    long offset = z*NBins*NBins 
	      + y*NBins + x;
	    scalars->InsertTuple1(offset,FA.Field[x][y][z]/scale);
	  }
    
    return scalars;
  }

  template<class T>
  inline vtkIntArray* getVTKIntField(CFieldArray<T>& FA,const char *fieldName, double scale = 1.0)
  {
    vtkIntArray *scalars = vtkIntArray::New();
    scalars->SetName(fieldName);
    
    for (long z = 0; z < NBins; z++)
      for (long y = 0; y < NBins; y++)
	for (long x = 0; x < NBins; x++)
	  {
	    long offset = z*NBins*NBins 
	      + y*NBins + x;
	    scalars->InsertTuple1(offset, FA.Field[x][y][z]/scale);
	  }
    
    return scalars;
  }

  template<class T>
  vtkFloatArray* getVTKField(CFieldArray<CVector<T> >& FA, const char *fieldName, double scale = 1.0)
  {
    vtkFloatArray *vectors = vtkFloatArray::New();
    vectors->SetNumberOfComponents(3);
    vectors->SetNumberOfTuples(NBins*NBins*NBins);
    vectors->SetName(fieldName);
    
    for (long z = 0; z < NBins; z++)
      for (long y = 0; y < NBins; y++)
	for (long x = 0; x < NBins; x++)
	  {
	    long offset = z*NBins*NBins 
	      + y*NBins + x;
	    vectors->InsertTuple(offset,(FA.Field[x][y][z]/scale).data);
	 }
    
    return vectors;
  }

}
#endif

