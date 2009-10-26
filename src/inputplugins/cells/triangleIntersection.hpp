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

#pragma once
#include "cell.hpp"
#include "../../datatypes/vector.hpp"
#include <string>
#include <fstream>
#include <boost/array.hpp>
#include <boost/foreach.hpp>

struct CUTriangleIntersect: public CUCell
{
  Iflt _radius;
  std::string _fileName;
  typedef boost::array<Vector, NDIM+1> triangle_type;
  std::vector<triangle_type> _triangles;

  CUTriangleIntersect(CUCell* nextCell, Iflt radius, std::string fileName):
    CUCell(nextCell),
    _radius(radius),
    _fileName(fileName)
  {}

  virtual void initialise() 
  { 
    uc->initialise();
    std::ifstream input(_fileName.c_str());
    if (!input)
      D_throw() << "Could not open " << _fileName << " to load the triangles from";

    
    while (!input.eof())
      {
	triangle_type tmp;
	for (size_t i(0); i < NDIM; ++i)
	  for (size_t j(0); j < NDIM; ++j)
	    input >> tmp[i][j];
	_triangles.push_back(tmp);
      }

    std::cout << "\nCUTriangleIntersect :Loaded " << _triangles.size() << " triangles";

    BOOST_FOREACH(triangle_type& triangle, _triangles)
      {
	triangle[3] = ((triangle[1]-triangle[0]) ^ (triangle[2]-triangle[0]));
	triangle[3] /= triangle[3].nrm();
      }
      
  }

  virtual std::vector<Vector  > placeObjects(const Vector & centre)
  {
    std::vector<Vector> initval(uc->placeObjects(centre)), retval;
    
    
    return retval;
  }
  
  bool triangleIntersects(const Vector& sphere, const triangle_type& triangle) const
  {
    Iflt p = 1;
  }
};
