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

#pragma once
#include "cell.hpp"
#include "../../datatypes/vector.hpp"
#include <string>
#include <fstream>
#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include <boost/progress.hpp>

struct CUTriangleIntersect: public CUCell
{
  double _diameter, _diametersq;
  std::string _fileName;
  typedef boost::array<Vector, NDIM+1> triangle_type;
  std::vector<triangle_type> _triangles;

  CUTriangleIntersect(CUCell* nextCell, double diameter, std::string fileName):
    CUCell(nextCell),
    _diameter(diameter),
    _diametersq(diameter * diameter),
    _fileName(fileName)
  {}

  virtual void initialise() 
  { 
    uc->initialise();
    std::ifstream input(_fileName.c_str());
    if (!input)
      M_throw() << "Could not open " << _fileName << " to load the triangles from";

    
    while (!input.eof())
      {
	triangle_type tmp;
	for (size_t i(0); i < NDIM; ++i)
	  for (size_t j(0); j < NDIM; ++j)
	    input >> tmp[i][j];
	
	_triangles.push_back(tmp);
      }

    std::cout << "\nCUTriangleIntersect :Loaded " << _triangles.size() << " triangles";
    std::cout.flush();
    
    BOOST_FOREACH(triangle_type& triangle, _triangles)
      {
	//std::cout << "\nCUTriangleIntersect :Loaded " << triangle[0] << " " << triangle[1] << " " << triangle[2] <<" triangle";
	
	//triangle[3] is the triangle norm
	triangle[3] = ((triangle[1]-triangle[0]) ^ (triangle[2]-triangle[0]));
	triangle[3] /= triangle[3].nrm();	
	//The other two vectors are set to an origin at triangle[0]
	triangle[1] -= triangle[0];
	triangle[2] -= triangle[0];
      }
      
  }

  virtual std::vector<Vector  > placeObjects(const Vector & centre)
  {
    std::vector<Vector>  retval, initval(uc->placeObjects(centre));
    
    std::cout << "\nCUTriangleIntersect :Checking spheres\n";
    std::cout.flush();
    
    boost::progress_display prog(initval.size());


    BOOST_FOREACH(const Vector& sphere, initval)
      {
	BOOST_FOREACH(const triangle_type& triangle, _triangles)
	  if (triangleIntersects(sphere, triangle))
	    {
	      retval.push_back(sphere);
	      break;
	    }
	++prog;
      }
    
    return retval;
  }
  
  bool triangleIntersects(const Vector& sphere, const triangle_type& triangle) const
  {
    
    //Check the plane of the triangle and sphere intersect
    //Elevation of the sphere over the triangle
    double p = (sphere - triangle[0]) | triangle[3];

    if (fabs(p) > 0.5 * _diameter)
      return false;

    //Check if any of the verticies of the triangle are in the sphere
    if (_diametersq >= (triangle[0] - sphere).nrm2())
      return true;

    if (_diametersq >= (triangle[0] + triangle[1] - sphere).nrm2())
      return true;

    if (_diametersq >= (triangle[0] + triangle[2] - sphere).nrm2())
      return true;
    
    //Check if the triangle contains the sphere
    {
      //Generate the point of the sphere on the surface of the triangle
      Vector C = sphere - triangle[3] * p - triangle[0];
      
      //std::cerr << "\nCross product should be zero and is actually " << ((C ^ triangle[1]) | triangle[0]); 
            
      //!Taken from http://www.blackpawn.com/texts/pointinpoly/default.html
      //v0 = triangle[1], v1 = triangle[2], v2 = C
      double dot00 = triangle[1].nrm2();
      double dot11 = triangle[2].nrm2();
      double dot01 = triangle[1] | triangle[2];
      double dot02 = triangle[1] | C;
      double dot12 = triangle[2] | C;
    
      // Compute barycentric coordinates
      double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
      double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
      double v = (dot00 * dot12 - dot01 * dot02) * invDenom;
    	
    	// Check if point is in triangle
      if ((u > 0) && (v > 0) && (u + v < 1)) 
    	return true;
    }

    //Check if the edges of the triangle intersect the spheres
    //Just check for real roots of the sphere sphere intersection test 
    if (sphereEdgeCheck(triangle[0], triangle[1], sphere)) return true;
    if (sphereEdgeCheck(triangle[0], triangle[2], sphere)) return true;
    if (sphereEdgeCheck(triangle[0] + triangle[1], triangle[2] - triangle[1], sphere)) return true;

    return false;
  }

  bool sphereEdgeCheck(const Vector& linecentre, const Vector& edge, const Vector& sphere) const
  {
    //Taken from http://jgt.akpeters.com/papers/KarabassiEtAl99/collision.html
    //n = edge,  
    double a = edge | linecentre;
    double b = edge | sphere;
    double c = edge | edge;

    double t = (b - a) / c;

    Vector q = linecentre + t * edge;
    
    Vector ijk = sphere - q;
    double d = ijk.nrm2();
    if (d > _diametersq) return false;

    if (d == _diametersq)
      {
	if ((t >= 0) && (t <= 1))
	  return true;
	else
	  return false;
      }

    double r = std::sqrt((_diametersq-d) / c);
    double t1 = t - r;
    double t2 = t + r;

    if ( ( (t1<0) && (t2<0) )
	 ||( (t1>1) && (t2>1) ) )
      return false;

    return true;

//    const Vector r0(linecentre - sphere);
//    double a = edge.nrm2();
//    double b = (r0 | edge);
//    double c = r0.nrm2() - (_diametersq * 0.25);
//    double arg = b * b - a * c;
//    //They are approaching? Previous tests ensured the vertex is not in the sphere already
//    if (b < 0)
//	if (arg >= 0)
//	  {
//	    double x = c / (b - std::sqrt(arg));
//	    if ((x > 0) && (x < 1)) return true;
//	  }
//
//    return false;
  }
};
