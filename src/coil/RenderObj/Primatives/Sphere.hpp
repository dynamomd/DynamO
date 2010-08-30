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
#pragma once
#include <vector>

class Sphere {
public:
  enum SphereType {
    tetrahedron,
    octahedron,
    icosahedron
  };
  
  Sphere(SphereType type, size_t Order);
  ~Sphere();

  int n_vertices;
  int n_faces;
  int n_edges;
  float *vertices;
  int *faces; 
  
  int edge_walk; 
  int *start; 
  int *end; 
  int *midpoint;

private:

  void refine();
 
  int search_midpoint(int index_start, int index_end);
};
