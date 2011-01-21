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
#include <vector>
#include <cmath>

#include <string.h>
#include <stdlib.h>
#include <coil/extcode/vector2.hpp>
#include <magnet/exception.hpp>

namespace magnet {
  namespace GL {
    namespace primatives {
      class Sphere {
      public:
	enum SphereType {
	  tetrahedron,
	  octahedron,
	  icosahedron
	};
	
	inline Sphere(SphereType type, size_t order):
	  vertices(NULL),
	  faces(NULL)
	{
	  switch (type)
	    {
	    case tetrahedron:
	      { 
		float sqrt3 = 1 / std::sqrt(3.0);
		float tetrahedron_vertices[] = {sqrt3, sqrt3, sqrt3,
						-sqrt3, -sqrt3, sqrt3,
						-sqrt3, sqrt3, -sqrt3,
						sqrt3, -sqrt3, -sqrt3}; 
		int tetrahedron_faces[] = {0, 2, 1, 0, 1, 3, 2, 3, 1, 3, 2, 0};
	
		n_vertices = 4; 
		n_faces = 4; 
		n_edges = 6; 
		vertices = (float*)malloc(3*n_vertices*sizeof(float)); 
		faces = (int*)malloc(3*n_faces*sizeof(int)); 
		memcpy ((void*)vertices, (void*)tetrahedron_vertices, 3*n_vertices*sizeof(float)); 
		memcpy ((void*)faces, (void*)tetrahedron_faces, 3*n_faces*sizeof(int)); 
	      } 
	      break;
	    case octahedron:
	      { 
		float octahedron_vertices[] = {0.0, 0.0, -1.0,
					       1.0, 0.0, 0.0,
					       0.0, -1.0, 0.0,
					       -1.0, 0.0, 0.0,
					       0.0, 1.0, 0.0,
					       0.0, 0.0, 1.0}; 
		int octahedron_faces[] = {0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 1, 5, 2, 1, 5, 3, 2, 5, 4, 3, 5, 1, 4}; 
	
		n_vertices = 6; 
		n_faces = 8;
		n_edges = 12; 
		vertices = (float*)malloc(3*n_vertices*sizeof(float)); 
		faces = (int*)malloc(3*n_faces*sizeof(int)); 
		memcpy ((void*)vertices, (void*)octahedron_vertices, 3*n_vertices*sizeof(float)); 
		memcpy ((void*)faces, (void*)octahedron_faces, 3*n_faces*sizeof(int)); 
	      } 
	      break;
	    case icosahedron:
	      { 
		float t = (1+sqrt(5))/2;
		float tau = t/sqrt(1+t*t);
		float one = 1/sqrt(1+t*t);
	
		float icosahedron_vertices[] = {tau, one, 0.0,
						-tau, one, 0.0,
						-tau, -one, 0.0,
						tau, -one, 0.0,
						one, 0.0 ,  tau,
						one, 0.0 , -tau,
						-one, 0.0 , -tau,
						-one, 0.0 , tau,
						0.0 , tau, one,
						0.0 , -tau, one,
						0.0 , -tau, -one,
						0.0 , tau, -one};

		int icosahedron_faces[] = {4, 8, 7,
					   4, 7, 9,
					   5, 6, 11,
					   5, 10, 6,
					   0, 4, 3,
					   0, 3, 5,
					   2, 7, 1,
					   2, 1, 6,
					   8, 0, 11,
					   8, 11, 1,
					   9, 10, 3,
					   9, 2, 10,
					   8, 4, 0,
					   11, 0, 5,
					   4, 9, 3,
					   5, 3, 10,
					   7, 8, 1,
					   6, 1, 11,
					   7, 2, 9,
					   6, 10, 2};
	
		n_vertices = 12; 
		n_faces = 20;
		n_edges = 30;
		vertices = (float*)malloc(3*n_vertices*sizeof(float)); 
		faces = (int*)malloc(3*n_faces*sizeof(int)); 
		memcpy ((void*)vertices, (void*)icosahedron_vertices, 3*n_vertices*sizeof(float)); 
		memcpy ((void*)faces, (void*)icosahedron_faces, 3*n_faces*sizeof(int)); 
	      }
	      break;
	    default:
	      M_throw() << "Unknown Sphere Type specified";
	    }

	  for (size_t i(0); i < order; ++i)
	    refine();

  
	  //Renormalize to the volume of a sphere with radius 1!
	  double spherevol = 4.0 * M_PI / 3.0;

	  //Now calculate the current volume
	  double volumeSum = 0;

	  //Iterate over each surface triangle
	  for (int i(0); i < n_faces; ++i)
	    {
	      Vector a(vertices[3 * faces[3*i + 0] + 0],
		       vertices[3 * faces[3*i + 0] + 1],
		       vertices[3 * faces[3*i + 0] + 2]);

	      Vector b(vertices[3 * faces[3*i + 1] + 0],
		       vertices[3 * faces[3*i + 1] + 1],
		       vertices[3 * faces[3*i + 1] + 2]);

	      Vector c(vertices[3 * faces[3*i + 2] + 0],
		       vertices[3 * faces[3*i + 2] + 1],
		       vertices[3 * faces[3*i + 2] + 2]);

	      volumeSum += (a | (b ^ c));
	    }
  
	  //Calculate the ratio of the two volumes
	  volumeSum /= 6;

	  double lengthscale = std::pow(spherevol / volumeSum, 1.0 / 3.0);

	  for (int i(0); i < n_vertices; ++i)
	    for (size_t j(0); j < 3; ++j)
	      vertices[i * 3 + j] *= lengthscale;

	}

	inline ~Sphere()
	{
	  free (vertices);
	  free (faces);
	}

	inline Sphere(const Sphere& other)
	{
	  n_vertices = other.n_vertices; 
	  n_faces = other.n_faces;
	  n_edges = other.n_edges;
	  vertices = (float*)malloc(3*n_vertices*sizeof(float)); 
	  faces = (int*)malloc(3*n_faces*sizeof(int));
	  memcpy ((void*)vertices, (void*)other.vertices, 3*n_vertices*sizeof(float)); 
	  memcpy ((void*)faces, (void*)other.faces, 3*n_faces*sizeof(int)); 
	}

	inline const int& getVertexCount() const {return n_vertices; }
	inline const int& getFaceCount() const {return n_faces; }
	inline const int& getEdgeCount() const {return n_edges; }
	inline float* getVertices() {return vertices; }
	inline const float* getVertices() const {return vertices; }
	inline const int* getFaces() const {return faces; }
	
      private:

	int n_vertices;
	int n_faces;
	int n_edges;

	float *vertices;
	int *faces; 
	
	int edge_walk; 
	
	int* start;
	int* midpoint;
	int* end;

	inline void refine()
	{
	  int n_vertices_new = n_vertices+2*n_edges; 
	  int n_faces_new = 4*n_faces; 
	  int i; 
  
	  edge_walk = 0; 
	  n_edges = 2*n_vertices + 3*n_faces; 
	  start = (int*)malloc(n_edges*sizeof (int)); 
	  end = (int*)malloc(n_edges*sizeof (int)); 
	  midpoint = (int*)malloc(n_edges*sizeof (int)); 
  
	  int *faces_old = (int*)malloc (3*n_faces*sizeof(int)); 
	  faces_old = (int*)memcpy((void*)faces_old, (void*)faces, 3*n_faces*sizeof(int)); 
	  vertices = (float*)realloc ((void*)vertices, 3*n_vertices_new*sizeof(float)); 
	  faces = (int*)realloc ((void*)faces, 3*n_faces_new*sizeof(int)); 
	  n_faces_new = 0; 
  
	  for (i=0; i<n_faces; i++) 
	    { 
	      int a = faces_old[3*i]; 
	      int b = faces_old[3*i+1]; 
	      int c = faces_old[3*i+2]; 
      
	      int ab_midpoint = search_midpoint (b, a); 
	      int bc_midpoint = search_midpoint (c, b); 
	      int ca_midpoint = search_midpoint (a, c); 
      
	      faces[3*n_faces_new] = a; 
	      faces[3*n_faces_new+1] = ab_midpoint; 
	      faces[3*n_faces_new+2] = ca_midpoint; 
	      n_faces_new++; 
	      faces[3*n_faces_new] = ca_midpoint; 
	      faces[3*n_faces_new+1] = ab_midpoint; 
	      faces[3*n_faces_new+2] = bc_midpoint; 
	      n_faces_new++; 
	      faces[3*n_faces_new] = ca_midpoint; 
	      faces[3*n_faces_new+1] = bc_midpoint; 
	      faces[3*n_faces_new+2] = c; 
	      n_faces_new++; 
	      faces[3*n_faces_new] = ab_midpoint; 
	      faces[3*n_faces_new+1] = b; 
	      faces[3*n_faces_new+2] = bc_midpoint; 
	      n_faces_new++; 
	    } 
	  n_faces = n_faces_new; 
	  free (start); 
	  free (end); 
	  free (midpoint); 
	  free (faces_old); 
	}

	
	inline int search_midpoint (int index_start, int index_end)
	{ 
	  int i;
	  for (i=0; i<edge_walk; i++) 
	    if ((start[i] == index_start && end[i] == index_end) || 
		(start[i] == index_end && end[i] == index_start)) 
	      {
		int res = midpoint[i];

		/* update the arrays */
		start[i]    = start[edge_walk-1];
		end[i]      = end[edge_walk-1];
		midpoint[i] = midpoint[edge_walk-1];
		edge_walk--;
	
		return res; 
	      }

	  /* vertex not in the list, so we add it */
	  start[edge_walk] = index_start;
	  end[edge_walk] = index_end; 
	  midpoint[edge_walk] = n_vertices; 
  
	  /* create new vertex */ 
	  vertices[3*n_vertices]   = (vertices[3*index_start] + vertices[3*index_end]) / 2.0;
	  vertices[3*n_vertices+1] = (vertices[3*index_start+1] + vertices[3*index_end+1]) / 2.0;
	  vertices[3*n_vertices+2] = (vertices[3*index_start+2] + vertices[3*index_end+2]) / 2.0;
  
	  /* normalize the new vertex */ 
	  float length = sqrt (vertices[3*n_vertices] * vertices[3*n_vertices] +
			       vertices[3*n_vertices+1] * vertices[3*n_vertices+1] +
			       vertices[3*n_vertices+2] * vertices[3*n_vertices+2]);
	  length = 1/length;
	  vertices[3*n_vertices] *= length;
	  vertices[3*n_vertices+1] *= length;
	  vertices[3*n_vertices+2] *= length;
  
	  n_vertices++;
	  edge_walk++;
	  return midpoint[edge_walk-1];
	} 
      };
    }
  }
}
