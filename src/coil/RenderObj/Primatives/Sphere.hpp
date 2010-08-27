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
