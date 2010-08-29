#pragma once
#include "Triangles.hpp"
#include <time.h>
#include "Primatives/Sphere.hpp"
#include "../clWindow.hpp"

class RTSpheres : public RTriangles
{
public:
  RTSpheres(cl::CommandQueue& CmdQ, cl::Context& Context, cl::Device& Device, bool hostTransfers,
	    const float& cameraX, const float& cameraY, const float& cameraZ,
	    size_t N, 
	    Sphere::SphereType type1, 
	    size_t order1,
	    Sphere::SphereType type2, 
	    size_t order2,
	    size_t nSphere1);

  virtual void clTick(cl::CommandQueue& CmdQ, cl::Context& Context);

protected:
  cl::Kernel _renderKernel;
  cl::Kernel _sortDataKernel;
  cl::Kernel _sortKernel;

  cl_uint _N;
  cl_uint _numStages, _powerOfTwo;

  static const std::string kernelsrc;

  Sphere primSphere1;
  Sphere primSphere2;

  cl::Buffer _spherePositions;
  cl::Buffer _primativeVertices1;
  cl::Buffer _primativeVertices2;
  cl::Buffer _sortData;
  size_t _workgroupsize;
  size_t _globalsize;
  cl_uint _nSpheres1;

  const float& _cameraX,&_cameraY,&_cameraZ;
};

template<>
inline void CLGLWindow::addRenderObj<RTSpheres, size_t, Sphere::SphereType, size_t, Sphere::SphereType, size_t, size_t>
(size_t N, Sphere::SphereType type1, size_t order1, Sphere::SphereType type2, size_t order2, size_t nSphere1)
{
  RenderObjects.push_back(new RTSpheres(_clcmdq, _clcontext, _cldevice, _hostTransfers,
					_cameraX, _cameraY, _cameraZ,
					N, type1, order1, type2, order2, nSphere1));
}
