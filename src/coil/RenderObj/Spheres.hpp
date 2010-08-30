#pragma once
#include "Triangles.hpp"
#include <time.h>
#include "Primatives/Sphere.hpp"
#include "../clWindow.hpp"

class RTSpheres : public RTriangles
{
public:
  struct SphereDetails
  {
    inline SphereDetails(Sphere::SphereType type, size_t order, size_t n):
      _type(type, order),
      _nSpheres(n)
    {}

    Sphere _type;
    cl_uint _nSpheres;

    void setupCLBuffers(cl::Context& Context)
    {
      _primativeVertices = cl::Buffer(Context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
				      sizeof(cl_float) * 3 * _type.n_vertices,
				      _type.vertices);
    }

    cl::Buffer _primativeVertices;
  };

  RTSpheres(cl::CommandQueue& CmdQ, cl::Context& Context, cl::Device& Device, bool hostTransfers,
	    const float& cameraX, const float& cameraY, const float& cameraZ,
	    size_t N, const std::vector<SphereDetails>& renderDetailLevels);

  virtual void clTick(cl::CommandQueue& CmdQ, cl::Context& Context);

protected:
  cl::Kernel _renderKernel;
  cl::Kernel _sortDataKernel;
  cl::Kernel _sortKernel;

  cl_uint _N;
  cl_uint _numStages, _powerOfTwo;

  static const std::string kernelsrc;

  std::vector<SphereDetails> _renderDetailLevels;

  cl::Buffer _spherePositions;
  cl::Buffer _sortData;
  size_t _workgroupsize;
  size_t _globalsize;

  const float& _cameraX,&_cameraY,&_cameraZ;
};

template<>
inline void CLGLWindow::addRenderObj<RTSpheres, size_t, std::vector<RTSpheres::SphereDetails> >
(size_t N, std::vector<RTSpheres::SphereDetails> renderDetailLevels)
{
  RenderObjects.push_back(new RTSpheres(_clcmdq, _clcontext, _cldevice, _hostTransfers,
					_cameraX, _cameraY, _cameraZ,
					N, renderDetailLevels));
}
