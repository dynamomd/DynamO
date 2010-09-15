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
#include "Triangles.hpp"
#include <time.h>
#include "Primatives/Sphere.hpp"
#include "../clWindow.hpp"
#include <pthread.h>

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

  RTSpheres(const magnet::GL::viewPort& viewPortInfo,
	    size_t N, const std::vector<SphereDetails>& renderDetailLevels);

  ~RTSpheres();

  cl_float4* writePositionData(cl::CommandQueue& cmdq);
  void returnPositionData(cl::CommandQueue& cmdq, cl_float4* clBufPointer);

  virtual void clTick(cl::CommandQueue& CmdQ, cl::Context& Context);
  void sortTick(cl::CommandQueue& CmdQ, cl::Context& Context);

  virtual void initOpenGL() {}
  virtual void initOpenCL(cl::CommandQueue& CmdQ, cl::Context& Context, cl::Device& Device, bool hostTransfers);
  
protected:

  void clTick_no_sort_or_locking(cl::CommandQueue& CmdQ, cl::Context& Context);

  cl::Kernel _renderKernel;
  cl::Kernel _sortDataKernel;

  cl_uint _N;
  static const std::string kernelsrc;

  std::vector<SphereDetails> _renderDetailLevels;

  cl::Buffer _spherePositions;
  cl::Buffer _sortKeys, _sortData;

  size_t _frameCount;
  size_t _sortFrequency;
  size_t _workgroupsize;
  size_t _globalsize;

  const magnet::GL::viewPort & _viewPortInfo;

  pthread_mutex_t _sphereDataLock;
};

template<>
inline RTSpheres& CLGLWindow::addRenderObj<RTSpheres, size_t, std::vector<RTSpheres::SphereDetails> >
(size_t N, std::vector<RTSpheres::SphereDetails> renderDetailLevels)
{
  RenderObjects.push_back(new RTSpheres(_viewPortInfo, N, renderDetailLevels));
  return *((RTSpheres*)RenderObjects.back());
}
