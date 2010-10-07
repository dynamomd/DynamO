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
#include "../clWindow.hpp"
#include "Triangles.hpp"
#include <time.h>
#include <magnet/GL/primatives/Sphere.hpp>
#include <magnet/thread/mutex.hpp>
#include <magnet/CL/radixsort.hpp>
#include <magnet/CL/heapSort.hpp>

class RTSpheres : public RTriangles
{
public:
  struct SphereDetails
  {
    inline SphereDetails(magnet::GL::primatives::Sphere::SphereType type, size_t order, size_t n):
      _type(type, order),
      _nSpheres(n)
    {}

    magnet::GL::primatives::Sphere _type;
    cl_uint _nSpheres;

    void setupCLBuffers(magnet::CL::CLGLState& CLState)
    {
      _primativeVertices = cl::Buffer(CLState.getContext(), CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
				      sizeof(cl_float) * 3 * _type.getVertexCount(),
				      _type.getVertices());
    }

    cl::Buffer _primativeVertices;
  };

  RTSpheres(const magnet::GL::viewPort& viewPortInfo,
	    size_t N, const std::vector<SphereDetails>& renderDetailLevels);

  ~RTSpheres();

  cl_float4* writePositionData(magnet::CL::CLGLState&);
  void returnPositionData(magnet::CL::CLGLState&, cl_float4* );

  virtual void clTick(magnet::CL::CLGLState&);
  void sortTick(magnet::CL::CLGLState&);

  virtual void initOpenGL() {}
  virtual void initOpenCL(magnet::CL::CLGLState&);
  
protected:

  void clTick_no_sort_or_locking(magnet::CL::CLGLState&);

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


  magnet::CL::radixSort<cl_float> sortFunctor;
  magnet::CL::heapSort<cl_float> CPUsortFunctor;

  const magnet::GL::viewPort & _viewPortInfo;

  magnet::thread::Mutex _sphereDataLock;
};

template<>
inline RTSpheres& CLGLWindow::addRenderObj<RTSpheres, size_t, std::vector<RTSpheres::SphereDetails> >
(size_t N, std::vector<RTSpheres::SphereDetails> renderDetailLevels)
{
  RenderObjects.push_back(new RTSpheres(_viewPortInfo, N, renderDetailLevels));
  return *((RTSpheres*)RenderObjects.back());
}
